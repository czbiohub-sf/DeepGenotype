# Consensus grouping for DeepGenotype
# Clusters similar reads (within edit distance threshold) and forms consensus sequences
# via majority vote, producing a new allele frequency table in CRISPResso format.

import argparse
import os
import sys
import zipfile
import logging
from collections import defaultdict

log = logging.getLogger("consensus_grouping.py")
log.setLevel(logging.INFO)


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)


def parse_args():
    parser = MyParser(
        description="Cluster similar reads and form consensus sequences for DeepGenotype"
    )
    parser.add_argument(
        "--input_zip",
        required=True,
        type=str,
        help="path to CRISPResso Alleles_frequency_table.zip",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=str,
        help="directory to write consensus allele freq table and diagnostic file",
    )
    parser.add_argument(
        "--max_edit_distance",
        default=3,
        type=int,
        help="maximum edit distance for clustering reads [default=3]",
    )
    config = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return config


def hamming_distance(s1, s2, max_dist):
    """Fast Hamming distance for same-length strings with early termination.

    Returns the number of mismatches if <= max_dist, otherwise max_dist + 1.
    Only valid when len(s1) == len(s2).
    """
    dist = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            dist += 1
            if dist > max_dist:
                return max_dist + 1
    return dist


def bounded_levenshtein(s1, s2, max_dist):
    """Compute Levenshtein distance using banded DP with early termination.

    Returns the actual distance if <= max_dist, otherwise returns max_dist + 1.
    Uses a diagonal band of width 2*max_dist+1 for O(n * max_dist) performance
    instead of O(n * m).
    """
    len1, len2 = len(s1), len(s2)

    # Quick length-difference check
    if abs(len1 - len2) > max_dist:
        return max_dist + 1

    # Same-length fast path: use Hamming distance
    if len1 == len2:
        return hamming_distance(s1, s2, max_dist)

    # Ensure s1 is the shorter string
    if len1 > len2:
        s1, s2 = s2, s1
        len1, len2 = len2, len1

    # Banded Levenshtein DP
    # Only compute cells within max_dist of the diagonal
    prev_row = [0] * (len1 + 1)
    for i in range(len1 + 1):
        prev_row[i] = i

    for j in range(1, len2 + 1):
        curr_row = [0] * (len1 + 1)
        curr_row[0] = j

        # Band limits for this row
        band_lo = max(1, j - max_dist)
        band_hi = min(len1, j + max_dist)

        # Set out-of-band cells to a large value
        if band_lo > 1:
            curr_row[band_lo - 1] = max_dist + 1

        row_min = max_dist + 1
        for i in range(band_lo, band_hi + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            curr_row[i] = min(
                curr_row[i - 1] + 1,  # insertion
                prev_row[i] + 1,  # deletion
                prev_row[i - 1] + cost,  # substitution
            )
            if curr_row[i] < row_min:
                row_min = curr_row[i]

        if row_min > max_dist:
            return max_dist + 1

        prev_row = curr_row

    return prev_row[len1] if prev_row[len1] <= max_dist else max_dist + 1


def strip_gaps(seq):
    """Remove gap characters from a sequence."""
    return seq.replace("-", "")


def parse_allele_frequency_table(input_zip_path):
    """Read all rows from the CRISPResso Alleles_frequency_table.txt inside a zip.

    Returns a list of dicts, each representing one row.
    """
    rows = []
    with zipfile.ZipFile(input_zip_path, "r") as myzip:
        with myzip.open("Alleles_frequency_table.txt") as fh:
            header = next(fh).decode().rstrip("\n").rstrip("\r")
            header_fields = header.split("\t")
            for line in fh:
                line_deco = line.decode().rstrip("\n").rstrip("\r")
                fields = line_deco.rstrip().split("\t")
                row = {
                    "Aligned_Sequence": fields[0],
                    "Reference_Sequence": fields[1],
                    "Reference_Name": fields[2],
                    "Read_Status": fields[3],
                    "n_deleted": int(fields[4]),
                    "n_inserted": int(fields[5]),
                    "n_mutated": int(fields[6]),
                    "n_Reads": int(fields[7]),
                    "perc_Reads": float(fields[8]),
                    "raw_line": line_deco,
                }
                rows.append(row)
    return rows


def separate_by_reference(rows):
    """Separate rows into HDR, non-HDR (WT), and AMBIGUOUS groups."""
    hdr_rows = []
    wt_rows = []
    ambiguous_rows = []

    for row in rows:
        ref_name = row["Reference_Name"]
        if ref_name == "HDR":
            hdr_rows.append(row)
        elif ref_name.startswith("AMBIGUOUS"):
            ambiguous_rows.append(row)
        else:
            wt_rows.append(row)

    return hdr_rows, wt_rows, ambiguous_rows


def greedy_cluster(rows, max_edit_distance):
    """Cluster rows using greedy edit-distance on gap-stripped aligned sequences.

    Rows should be pre-sorted by n_Reads descending.
    Returns a list of clusters, where each cluster is a list of rows.
    The first row in each cluster is the seed (most abundant).
    """
    if max_edit_distance <= 0:
        # No clustering: each row is its own cluster
        return [[row] for row in rows]

    # Pre-compute raw sequences and lengths
    for row in rows:
        row["_raw_seq"] = strip_gaps(row["Aligned_Sequence"])
        row["_raw_len"] = len(row["_raw_seq"])

    clusters = []  # list of (seed_row, [member_rows])
    assigned = [False] * len(rows)

    for i, row in enumerate(rows):
        if assigned[i]:
            continue

        # Start a new cluster with this row as seed
        cluster = [row]
        assigned[i] = True
        seed_raw = row["_raw_seq"]
        seed_len = row["_raw_len"]

        # Try to assign remaining unassigned rows to this cluster
        for j in range(i + 1, len(rows)):
            if assigned[j]:
                continue

            candidate = rows[j]

            # Length pre-filter: skip if length difference exceeds threshold
            if abs(candidate["_raw_len"] - seed_len) > max_edit_distance:
                continue

            dist = bounded_levenshtein(
                seed_raw, candidate["_raw_seq"], max_edit_distance
            )
            if dist <= max_edit_distance:
                cluster.append(candidate)
                assigned[j] = True

        clusters.append(cluster)

    return clusters


def form_consensus(cluster):
    """Form a consensus entry from a cluster of rows.

    If all members have same aligned-sequence length: weighted majority vote.
    Otherwise: use the seed's (first member's) aligned sequence.

    Returns a dict with the same fields as a parsed row.
    """
    seed = cluster[0]

    # Sum total reads
    total_reads = sum(row["n_Reads"] for row in cluster)

    # Check if all aligned sequences have the same length
    aligned_lengths = set(len(row["Aligned_Sequence"]) for row in cluster)

    if len(aligned_lengths) == 1:
        # All same length: do weighted majority vote
        seq_len = aligned_lengths.pop()
        consensus_chars = []
        ref_chars = list(seed["Reference_Sequence"])

        for pos in range(seq_len):
            # Count each character at this position, weighted by read count
            char_counts = defaultdict(int)
            for row in cluster:
                char_counts[row["Aligned_Sequence"][pos]] += row["n_Reads"]

            # Pick the character with highest count; break ties with seed
            max_count = max(char_counts.values())
            candidates = [c for c, cnt in char_counts.items() if cnt == max_count]

            if len(candidates) == 1:
                consensus_chars.append(candidates[0])
            else:
                # Tie: use seed's character
                seed_char = seed["Aligned_Sequence"][pos]
                if seed_char in candidates:
                    consensus_chars.append(seed_char)
                else:
                    consensus_chars.append(candidates[0])

        consensus_aligned = "".join(consensus_chars)
        consensus_ref = seed["Reference_Sequence"]
    else:
        # Different lengths: use seed's alignment
        consensus_aligned = seed["Aligned_Sequence"]
        consensus_ref = seed["Reference_Sequence"]

    # Recalculate n_deleted, n_inserted, n_mutated from the consensus alignment
    n_deleted = 0
    n_inserted = 0
    n_mutated = 0
    for c_read, c_ref in zip(consensus_aligned, consensus_ref):
        if c_read == "-" and c_ref != "-":
            n_deleted += 1
        elif c_read != "-" and c_ref == "-":
            n_inserted += 1
        elif c_read != c_ref and c_read != "-" and c_ref != "-":
            n_mutated += 1

    # Determine Read_Status
    if n_deleted == 0 and n_inserted == 0 and n_mutated == 0:
        read_status = "UNMODIFIED"
    else:
        read_status = "MODIFIED"

    return {
        "Aligned_Sequence": consensus_aligned,
        "Reference_Sequence": consensus_ref,
        "Reference_Name": seed["Reference_Name"],
        "Read_Status": read_status,
        "n_deleted": n_deleted,
        "n_inserted": n_inserted,
        "n_mutated": n_mutated,
        "n_Reads": total_reads,
        "perc_Reads": 0.0,  # will be recalculated
    }


def write_allele_frequency_table(consensus_rows, output_dir):
    """Write the consensus allele frequency table in CRISPResso format.

    Recalculates %Reads from the total reads across all consensus rows.
    Writes Alleles_frequency_table.txt then packs into .zip.
    """
    # Recalculate %Reads
    total_reads = sum(row["n_Reads"] for row in consensus_rows)
    if total_reads > 0:
        for row in consensus_rows:
            row["perc_Reads"] = row["n_Reads"] / total_reads * 100.0
    else:
        for row in consensus_rows:
            row["perc_Reads"] = 0.0

    # Sort by #Reads descending (as CRISPResso does)
    consensus_rows.sort(key=lambda r: r["n_Reads"], reverse=True)

    # Write the text file
    txt_path = os.path.join(output_dir, "Alleles_frequency_table.txt")
    with open(txt_path, "w") as fh:
        fh.write(
            "Aligned_Sequence\tReference_Sequence\tReference_Name\t"
            "Read_Status\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\n"
        )
        for row in consensus_rows:
            fh.write(
                f"{row['Aligned_Sequence']}\t"
                f"{row['Reference_Sequence']}\t"
                f"{row['Reference_Name']}\t"
                f"{row['Read_Status']}\t"
                f"{row['n_deleted']}\t"
                f"{row['n_inserted']}\t"
                f"{row['n_mutated']}\t"
                f"{row['n_Reads']}\t"
                f"{row['perc_Reads']}\n"
            )

    # Pack into zip
    zip_path = os.path.join(output_dir, "Alleles_frequency_table.zip")
    if os.path.exists(zip_path):
        os.remove(zip_path)
    with zipfile.ZipFile(zip_path, mode="x", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(txt_path, "Alleles_frequency_table.txt")

    # Remove the uncompressed text file
    os.remove(txt_path)


def write_diagnostic_file(clusters_by_group, output_dir):
    """Write a diagnostic TSV showing consensus groups and their members.

    clusters_by_group is a dict mapping group_label -> list of clusters.
    Each cluster is a list of rows (first is seed).
    """
    diag_path = os.path.join(output_dir, "consensus_groups_diagnostic.tsv")
    with open(diag_path, "w") as fh:
        fh.write(
            "Consensus_Group_ID\tReference_Name\tNum_Members\tTotal_Reads\t"
            "Seed_RawSequence\tConsensus_RawSequence\t"
            "Member_ReadCounts\n"
        )
        group_id = 1
        for group_label, clusters in clusters_by_group.items():
            for cluster_info in clusters:
                cluster = cluster_info["cluster"]
                consensus_row = cluster_info["consensus"]
                seed = cluster[0]
                seed_raw = strip_gaps(seed["Aligned_Sequence"])
                consensus_raw = strip_gaps(consensus_row["Aligned_Sequence"])
                member_counts = [str(row["n_Reads"]) for row in cluster]
                total_reads = sum(row["n_Reads"] for row in cluster)

                fh.write(
                    f"{group_id}\t"
                    f"{seed['Reference_Name']}\t"
                    f"{len(cluster)}\t"
                    f"{total_reads}\t"
                    f"{seed_raw}\t"
                    f"{consensus_raw}\t"
                    f"{','.join(member_counts)}\n"
                )
                group_id += 1


def run_consensus_grouping(input_zip, output_dir, max_edit_distance):
    """Main consensus grouping logic."""
    # Parse the allele frequency table
    log.info(f"Reading allele frequency table from {input_zip}")
    rows = parse_allele_frequency_table(input_zip)
    log.info(f"Read {len(rows)} unique allele entries")

    # Separate by reference
    hdr_rows, wt_rows, ambiguous_rows = separate_by_reference(rows)
    log.info(
        f"Separated into {len(hdr_rows)} HDR, {len(wt_rows)} WT, "
        f"{len(ambiguous_rows)} AMBIGUOUS rows"
    )

    # Track clusters for diagnostic output
    clusters_by_group = {}
    all_consensus_rows = []

    # Process each group
    for group_label, group_rows in [("HDR", hdr_rows), ("WT", wt_rows)]:
        if not group_rows:
            continue

        # Sort by read count descending
        group_rows.sort(key=lambda r: r["n_Reads"], reverse=True)

        # Cluster
        log.info(
            f"Clustering {len(group_rows)} {group_label} alleles "
            f"(max edit distance={max_edit_distance})"
        )
        clusters = greedy_cluster(group_rows, max_edit_distance)
        log.info(
            f"{group_label}: {len(group_rows)} alleles -> {len(clusters)} consensus groups"
        )

        # Form consensus for each cluster
        cluster_infos = []
        for cluster in clusters:
            consensus = form_consensus(cluster)
            all_consensus_rows.append(consensus)
            cluster_infos.append({"cluster": cluster, "consensus": consensus})

        clusters_by_group[group_label] = cluster_infos

    # Pass through AMBIGUOUS rows unchanged
    for row in ambiguous_rows:
        all_consensus_rows.append(
            {
                "Aligned_Sequence": row["Aligned_Sequence"],
                "Reference_Sequence": row["Reference_Sequence"],
                "Reference_Name": row["Reference_Name"],
                "Read_Status": row["Read_Status"],
                "n_deleted": row["n_deleted"],
                "n_inserted": row["n_inserted"],
                "n_mutated": row["n_mutated"],
                "n_Reads": row["n_Reads"],
                "perc_Reads": 0.0,
            }
        )

    # Write output
    os.makedirs(output_dir, exist_ok=True)
    log.info(f"Writing {len(all_consensus_rows)} consensus entries to {output_dir}")
    write_allele_frequency_table(all_consensus_rows, output_dir)
    write_diagnostic_file(clusters_by_group, output_dir)
    log.info("Consensus grouping complete")


def main():
    config = vars(parse_args())

    # Set up logging
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("[%(name)s] %(message)s"))
    log.addHandler(handler)

    input_zip = config["input_zip"]
    output_dir = config["output_dir"]
    max_edit_distance = config["max_edit_distance"]

    if not os.path.isfile(input_zip):
        log.error(f"Input zip file not found: {input_zip}")
        sys.exit(1)

    run_consensus_grouping(input_zip, output_dir, max_edit_distance)


if __name__ == "__main__":
    main()
