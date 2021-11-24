#load packages
import argparse
import sys
import linecache
import pandas as pd
import os
import re
import zipfile
import requests
import logging
from Bio.pairwise2 import align
from Bio.pairwise2 import format_alignment
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import warnings
warnings.filterwarnings('ignore')
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
# supporting functions
def translate(seq, strand, frame):  # translate DNA sequence according to strand and frame
    coding_dna = seq.replace("-", "")
    if strand == "-":
        coding_dna = str(Seq(coding_dna).reverse_complement())  # revcom
        coding_dna = adjust_frame(coding_dna, frame)  # adjust frame
        trans = str(Seq(coding_dna).translate())  # translate
        return (trans)
    else:
        coding_dna = adjust_frame(coding_dna, frame)  # adjust frame
        trans = str(Seq(coding_dna).translate())  # translate
        return (trans)


def adjust_frame(seq, frame):  # adjusts the frame
    if frame == 1:
        return (seq)
    elif frame == 2:
        return (seq[1:])
    elif frame == 3:
        return (seq[2:])
    else:
        return ("invalid frame!")


def translate2(seq, start, end, strand, frame):  # adjusts the frame
    coding_dna = seq[start:end].replace("-", "")
    if strand == "-":
        coding_dna = str(Seq(coding_dna).reverse_complement())  # revcom
        coding_dna = adjust_frame(coding_dna, frame)  # adjust frame
        trans = str(Seq(coding_dna).translate())  # translate
        return (trans)
    else:
        coding_dna = adjust_frame(coding_dna, frame)  # adjust frame
        trans = str(Seq(coding_dna).translate())  # translate
        return (trans)


def infer_strand(tag_for, tag_rev,
                 HDR_amp):  # use the provided tag sequence to infer the strand and frame of the coding sequence
    mFor = re.finditer(tag_for, HDR_amp)
    mRev = re.finditer(tag_rev, HDR_amp)
    if sum(1 for _ in mFor) == 1:  # found tag in sense strand
        for m in re.finditer(tag_for, HDR_amp):
            st, en = m.span()
            frame = st % 3 + 1
            # print(f"{row['gene_name']} {st} {en} + {frame}")
            return (["+", frame, st, en])
    else:
        for m in re.finditer(tag_rev, HDR_amp):
            st, en = m.span()
            frame = (len(HDR_amp) - en) % 3 + 1
            # print(f"{row['gene_name']} {st} {en} - {frame}")
            return (["-", frame, st, en])

        # get cds coordinate in amplicons


# align each cds to amplicon and mark cds segments
def get_cds_coord_in_amplicon(cds_seqs, amplicon):
    # log.info(f"Searching for cds segments in amplicon...")
    # add revcom cds
    cds_seqs_revcom = []
    for cds in cds_seqs:
        cds_seqs_revcom.append(str(Seq(cds).reverse_complement()))

    cds_coord_in_amp = []

    for cds in cds_seqs:
        aln = align.localms(amplicon, cds, 2, -1, -.5, -.1)
        # print(format_alignment(*aln[0]))
        # print((aln[0].start,aln[0].end))
        # decide if alignment is legit
        match_count = format_alignment(*aln[0]).split("\n")[1].count('|')
        total_len = len(format_alignment(*aln[0]).split("\n")[1])
        if match_count / total_len >= 0.8:
            # print(format_alignment(*aln[0]))
            # print((aln[0].start,aln[0].end))

            # get frame
            aln2 = align.localms(cds, amplicon, 2, -1, -.5, -.1)
            cds_match_start = format_alignment(*aln2[0]).split("\n")[0].split()[0]
            frame = int(cds_match_start) % 3
            if frame == 2: frame = 3
            if frame == 0: frame = 2

            cds_coord_in_amp.append((aln[0].start, aln[0].end, "+", frame))
            # log.info(f"Found a cds segment in the amplicon")

    for cds in cds_seqs_revcom:
        aln = align.localms(amplicon, cds, 2, -1, -.5, -.1)
        # print(format_alignment(*aln[0]))
        # print((aln[0].start,aln[0].end))
        # decide if alignment is legit
        match_count = format_alignment(*aln[0]).split("\n")[1].count('|')
        total_len = len(format_alignment(*aln[0]).split("\n")[1])
        if match_count / total_len >= 0.8:
            # print(format_alignment(*aln[0]))
            # print((aln[0].start,aln[0].end))

            # get frame
            aln2 = align.localms(cds, amplicon, 2, -1, -.5, -.1)
            cds_match_start = format_alignment(*aln2[0]).split("\n")[0].split()[0]
            frame = int(cds_match_start) % 3
            if frame == 2: frame = 3
            if frame == 0: frame = 2

            cds_coord_in_amp.append((aln[0].start, aln[0].end, "-", frame))
            # log.info(f"Found a cds segment (revcom) in the amplicon")

    # log.info(f"Done searching cds segments in the amplicon")
    return (cds_coord_in_amp)


# get cds seq
def get_cds_seq_in_transcript(mytranscript):
    '''
    return a list of tuples, each tuple is a set of start,end of a segment of cds
    '''
    wholeSeq = str(mytranscript.seq)
    wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
    cds_seqs = []
    for feat in mytranscript.features:
        if feat.type == "cds":
            st = feat.location.start.position
            en = feat.location.end.position
            if (feat.strand == -1):  # neg strand
                cds_seqs.append(str(Seq(wholeSeq_rc[
                                        st:en]).reverse_complement()))  # the coord are respective to the revcom of the retrieved seq (weird)
            else:  # pos strand
                cds_seqs.append(wholeSeq[st:en])
    return (cds_seqs)


def cal_gap_positions(seq):
    gap = re.compile("-+")
    # find gaps
    gaps = []
    for m in gap.finditer(seq):
        gaps.append(m.span())
    # go through all gaps and project each gap to a point in the ungapped seq
    gap_positions = []
    gap_accu_len = 0
    for gap_coord in gaps:
        gap_st = gap_coord[0]
        gap_en = gap_coord[1]
        adjusted_gap_st = gap_st - gap_accu_len
        gap_positions.append(adjusted_gap_st)
        gap_accu_len += gap_en - gap_st
    return (gap_positions)


def cal_gap_win(seq):
    gap = re.compile("-+")
    # find gaps
    gaps = []
    for m in gap.finditer(seq):
        gaps.append(m.span())
    return (gaps)


def cal_leading_gaps(seq):
    gap = re.compile("^-+")
    # find gaps
    gaps = []
    for m in gap.finditer(seq):
        gaps.append(m.span())
    return (gaps)


def cal_trailing_gaps(seq):
    gap = re.compile("-+$")
    # find gaps
    gaps = []
    for m in gap.finditer(seq):
        gaps.append(m.span())
    return (gaps)


def check_deletion_in_cds(read, amp_cds_coords):
    '''
    input: read, amp_cds_coords
    this function will first map the gap in the reads to the ref
    then check if any of the gap(deletions) are in the cds
    returns a flag indicating if deletions are found in the cds
    '''
    read_gap_loc = cal_gap_positions(read)
    # check gap in the read (deletion in the reads)
    deletion_in_cds_flag = False
    for coord_set in amp_cds_coords:
        for read_gap in read_gap_loc:
            if read_gap >= coord_set[0] and read_gap <= coord_set[1]:
                deletion_in_cds_flag = True
    return ({"deletion_in_cds_flag": deletion_in_cds_flag,
             "read_gap_loc": read_gap_loc})


def check_insertion_in_cds(ref, amp_cds_coords):
    '''
    input: ref, amp_cds_coords
    this function will first map the gap(insertions) in the ref to the original ref
    then check if any of the gap(insertions) are in the cds
    returns a flag indicating if deletions are found in the cds
    '''
    # input: ref, wt_amp_cds_coords
    ref_gap_loc = cal_gap_positions(ref)
    # check gap in the ref (insertion in the reads)
    insertion_in_cds_flag = False
    for coord_set in amp_cds_coords:
        for ref_gap in ref_gap_loc:
            if ref_gap >= coord_set[0] and ref_gap <= coord_set[1]:
                insertion_in_cds_flag = True
    return ({"insertion_in_cds_flag": insertion_in_cds_flag,
             "ref_gap_loc": ref_gap_loc, })


# fetch ensembl transcript
def fetch_ensembl_transcript(
        ensembl_transcript_id,
        exon_annot=False):
    '''
    Fetch the requested Ensembl transcript.
    Get the requested Ensembl transcript, together with exon and
    coding region (CDS) boundaries.
    Parameters
    ----------
    ensembl_transcript_id : str
    the ensembl transcript id, of the form ENST...

    Returns
    -------
    `Bio.SeqRecord`
    The requested transcript sequence, in 5' -> 3' order, together
    with exon and CDS features. The coordinates of exons and CDS
    features are relative to the sequence fragment.
    '''

    base_url = "http://rest.ensembl.org"

    # First, fetch the transcript sequence
    url = base_url + f"/sequence/id/{ensembl_transcript_id}"

    # log.info(f"Querying Ensembl for sequence of {ensembl_transcript_id}")
    response = requests.get(url, {"type": "genomic",
                                  "content-type": "application/json"})

    try:
        response.raise_for_status()
    except requests.HTTPError:
        log.error("Ensembl sequence REST query returned error "
                  "{}".format(response.text))
        raise ValueError(reponse.text)

    response_data = response.json()

    try:
        description = response_data['desc'].split(':')
        species = description[1]
        try:
            chromosome_number = int(description[2])
        except:
            chromosome_number = str(description[2])
        sequence_left = int(description[3])
        sequence_right = int(description[4])
        transcript_strand = int(description[5])

        if sequence_left > sequence_right:
            raise ValueError(f"Expected left sequence boundary {sequence_left} "
                             f"<= right sequence boundary {sequence_right}: did "
                             "the format of the Ensembl REST response change?")

        sequence_id = response_data['id']

        seq_str = response_data['seq']

        # log.info(f"Retrieved sequence {response_data['desc']} of length "
        #       f"{sequence_right - sequence_left} for species {species} on "
        #       f"strand {transcript_strand}")
    except (KeyError, ValueError) as e:
        log.error(e)
        log.error('Error parsing sequence metadata from Ensembl REST response - '
                  'did the format of the response change?')
        raise ValueError(e)

    seq = Seq(seq_str)

    record = SeqRecord(seq, id=sequence_id,
                       description=":".join(description))
    if exon_annot:

        url = base_url + f"/overlap/id/{ensembl_transcript_id}"

        # log.info(f"Querying Ensembl for overlaps of {ensembl_transcript_id}")
        response = requests.get(url, {"feature": ["cds", "exon"],
                                      "content-type": "application/json"})
        try:
            response.raise_for_status()
        except requests.HTTPError:
            log.error("Ensembl sequence REST query returned error "
                      "{}".format(response.text))
            raise ValueError(reponse.text)

        response_data = response.json()

        try:
            # Handle the unlikely event of a single piece of information overlapping a lonely transcript

            if not hasattr(response_data, '__iter__'):
                response_data = [response_data]

            for response_datum in response_data:
                if response_datum['Parent'] != ensembl_transcript_id:
                    continue

                if response_datum['assembly_name'] != species:
                    continue

                # Store feature locations 0-indexed from the left-most sequence boundary

                record.features.append(SeqFeature(
                    location=FeatureLocation(
                        int(response_datum['start']) - sequence_left,
                        int(response_datum['end']) - sequence_left + 1,
                        strand=int(response_datum['strand'])),
                    type=response_datum['feature_type']))
            num_exon_boundaries = len([f for f in record.features
                                       if f.type == 'exon'])

            num_cds_boundaries = len([f for f in record.features
                                      if f.type == 'cds'])

            # log.info(f"Retrieved {num_exon_boundaries} exons and "
            #        f"{num_cds_boundaries} coding regions for transcript "
            #        f"{ensembl_transcript_id}")
        except (KeyError, ValueError) as e:
            log.error(e)
            log.error('Error parsing overlap metadata from Ensembl REST response - '
                      'did the format of the response change?')
            raise ValueError(e)

    record.annotations['reference_species'] = species
    record.annotations['reference_chromosome_number'] = chromosome_number
    record.annotations['reference_left_index'] = sequence_left
    record.annotations['reference_right_index'] = sequence_right
    record.annotations['transcript_strand'] = transcript_strand

    # Finally, sort features by their start locations
    record.features.sort(key=lambda f: f.location.start)

    return record