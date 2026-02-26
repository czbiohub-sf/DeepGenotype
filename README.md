# DeepGenotype
Calculates the frequencies of protein-level mutations from deep-sequencing reads of CRISPR-edited cells

&nbsp;
## Features
- Compute genotypes for protein/payload expressibility and amino acid sequence correctness (automatically detects coding regions)
- Sophisticated read error correction using consensus grouping
- Multiple modes:
  - Edit types: tagging/insertion and SNP/base-editing
  - Sequencing: short (Illumina) and long (PacBio) reads
- Batch processing a large list of samples/clones

  
&nbsp;
## Installation

**NOTE**: *if you installed DeepGenotype before 2025-01-15, please reinstall DeepGenotype to update CRISPResso2 to 2.3.1 to enable read quality-trimming.*

### Option A: Install using environment.yml (recommended)
```shell
module load anaconda # if on the hpc
git clone https://github.com/czbiohub-sf/DeepGenotype
cd DeepGenotype
conda env create -f environment.yml
conda activate DeepGenotype
pip install .
```

### Option B: Manual installation
create a conda environment and activate it
```shell
module load anaconda # if on the hpc
conda create -n DeepGenotype python=3.9
conda activate DeepGenotype
```
install CRISPResso2 and BBTools (includes bbduk.sh)
```shell
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install CRISPResso2==2.3.1 bbmap
```
clone DeepGenotype repo and install dependencies
```shell
git clone https://github.com/czbiohub-sf/DeepGenotype
cd DeepGenotype
pip install . # or pip install biopython==1.78 pandas requests openpyxl==3.1.2

```
### Verify installation
verify CRISPResso2 installation
```shell
CRISPResso -h
```
verify DeepGenotype installation
```shell
cd DeepGenotype # must be in the DeepGenotype/DeepGenotype directory
python DeepGenotype.py
```

&nbsp;
## Usage:
```shell
cd DeepGenotype # cd into the DeepGenotype/DeepGenotype directory
python DeepGenotype.py --path2csv example_csv/test.csv --path2workDir test_dir/ --path2fastqDir test_dir/fastq_dir/
```
All paths are relative to `DeepGenotype.py`, and please make sure the following two python scripts are in the same directory as DeepGenotype.py:  
`process_alleles_freq_table_INS.py`  
`process_alleles_freq_table_SNP.py`  


&nbsp;
## Read quality trimming options
DeepGenotype supports two mutually exclusive approaches for read quality filtering. Only one is active per run.

#### Option 1: BBDuk preprocessing (default)

By default, [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) (part of BBTools) runs **before** CRISPResso on each sample's fastq files. It performs adapter trimming, quality trimming, and short-read filtering in a single pass. CRISPResso's built-in fastp trimming is skipped to avoid double-trimming.

Two parameter presets are available (`--bbduk short` is the default):

| Parameter | `--bbduk short` (Illumina/MiSeq) | `--bbduk long` (PacBio) |
|---|---|---|
| Adapter reference | TruSeq, Nextera & PhiX (`ref=adapters,phix`) | TruSeq, Nextera & PhiX (`ref=adapters,phix`) |
| Trim direction | right end (`ktrim=r`) | right end (`ktrim=r`) |
| Kmer length | 27 (`k=27`) | 21 (`k=21`) |
| Max substitutions | 1 (`hdist=1`) | 1 (`hdist=1`) |
| Max substitutions+indels | 0 (`edist=0`) | 0 (`edist=0`) |
| Quality trim | both ends (`qtrim=rl`) | both ends (`qtrim=rl`) |
| Min quality | 20 (`trimq=20`) | 10 (`trimq=10`) |
| Min read length | 220 bp (`minlen=220`) | 500 bp (`minlen=500`) |

&nbsp;
#### Option 2: fastp via CRISPResso (activate with `--fastp`)

Pass `--fastp` to disable BBDuk and use CRISPResso2's built-in [fastp](https://github.com/OpenGene/fastp) integration for quality trimming instead.

Default fastp options: `--cut_front --cut_tail --cut_mean_quality 30 --cut_window_size 30`
(quality trimming from both ends, sliding window of 30 bp, mean quality threshold of 30)

```shell
python DeepGenotype.py \
--path2csv example_csv/test_INS.csv \
--path2workDir test_MiSeq_INS \
--path2fastqDir test_MiSeq_INS/fastq \
--fastp
```

You can customize fastp parameters with `--fastp_options_string`:
```shell
python DeepGenotype.py \
--path2csv example_csv/test_INS.csv \
--path2workDir test_MiSeq_INS \
--path2fastqDir test_MiSeq_INS/fastq \
--fastp \
--fastp_options_string "--cut_front --cut_tail --cut_mean_quality 25 --cut_window_size 20"
```

**Automatic retry logic**: if CRISPResso returns fewer than `--min_reads_post_filter` reads (default 50), it automatically retries with progressively less stringent fastp settings:
1. `--cut_mean_quality 30 --cut_window_size 20`
2. `--cut_mean_quality 20 --cut_window_size 10`
3. `--cut_mean_quality 20 --cut_window_size 4`

This retry logic applies **only** in fastp mode.

&nbsp;
## Consensus grouping

By default, DeepGenotype performs **consensus grouping** to correct sequencing errors before genotyping. Individual reads with minor imperfections (e.g., single-base sequencing errors) that support the same underlying biological sequence are clustered into **consensus groups**. A consensus sequence is derived for each group via weighted majority vote, and genotyping is performed on the consensus sequences.

This produces a second set of **consensus-adjusted genotype frequencies** alongside the standard results, allowing side-by-side comparison.

#### How it works

1. After CRISPResso alignment, reads are separated by reference (WT-aligned and HDR-aligned reads are clustered independently)
2. Within each reference group, reads are clustered using a greedy edit-distance algorithm (seeded from the most abundant read)
3. Reads within the edit distance threshold (default: 3) of a cluster seed are merged into that cluster
4. For each cluster, a consensus sequence is formed by weighted majority vote at each position
5. The consensus sequences are genotyped using the same protein-level classification logic as the standard pipeline

#### Output

Consensus grouping produces:
- A separate CSV file: `*_consensus_genotype_freq.csv` with consensus-adjusted genotype frequencies
- A second sheet ("Consensus Genotypes") in the output XLSX file
- Per-sample diagnostic files: `*_consensus_groups_diagnostic.tsv` showing the composition of each consensus group (member reads, read counts, consensus sequence)

#### Configuration

| Argument | Default | Description |
|---|---|---|
| `--skip_consensus` | `False` | Skip consensus grouping entirely |
| `--consensus_max_edit_distance` | `3` | Maximum edit distance for clustering reads. Higher values merge more aggressively; lower values are more conservative |

Example — disable consensus grouping:
```shell
python DeepGenotype.py \
--path2csv example_csv/test_INS.csv \
--path2workDir test_MiSeq_INS \
--path2fastqDir test_MiSeq_INS/fastq \
--skip_consensus
```

Example — use a stricter edit distance threshold:
```shell
python DeepGenotype.py \
--path2csv example_csv/test_INS.csv \
--path2workDir test_MiSeq_INS \
--path2fastqDir test_MiSeq_INS/fastq \
--consensus_max_edit_distance 1
```


&nbsp;
## Notable optional arguments
--fastq_R1_suffix &nbsp;&nbsp; (default "_R1_001.fastq.gz")  
--fastq_R2_suffix &nbsp;&nbsp; (default "_R2_001.fastq.gz")  
--single_fastq_suffix &nbsp;&nbsp; (use this option for **single-ended** reads as well as **pacbio** reads, need to specific the suffix, e.g.: fastq.gz)  
--quantification_window_size &nbsp;&nbsp; (default 50, which overrides CRISPResso2's default of 1)  
--bbduk &nbsp;&nbsp; BBDuk preprocessing mode, accepts `short` or `long` (see [Read filtering options](#read-filtering-options) below [default=short]  
--fastp &nbsp;&nbsp; use fastp (via CRISPResso) for read filtering instead of BBDuk (see [Read filtering options](#read-filtering-options) below) [default=False]  
--fastp_options_string &nbsp;&nbsp; options to pass to fastp (only used with `--fastp`), [default = '--cut_front --cut_tail --cut_mean_quality 30 --cut_window_size 30'] see [Read filtering options](#read-filtering-options)  
--skip_consensus &nbsp;&nbsp; skip the consensus grouping step [default=False] (see [Consensus grouping](#consensus-grouping) below)  
--consensus_max_edit_distance &nbsp;&nbsp; maximum edit distance for consensus read grouping [default=3] (see [Consensus grouping](#consensus-grouping) below)  


&nbsp;
## Inputs
There are *two* required input files:
- Fastq files (can be gzipped or not)
- A csv file (examples provided in `example_csv`), explanation of the columns is below

&nbsp;
## Outputs:
- A result table in the format of a csv file and a xlsx file, the table contains **sample-wise** information of:
  - Protein-level genotype frequencies
  - Two metrics that quantify DNA-level mismatches in edited alleles
    - weighted average of the percent identity of the reads (that aligned to the HDR amplicon)
    - weighted average of the number of mismatches of the reads (that aligned to the HDR amplicon)
- Consensus-adjusted genotype frequencies (unless `--skip_consensus` is used):
  - A separate csv file (`*_consensus_genotype_freq.csv`) with consensus-adjusted genotype frequencies
  - A "Consensus Genotypes" sheet in the xlsx file
  - Per-sample diagnostic files (`*_consensus_groups_diagnostic.tsv`) showing consensus group membership
- CRISPResso2 output that includes (and not limited to) the following:
  - Read aligning rate
  - Sequence-level genotype frequencies table
  - read-to-genotype assignments information


&nbsp;
## Example 1: To run Pacbio test dataset (insertion mode)
load conda, and activate the DeepGenotype conda environment
```shell
module load anaconda
conda activate DeepGenotype
```
run DeepGenotype
```shell
# in DeepGenotype/DeepGenotype directory
python DeepGenotype.py \
--path2csv example_csv/test_pacbio.csv \
--path2workDir test_PacBio \
--path2fastqDir test_PacBio/fastq \
--single_fastq_suffix .fastq
```
***NOTE***: to run DeepGenotype in the background (and thus safe to close the terminal), preprend `nohup` and append `&` to the command (or use screen, tmux, etc. instructions not listed here):
```shell
nohup python DeepGenotype.py \
--path2csv example_csv/test_pacbio.csv \
--path2workDir test_PacBio \
--path2fastqDir test_PacBio/fastq \
--single_fastq_suffix .fastq &
```
To check the terminal output (while running in the background)
```shell
cat nohup.out
```
The completed `nohup.out` should look like this
```
[DeepGenotype.py][INFO]  Genome edit type: INS
[DeepGenotype.py][INFO]  Processing sample: HEK-nocap-CLTA-R1_ccs.lbc89--lbc89.lbc89--lbc89
[DeepGenotype.py][INFO]  ...running CRISPResso
[DeepGenotype.py][INFO]  ...parsing allele frequency table and re-calculating allele frequencies
[DeepGenotype.py][INFO]  ...done
[DeepGenotype.py][INFO]  Processing sample: HEK-nocap-CLTA-R2_ccs.lbc90--lbc90.lbc90--lbc90
[DeepGenotype.py][INFO]  ...running CRISPResso
[DeepGenotype.py][INFO]  ...parsing allele frequency table and re-calculating allele frequencies
[DeepGenotype.py][INFO]  ...done
[DeepGenotype.py][INFO]  Processing sample: HEK-nocap-CLTA-R3_ccs.lbc91--lbc91.lbc91--lbc91
[DeepGenotype.py][INFO]  ...running CRISPResso
[DeepGenotype.py][INFO]  ...parsing allele frequency table and re-calculating allele frequencies
[DeepGenotype.py][INFO]  ...done
[DeepGenotype.py][INFO]  Done processing all samples in the csv file
```

&nbsp;
## Example 2: To run MiSeq test dataset (insertion mode)
load conda, and activate the DeepGenotype conda environment
```shell
module load anaconda
conda activate DeepGenotype
```
run DeepGenotype in the background (and thus safe to close the terminal)
```shell
nohup python DeepGenotype.py \
--path2csv example_csv/test_INS.csv \
--path2workDir test_MiSeq_INS \
--path2fastqDir test_MiSeq_INS/fastq &
```
To check the terminal output (while running in the background
```shell
cat nohup.out
```
The completed `nohup.out` should look like this (only first 9 lines shown)
```
[DeepGenotype.py][INFO]  Genome edit type: INS
[DeepGenotype.py][INFO]  Processing sample: mNGplate19_sorted_A2_DDX6-C
[DeepGenotype.py][INFO]  ...running CRISPResso
[DeepGenotype.py][INFO]  ...parsing allele frequency table and re-calculating allele frequencies
[DeepGenotype.py][INFO]  ...done
[DeepGenotype.py][INFO]  Processing sample: mNGplate19_sorted_A3_LSM14A-N
[DeepGenotype.py][INFO]  ...running CRISPResso
[DeepGenotype.py][INFO]  ...parsing allele frequency table and re-calculating allele frequencies
[DeepGenotype.py][INFO]  ...done
...
```

&nbsp;
## Example 3: to run MiSeq test dataset (SNP mode)
load conda, and activate the DeepGenotype conda environment
```shell
module load anaconda
conda activate DeepGenotype
```
run DeepGenotype in the background (and thus safe to close the terminal)
```shell
nohup python DeepGenotype.py \
--path2csv example_csv/test_MiSeq_SNP.csv \
--path2workDir test_MiSeq_SNP \
--path2fastqDir test_MiSeq_SNP/fastq &
```
To check the terminal output (while running in the background
```shell
cat nohup.out
```
The completed `nohup.out` should look like this (only first 9 lines shown)
```
[DeepGenotype.py][INFO]  Genome edit type: INS
[DeepGenotype.py][INFO]  Processing sample: mNGplate19_sorted_A2_DDX6-C
[DeepGenotype.py][INFO]  ...running CRISPResso
[DeepGenotype.py][INFO]  ...parsing allele frequency table and re-calculating allele frequencies
[DeepGenotype.py][INFO]  ...done
[DeepGenotype.py][INFO]  Processing sample: mNGplate19_sorted_A3_LSM14A-N
[DeepGenotype.py][INFO]  ...running CRISPResso
[DeepGenotype.py][INFO]  ...parsing allele frequency table and re-calculating allele frequencies
[DeepGenotype.py][INFO]  ...done
...
```

&nbsp;
## Explanation of columns in the input csv file
The input csv should contain the following columns with the exact names
  - Sample_ID (e.g. mNGplate19_sorted_A2_DDX6-C)  
    ***Important note***: For paired-end sequencing, only one Sample_ID is needed. We automatically find both R1 and R2 fastq files.   
     Check fastq file suffix parameters `--fastq_R1_suffix` and `--fastq_R1_suffix` in the `Usage` section.  
     For single-ended reads, set `--single_fastq_suffix` to the suffix of the fastq file.  
     Also check if you need `Fastq_extra_suffix` (below)

  - gene_name (e.g. DDX6)  
  - ENST_id (e.g. ENST00000620157)  
  - WT_amplicon_sequence
  - HDR_amplicon_sequence
  - gRNA_sequence
  - edit_type (e.g. INS or SNP, note that deletions, DEL is not supported at this point)  
      &nbsp;&nbsp;&nbsp; INS = insertion, SNP = single nucleotide polymorphism, DEL = deletion  
  - payload_block_index (e.g. 1 or 2 ...)  
      &nbsp;&nbsp;&nbsp; Default is 1. This parameter is only needed when there are multiple blocks of SNPs or insertion/deletions between the wt and HDR amplicon.  
      &nbsp;&nbsp;&nbsp; This parameter defines which block of SNPs or insertion/deletions is the **payload**  
      &nbsp;&nbsp;&nbsp; Blocks are ordered from left to right in respect to the amplicon sequence.  
      &nbsp;&nbsp;&nbsp; For example: there are 2 blocks of SNPs, the first block is a recut SNP and the second block of SNPs are of interest (payload), then `payload_block_index` should be set to 2, and the first block of SNPs will be analyzed for protein-changing mutations if it is in the coding region. 
      
  - Fastq_extra_suffix (Optional) 
     
      &nbsp;&nbsp;&nbsp;Extra suffix needed for mapping Sample_ID to corresponding fastq file names.
 
      &nbsp;&nbsp;&nbsp;**Please note that the common (as opposed to extra) suffixes are the following values by default:**  
      &nbsp;&nbsp;&nbsp;"_R1_001.fastq.gz"    
      &nbsp;&nbsp;&nbsp;"_R2_001.fastq.gz"
      
      &nbsp;&nbsp;&nbsp;For example if your fastq file names are:  
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; "mNGplate19_sorted_A2_DDX6-C_S90_R1_001.fastq.gz"  
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; "mNGplate19_sorted_A2_DDX6-C_S90_R2_001.fastq.gz"  
      &nbsp;&nbsp;&nbsp; and your sample name is "mNGplate19_sorted_A2_DDX6-C", and "_S90_" is another variable part of the name  
      &nbsp;&nbsp;&nbsp;Then you should add "_S90_" to the "Fastq_extra_suffix" column  


&nbsp;
## Helper commands

### Paginated view of fastq files
compressed
```shell
gzip -c example.fastq.gz | less
```
uncompressed
```shell
less example.fastq
```

### Count the number of lines in a fastq file (divide by 4 you'll get the read count)
compressed
```shell
gzip -c example.fastq.gz | wc -l
```
uncompressed
```shell
cat example.fastq | wc -l
```

### Count the number of line with a specific sequence in a fastq file
compressed
```shell
gzip -c example.fastq.gz | grep -c 'your-sequence-here'
```
uncompressed
```shell
cat example.fastq | grep -c 'your-sequence-here'
```
for matching sequences at the beginning of the line, add a `^`:  `'^your-sequence-here'`  
for matching sequences at the end of the line, add a `$`:  `'your-sequence-here$'`



# License
This project is licensed under the BSD 3-Clause license - see the LICENSE file for details.
