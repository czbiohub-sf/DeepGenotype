# DeepGenotype
Calculates the frequencies of protein-level mutations from deep-sequencing reads of CRISPR-edited cells

## Features
- Calculates genotypes with respect to protein/payload expressiblity and correctness  
- Automatically processes a list of samples 
- Supports CRISPR editing types: tagging/insertion and SNP/base-editing  
- Works with both Illumina and PacBio reads  
- Automatically finds coding regions  
- Invokes CRISPResso2 to perform read quality-trimming, alignment, and DNA-level genotype calculation   

## Inputs
There are *two* required input files:
- Fastq files (can be gzipped or not)
- A csv file (examples provided in `example_csv`), explanation of the columns is below

## Outputs:
- Protein-level genotype frequencies table in a csv file
- CRISPResso2 output that includes
  - Read aligning rate
  - Sequence-level genotype frequencies table
  - read-to-genotype assignments information

## Installation

**NOTE**: *if you installed DeepGenotype before 2025-01-15, please reinstall DeepGenotype to update CRISPResso2 to 2.3.1 to enable read quality-trimming.*

create a conda environment and activate it
```shell
module load anaconda # if on the hpc
conda create -n DeepGenotype python=3.9
conda activate DeepGenotype
```
install CRISPResso2
```shell
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install CRISPResso2==2.3.1
```
verify CRISPResso2 installation
```shell
CRISPResso -h
```
clone DeepGenotype repo and install dependencies
```shell
git clone https://github.com/czbiohub-sf/DeepGenotype
cd DeepGenotype
pip install . # or pip install biopython==1.78 pandas requests openpyxl==3.1.2

```
verify DeepGenotype installation
```shell
cd DeepGenotype # must be in the DeepGenotype/DeepGenotype directory
python DeepGenotype.py
```

&nbsp;
## Usage:
```shell
cd DeepGenotype # must be in the DeepGenotype/DeepGenotype directory
python DeepGenotype.py --path2csv example_csv/test.csv --path2workDir test_dir/ --path2fastqDir test_dir/fastq_dir/
```
All paths are relative to `DeepGenotype.py`  
Please make sure the following two python scripts are in the same directory as DeepGenotype.py:  
 &nbsp;&nbsp;&nbsp; process_alleles_freq_table_INS.py  
 &nbsp;&nbsp;&nbsp; process_alleles_freq_table_SNP.py  


#### Optional arguments
--fastq_R1_suffix &nbsp;&nbsp; (default "_R1_001.fastq.gz")  
--fastq_R2_suffix &nbsp;&nbsp; (default "_R2_001.fastq.gz")  
--single_fastq_suffix &nbsp;&nbsp; (use this option for **single-ended** reads as well as **pacbio** reads, need to specific the suffix, e.g.: fastq.gz)  
--quantification_window_size &nbsp;&nbsp; (default 50, which overrides CRISPResso2's default of 1)   


&nbsp;
## To run pacbio test dataset (insertion mode)
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
## Example 1: To run MiSeq test dataset (insertion mode)
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
## Example 2: to run MiSeq test dataset (SNP mode)
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
