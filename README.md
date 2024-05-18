# DeepGenotype
Calculates the frequencies of protein-level mutations from deep-sequencing reads of CRISPR-edited cells

## Features
- Calculates genotypes with respect to protein/payload expressiblity and correctness  
- Automatically processes a list of samples 
- Supports CRISPR editing types: tagging/insertion and SNP/base-editing  
- Works with both Illumina and PacBio reads  
- Automatically finds coding regions  
- Invokes CRISPResso2 to perform read processing and alignment   

## Inputs
There are *two* required input files:
- Fastq files (can be gzipped or not)
- A csv file (examples provided in `example_csv`), see below.

### csv file columns
The input csv should contain the following columns with the exact names
  - Sample_ID (e.g. mNGplate19_sorted_A2_DDX6-C)  
    ***Important note***: For paired-end sequencing, only one Sample_ID is needed. We automatically find both R1 and R2 fastq files.   
     Check fastq file suffix parameters `--fastq_R1_suffix` and `--fastq_R1_suffix` in the `Usage` section.  
    Also check if you need `Fastq_extra_suffix` (below)

  - gene_name (e.g. DDX6)  
  - ENST_id (e.g. ENST00000620157)  
  - WT_amplicon_sequence
  - HDR_amplicon_sequence
  - gRNA_sequence
  - edit_type (e.g. INS or SNP, note that deletions, DEL is not supported at this point)  
      &nbsp;&nbsp;&nbsp; INS = insertion, SNP = single nucleotide polymorphism, DEL = deletion  
  - SNP_payload_cluster (e.g. 1 or 2 ...)
      &nbsp;&nbsp;&nbsp; Only needed when edit_type = SNP   
      &nbsp;&nbsp;&nbsp; This defines which cluster of SNPs is the **payload** (TODO: clarify how the non-payload SNPs will affect the perfect HDR, and wt-protein rates)  
      &nbsp;&nbsp;&nbsp; Clusters are ordered from left to right in respect to the amplicon sequence.  
      &nbsp;&nbsp;&nbsp; For examplp: there are 2 clusters of SNPs, the first cluster is a recut SNP and the second cluster of SNPs is of interest (payload), SNP_cluster should be set to 2, and the first cluster of SNP will only be anlayzed for mutations if its in the coding region. 
      
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
## Outputs:
- Protein-level genotype frequencies table in a csv file
- CRISPResso2 output that includes
  - Read aligning rate
  - Sequence-level genotype frequencies table
  - read-to-genotype assignments information

&nbsp;
## Installation

create a conda environment and activate it
```
conda create -n DeepGenotype python=3.9
conda activate DeepGenotype
```
install CRISPResso2
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install CRISPResso2==2.2.14
```
verify CRISPResso2 installation
```
CRISPResso -h
```
install Python packages and DeepGenotype
```
pip install biopython==1.78 pandas requests openpyxl==3.1.2
git clone https://github.com/czbiohub-sf/DeepGenotype
```
verify DeepGenotype installation
```
cd DeepGenotype
python DeepGenotype.py
```

&nbsp;
## Usage:
```
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
```
module load anaconda
conda activate DeepGenotype
```
run DeepGenotype
```
python DeepGenotype.py \
--path2csv example_csv/test_pacbio.csv \
--path2workDir test_PacBio \
--path2fastqDir test_PacBio/fastq \
--single_fastq_suffix .fastq
```
***NOTE***: to run DeepGenotype in the background (and thus safe to close the terminal), preprend `nohup` and append `&` to the command:
```
nohup python DeepGenotype.py \
--path2csv example_csv/test_pacbio.csv \
--path2workDir test_PacBio \
--path2fastqDir test_PacBio/fastq \
--single_fastq_suffix .fastq &
```
To check the terminal output (while running in the background
```
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
## To run MiSeq test dataset (insertion mode)
load conda, and activate the DeepGenotype conda environment
```
module load anaconda
conda activate DeepGenotype
```
run DeepGenotype in the background (and thus safe to close the terminal)
```
nohup python DeepGenotype.py \
--path2csv example_csv/test_INS.csv \
--path2workDir test_MiSeq_INS \
--path2fastqDir test_MiSeq_INS/fastq &
```
To check the terminal output (while running in the background
```
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
## To run MiSeq test dataset (SNP mode)
load conda, and activate the DeepGenotype conda environment
```
module load anaconda
conda activate DeepGenotype
```
run DeepGenotype in the background (and thus safe to close the terminal)
```
nohup python DeepGenotype.py \
--path2csv example_csv/test_MiSeq_SNP.csv \
--path2workDir test_MiSeq_SNP \
--path2fastqDir test_MiSeq_SNP/fastq &
```
To check the terminal output (while running in the background
```
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
## Helper commands

### Paginated view of fastq files
compressed
```
gzip -c example.fastq.gz | less
```
uncompressed
```
less example.fastq
```

### Count the number of lines in a fastq file (divide by 4 you'll get the read count)
compressed
```
gzip -c example.fastq.gz | wc -l
```
uncompressed
```
cat example.fastq | wc -l
```

### Count the number of line with a specific sequence in a fastq file
compressed
```
gzip -c example.fastq.gz | grep -c 'your-sequence-here'
```
uncompressed
```
cat example.fastq | grep -c 'your-sequence-here'
```
for matching sequences at the beginning of the line, add a `^`:  `'^your-sequence-here'`  
for matching sequences at the end of the line, add a `$`:  `'your-sequence-here$'`



# License
This project is licensed under the BSD 3-Clause license - see the LICENSE file for details.
