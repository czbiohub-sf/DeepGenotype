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
- Fastq files
- A csv file (examples provided) containing columns with the exact names:
  - Sample_ID (e.g. mNGplate19_sorted_A2_DDX6-C)
  - gene_name (e.g. DDX6)  
  - ENST_id (e.g. ENST00000620157)  
  - WT_amplicon_sequence
  - HDR_amplicon_sequence
  - gRNA_sequence
  - edit_type (e.g. INS or SNP, note that deletions, DEL is not supported at this point)  
      &nbsp;&nbsp;&nbsp; INS = insertion, SNP = single nucleotide polymorphism, DEL = deletion  
  - SNP_payload_cluster (e.g. 1 or 2 ...) (TODO: check if this setting is really necessary)  
      &nbsp;&nbsp;&nbsp; Only needed when edit_type = SNP   
      &nbsp;&nbsp;&nbsp; This defines which cluster of SNPs is the **payload**  
      &nbsp;&nbsp;&nbsp; Clusters are ordered from left to right in respect to the amplicon sequence.  
      &nbsp;&nbsp;&nbsp; For examplp: there are 2 clusters of SNPs, the first cluster is a recut SNP and the second cluster of SNPs is of interest (payload), SNP_cluster should be set to 2, and the first cluster of SNP will only be anlayzed for mutations if its in the coding region. 
      
  - Fastq_extra_suffix (Optional) 
     
      &nbsp;&nbsp;&nbsp;Extra suffix needed for mapping Sample_ID to corresponding fastq files
 
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
- Genotype frequencies for each sample (.csv file)
- Alleles frequencies table (A folder containing a table of read-to-genotype assignments for each sample)

&nbsp;
## Usage:
```
python DeepGenotype.py --path2csv test_dir/test.csv --path2workDir test_dir/ --path2fastqDir test_dir/fastq_dir/
```
Please make sure the following two python scripts are in the same directory as DeepGenotype.py:  
 &nbsp;&nbsp;&nbsp; process_alleles_freq_table_INS.py  
 &nbsp;&nbsp;&nbsp; process_alleles_freq_table_SNP.py  

#### Optional aruments
--fastq_R1_suffix &nbsp;&nbsp; (default "_R1_001.fastq.gz")  
--fastq_R2_suffix &nbsp;&nbsp; (default "_R2_001.fastq.gz")  
--single_fastq_suffix &nbsp;&nbsp; (use this option for **single-ended** reads as well as **pacbio** reads, need to specific the suffix, e.g.: fastq.gz)  
--quantification_window_size &nbsp;&nbsp; (default 50, which overrides CRISPResso2's default of 1)   

&nbsp;
## Dependencies
- python >= 3.7

#### python packages  
- biopython  (v1.78)  
- pandas
- requests

To create an environment and install python + dependencies:
```
conda create -n DeepGenotype python=3.8
conda activate DeepGenotype
pip install biopython==1.78 pandas requests
```

#### python-based standalone software  
- CRISPResso2  
&nbsp;&nbsp;&nbsp;For installation, follow instructions here https://github.com/pinellolab/CRISPResso2  
&nbsp;&nbsp;&nbsp;It is recommended to create a conda environment to install CRISPREesso2 and other dependencies listed below 


&nbsp;
## Developed by:
Duo Peng&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (Development)  
Abigail Glascock &nbsp; &nbsp;                                         (Testing)  
Joan Wong &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;      (Supervision)  
