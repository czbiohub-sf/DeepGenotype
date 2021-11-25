# DeepGenotype
Calculate genotype frequencies in CRISPR genome editing experiments using deep sequencing reads

## Inputs
- Fastq files
- A csv file containing columns with the exact names:
  - Sample_ID (e.g. mNGplate19_sorted_A2_DDX6-C)
  - gene_name (e.g. DDX6)  
  - ENST_id (e.g. ENST00000620157)  
  - WT_amplicon_sequence
  - HDR_amplicon_sequence
  - gRNA_sequence
  - Fastq_extra_suffix (Optional)     
        
      This column contains extra suffixs that help to map Sample_ID to corresponding fastq files
      in addition to the two standard fastq suffixes below:  
      "_R1_001.fastq.gz"  
      "_R2_001.fastq.gz"  
      
      For example if your fastq file names are:  
      &nbsp;&nbsp;&nbsp; mNGplate19_sorted_A2_DDX6-C_S90_R1_001.fastq.gz  
      &nbsp;&nbsp;&nbsp; mNGplate19_sorted_A2_DDX6-C_S90_R2_001.fastq.gz  
      Then you should add "_S90_" to the "Sample_name_addon" column

## Usuage:
```
Python DeepGenotype.py --path2csv test_dir/test.csv --path2workDir test_dir/ --path2fastqDir test_dir/fastq_dir/ --sample_name_addon addon
```

## Dependencies
Python > 3.7  
Bio   
CRISPRESSO2
  
