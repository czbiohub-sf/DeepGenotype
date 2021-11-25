# DeepGenotype
Calculate genotype frequencies in CRISPR genome editing experiments using deep sequencing reads

## Inputs
- Fastq files
- A csv file containing the following columes:
  - Sample_ID (e.g. mNGplate19_sorted_A2_DDX6-C )
  - gene_name (e.g. DDX6)
  - WT_amplicon_sequence
  - HDR_amplicon_sequence
  - gRNA_sequence
  - (Optional but important) Sample_name_addon    
    
      This columns helps to map Sample_ID to corresponding fastq files in addition to the two standard fastq suffixes below :
      "R1_001.fastq.gz"
      "R2_001.fastq.gz"
      
      For example if your fastq file names are:  
      &nbsp;&nbsp;&nbsp; mNGplate19_sorted_A2_DDX6-C_S90_R1_001.fastq.gz  
      &nbsp;&nbsp;&nbsp; mNGplate19_sorted_A2_DDX6-C_S90_R2_001.fastq.gz  
      Then you should add "_S90_" to the "Sample_name_addon" column

## Usuage:
Python DeepGenotype.
--path2csv "/Users/duo.peng/My Drive/protospaceX_collaboration-20211111T220826Z-001/protospaceX_collaboration/example_mNG11_data/mNGplate19_MiSeq_analysis_sample_sheet.csv" 
--path2workDir "/Users/duo.peng/My Drive/genotyping_pipeline/testing2" 
--path2fastqDir "/Users/duo.peng/My Drive/protospaceX_collaboration-20211111T220826Z-001/protospaceX_collaboration/example_mNG11_data/plate19fastqs" 
--sample_name_addon _S_
## Dependencies
Bio   
CRISPRESSO2
  
