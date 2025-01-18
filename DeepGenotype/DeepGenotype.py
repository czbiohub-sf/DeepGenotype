import argparse
import sys
import linecache
import pandas as pd
import os
from subprocess import Popen
import shutil
import warnings
import logging
import re
import unicodedata
warnings.filterwarnings('ignore')

#################
#custom logging #
#################
BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

#The background is set with 40 plus the number of the color, and the foreground with 30

#These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

def formatter_message(message, use_color = True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message

COLORS = {
    'WARNING': YELLOW,
    'INFO': WHITE,
    'DEBUG': BLUE,
    'CRITICAL': YELLOW,
    'ERROR': RED
}

class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, use_color = True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)

# Custom logger class with multiple destinations
class ColoredLogger(logging.Logger):
    FORMAT = "[$BOLD%(name)-1s$RESET][%(levelname)-1s]  %(message)s " #($BOLD%(filename)s$RESET:%(lineno)d)
    COLOR_FORMAT = formatter_message(FORMAT, True)
    def __init__(self, name):
        logging.Logger.__init__(self, name, logging.DEBUG)

        color_formatter = ColoredFormatter(self.COLOR_FORMAT)

        console = logging.StreamHandler()
        console.setFormatter(color_formatter)

        self.addHandler(console)
        return

logging.setLoggerClass(ColoredLogger)
#logging.basicConfig()
log = logging.getLogger("DeepGenotype.py")
log.setLevel(logging.INFO) #set the level of warning displayed

############
#Arguments #
############
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser = MyParser(description='This script runs CRISPResso and calculate customized genotypes')
    parser.add_argument('--path2csv', default="", type=str,
                        help='path to a csv file containing sample information\n *required columns*: Sample_ID, gene_name, ENST_id, gRNA_sequence, WT_amplicon_sequence, HDR_amplicon_sequence ', metavar='')
    parser.add_argument('--path2workDir', default="", type=str, help='path to a directory for saving output from CRISPResso', metavar='')
    parser.add_argument('--path2fastqDir', default="", type=str, help='path to a directory containing the fastq files', metavar='')
    parser.add_argument('--quantification_win_size', default=50, type=int, help='quantification window size, default = 50', metavar='')
    parser.add_argument('--fastq_R1_suffix', default="_R1_001.fastq.gz", type=str, help='(optional) suffix to add to sample ID to map to fastq files, e.g. R1_001.fastq.gz', metavar='')
    parser.add_argument('--fastq_R2_suffix', default="_R2_001.fastq.gz", type=str, help='(optional) suffix to add to sample ID to map to fastq files, e.g. R2_001.fastq.gz', metavar='')
    parser.add_argument('--single_fastq_suffix', default="", type=str, help='(optional) suffix to add to sample ID to map to fastq files, e.g. R2_001.fastq.gz', metavar='')
    parser.add_argument('--fastp_options_string', default="--cut_front --cut_tail --cut_mean_quality 30 --cut_window_size 30", type=str, help='options to pass to fastp, the default is to do quality trimming from both ends of each read, using a slide window of 4 and a mean quality threshold of 20, see fastp documentation for more options', metavar='')
    parser.add_argument('--n_processes', default=1, type=int, help='number of cores to use for parallel processing, use with caution since increasing this parameter will significantly increase the memory required ', metavar='')
    parser.add_argument('--skip_crispresso', action='store_true', default=False, help='skip CRISPResso if results already exist')
    config = parser.parse_args()
    if len(sys.argv) == 1:  # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

config = vars(parse_args())
path2workDir = os.path.join(os.getcwd(), config["path2workDir"])
path2fastqDir = os.path.join(os.getcwd(), config['path2fastqDir'])
path2csv= os.path.join(os.getcwd(), config['path2csv'])
quantification_win_size= config['quantification_win_size']
fastq_R1_suffix= config['fastq_R1_suffix']
fastq_R2_suffix= config['fastq_R2_suffix']
single_fastq_suffix= config['single_fastq_suffix']
fastp_options_string= config['fastp_options_string']
path2_stdout = os.path.join(path2workDir, "CRISPResso_run_logs")
path2_CRISPResso_out = os.path.join(path2workDir, "CRISPResso_outputs")
path2_allelsFreqTabs = os.path.join(path2workDir, "Alleles_freq_tables_with_genotypes")

if os.path.isdir(path2_stdout):
    shutil.rmtree(path2_stdout)
os.makedirs(path2_stdout)
if os.path.isdir(path2_CRISPResso_out) and not config['skip_crispresso']:
    shutil.rmtree(path2_CRISPResso_out)
os.makedirs(path2_CRISPResso_out, exist_ok=True)
if os.path.isdir(path2_allelsFreqTabs):
    shutil.rmtree(path2_allelsFreqTabs)
os.makedirs(path2_allelsFreqTabs)

wd = os.getcwd()  # save current working dir
os.chdir(path2_CRISPResso_out)  # change to the the folder containng the file to be zipped


###################
#text manipulation#
###################

def slugify(value): #adapted from the Django project

    value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore')
    value = re.sub(rb'[^\w\s-]', b'_', value).strip()
    value = re.sub(rb'[-\s]+', b'-', value)

    return value.decode('utf-8')

#####################
##      main       ##
#####################
def main():
    try:
        # check input paths
        if not os.path.isfile(path2csv):
            log.error(f"Can not locate the input csv file at {path2csv}, please check the path")
            log.info(f"Please fix the path to the input csv file and try again")
            sys.exit()
        if not os.path.isdir(os.path.join(path2workDir)):
            log.error(f"Can not locate the working directory at {path2workDir}, please check the path")
            log.info(f"Please fix the path to the working directory  and try again")
            sys.exit()

        # read input csv file
        df = pd.read_csv(os.path.join(os.getcwd(),path2csv))

        # check edit type
        if len(set(df['edit_type'])) >= 2:
            sys.exit("There are multiple \"edit_type\" values, please include only one type of edits (either INS or SNP, not both).\nThe reason behind this is: SNP and INS have different format of the output allele frequency spreadsheet generated from all samples")
        edit_type = list(set(df['edit_type']))[0]
        log.info(f"Genome edit type: {edit_type}")
        if not edit_type in ['INS', 'SNP']:
            log.error(f"{edit_type} is not support, only INS and SNP (case sensitive) are supported edit types")
            log.info(f"Please fix the inputs and try again")
            sys.exit()            

        # map edit_type to python scripts that process alleles freq tables
        script_path=''
        if edit_type == "INS":
            script_path = os.path.join(wd, "process_alleles_freq_table_INS.py")
            if not os.path.isfile(script_path):  # check script file existence
                log.error(f"Python script not found, please place file \"process_alleles_freq_table_INS.py\" in the same directory as \"DeepGenotype.py\"")
                log.error(f"...{row['Sample_ID']} was not processed")
                sys.exit()
        elif edit_type == "SNP":
            script_path = os.path.join(wd, "process_alleles_freq_table_SNP.py")
            if not os.path.isfile(script_path):  # check script file existence
                log.error(f"Python script not found, please place file \"process_alleles_freq_table_INS.py\" in the same directory as \"DeepGenotype.py\"")
                log.error(f"...{row['Sample_ID']} was not processed")
                sys.exit()
        log.debug(f"script path: {script_path}")

        #check csv columns
        keys2check = set(['Sample_ID','gene_name','WT_amplicon_sequence','HDR_amplicon_sequence', 'gRNA_sequence','edit_type'])
        if not keys2check.issubset(df.columns):
            log.error(f"Missing columns in the input csv file\n Required columns:\"Sample_ID\", \"gene_name\", \"WT_amplicon_sequence\", \"HDR_amplicon_sequence\", \"gRNA_sequence\", \"edit_type\"")
            log.info(f"Please fix the input csv file and try again")
            sys.exit()

        if not "payload_block_index" in df.columns:
            log.warning(f"The \"payload_block_index\" column in not found in the input csv file\nDefaulting \"payload_block_index\" to 1.\nPlease verify that there is only one block of SNPs or insertion in the HDR_amplicon relative to the wt_amplicon.\nThe first block is assumed to be the payload block.")


        #start processing samples through CRISPResso and recalculate allele frequency
        out_basename = os.path.basename(path2csv).strip(r".csv")
        with open(os.path.join(path2workDir,f"{out_basename}_genotype_freq.csv"), "w", buffering=1) as writehandle:
            if edit_type == "INS":
                writehandle.write("Edit_status,unedited,edited,edited,edited,edited,edited,edited,edited,,,\n")
                writehandle.write("Perfect_edit_status,-,perfect edit,imperfect edit,imperfect edit,imperfect edit,imperfect edit,imperfect edit,imperfect edit,,,\n")
                writehandle.write(f"Protein-level_genotype,wt_allele,HDR_perfect,wtProt_noPL,wtProt_okPL,mutProt_noPL,mutProt_okPL,mutProt_mutPL,wtProt_mutPL,num_post_filter_reads,weighted_avg_percent_identity,weighted_avg_num_of_mismatches\n") #write header
                writehandle.write("Sample_ID,-,-,-,-,-,-,-,-,-,-,-\n")
            elif edit_type == "SNP":
                #THIS IS THE OLD OUTPUT HEADER: writehandle.write("Sample_ID,wt_allele,HDR_perfect,wtProt_wtSNP,wtProt_hdrSNP,mutProt_wtSNP,mutProt_hdrSNP,mutProt_mutSNP,wtProt_mutSNP, (Prot=protein; SNP=SNP-of-interest; mutProt=mutation-in-protein-exclusing-SNP-site; hdrSNP=intended-protein-sequence-change-by-SNP; mutSNP=unintended-protein-sequence-change-by-SNP; wtSNP=unchanged-DNA-sequence-at-SNP-site)\n")  # write header
                writehandle.write("Edit_status,unedited,edited,edited,edited,edited,edited,edited,edited,,,\n")
                writehandle.write("Perfect_edit_status,-,perfect edit,imperfect edit,imperfect edit,imperfect edit,imperfect edit,imperfect edit,imperfect edit,,,\n")
                writehandle.write("Amino_acid_change_status,-,a.a change(s),a.a change(s),a.a change(s),a.a change(s),a.a change(s),a.a change(s),no a.a. change,,,\n")
                writehandle.write("Intended_edit_status,-,intended a.a. change(s) only,intended a.a. change(s) + synomous change(s),intended a.a change + unintended a.a. change(s),unintended a.a. change(s),unintended a.a. change(s),unintended a.a. change(s),synomous changes only and no intended a.a change,,,\n")
                writehandle.write("Protein-level_genotype,wt_allele,HDR_perfect,wtProt_hdrSNP,mutProt_hdrSNP,wtProt_mutSNP,mutProt_wtSNP,mutProt_mutSNP,wtProt_wtSNP,num_post_filter_reads,weighted_avg_percent_identity,weighted_avg_num_of_mismatches\n")
                writehandle.write("Sample_ID,-,-,-,-,-,-,-,-,-,-,-\n")
            # run CRISPResso for each sample_ID
            for index, row in df.iterrows():
                fq_ex_suffix = ""
                if 'Fastq_extra_suffix' in df.columns: #get extra suffix in the fastq file name
                    fq_ex_suffix = row["Fastq_extra_suffix"]
                fastq_r1 = f"{path2fastqDir}/{row['Sample_ID']}{fq_ex_suffix}{fastq_R1_suffix}"
                fastq_r2 = f"{path2fastqDir}/{row['Sample_ID']}{fq_ex_suffix}{fastq_R2_suffix}"
                fastq = f"{path2fastqDir}/{row['Sample_ID']}{fq_ex_suffix}{single_fastq_suffix}"

                #check fastq file
                if single_fastq_suffix=="":
                    if not os.path.isfile(fastq_r1):
                        log.error(f"...Can not locate fastq file: {fastq_r1}")
                        log.error(f"...{row['Sample_ID']} was not processed")
                        continue
                    if not os.path.isfile(fastq_r2):
                        log.error(f"...Can not locate fastq file: {fastq_r2}")
                        log.error(f"...{row['Sample_ID']} was not processed")
                        continue
                else:
                    if not os.path.isfile(fastq):
                        log.error(f"...Can not locate fastq file: {fastq}")
                        log.error(f"...{row['Sample_ID']} was not processed")
                        continue

                #check payload_block_index
                if 'payload_block_index' in row:
                    if not str(row['payload_block_index']).isdigit():
                        log.error(f"...the parameter \"payload_block_index\" needs to take an integer number")
                        log.error(f"...{row['Sample_ID']} was not processed")
                        continue

                #create fastq list for building the CRISPResso command
                if single_fastq_suffix == "":
                    fastq_list = [f"--fastq_r1", f"{fastq_r1}",
                                 f"--fastq_r2", f"{fastq_r2}"]
                else:
                    fastq_list = [f"--fastq_r1", f"{fastq}"]
                log.debug(f"{fastq_list}")
                num_fastq_files = int(len(fastq_list)/2)


                #get sample name and predict CRISPResso output directory
                sample_name = row['Sample_ID']
                #sample_name = slugify(sample_name)
                current_CRISPResso_out_dir = os.path.join(path2_CRISPResso_out,f"CRISPResso_on_{sample_name}")

                #check if CRISPResso results already exist
                skip_crispresso_flag = config['skip_crispresso']
                crispresso_results_exist = os.path.isfile(os.path.join(current_CRISPResso_out_dir,"Alleles_frequency_table.zip"))
                if skip_crispresso_flag and crispresso_results_exist:
                    log.info(f"...CRISPResso results already exist for {row['Sample_ID']}, skipping CRISPResso")
                else:
                    # build CRISPResso command
                    command = [f"CRISPResso"] + fastq_list + [
                            f"--amplicon_seq", f"{row['WT_amplicon_sequence']}",
                            f"--expected_hdr_amplicon_seq", f"{row['HDR_amplicon_sequence']}",
                            f"--amplicon_name", f"{row['gene_name']}",
                            f"--guide_seq", f"{row['gRNA_sequence']}",
                            f"--name", f"{row['Sample_ID']}",
                            f"--quantification_window_size", f"{quantification_win_size}",
                            f"--n_processes", f"{config['n_processes']}",
                            f"--fastp_options_string", f"{config['fastp_options_string']}",
                            #f"--plot_window_size", "20", # reduce plot size
                            #f"--max_rows_alleles_around_cut_to_plot", "20" # reduce plot size
                            ]

                    log.info(f"Processing sample: {row['Sample_ID']}")
                    if num_fastq_files == 2:
                        log.info(f"...running CRISPResso on {num_fastq_files} paired fastq files")
                    else:
                        log.info(f"...running CRISPResso on {num_fastq_files} fastq file")

                    #check amplicon length
                    if edit_type == "SNP": # For SNP, check if length of wt amp and HDR amp are the same
                        if len(row['WT_amplicon_sequence']) != len(row['HDR_amplicon_sequence']):
                            log.error(f"...for SNP analysis, the length must match between the wt amplicon and the HDR amplicon")
                            log.error(f"...{row['Sample_ID']} was not processed")
                            continue

                    #run CRISPResso
                    path_to_stderr_file = os.path.join(path2_stdout, f"{row['Sample_ID']}.stderr.txt")
                    path_to_stdout_file = os.path.join(path2_stdout, f"{row['Sample_ID']}.stdout.txt")
                    mystdput = open(path_to_stdout_file, 'w+')
                    mystderr = open(path_to_stderr_file, 'w+')
                    p = Popen(command, stdout=mystdput, stderr=mystderr, universal_newlines=True)

                    p.communicate()  # wait for the commands to process

                #process CRISPresso result: allele frequency table

                #CRISPResso_out_dir
                #CRISPResso changes "." to "_" in the output folder, so we need to account for this behavior here
                #NOTE: upgrading CRISPResso from v2.2.6 to 2.2.14 no longer have the above behavior, thus commenting the following line out

                #build and execute the shell command
                if os.path.isfile(os.path.join(current_CRISPResso_out_dir,"Alleles_frequency_table.zip")):
                    command2 = [f"{sys.executable}",f"{script_path}",
                                f"--path", f"{current_CRISPResso_out_dir}",
                                f"--allele_freq_file", f"Alleles_frequency_table.zip",
                                f"--wt_amp", f"{row['WT_amplicon_sequence']}",
                                f"--HDR_amp", f"{row['HDR_amplicon_sequence']}",
                                f"--ENST_ID", f"{row['ENST_id']}"
                                ]
                    
                    if 'payload_block_index' in row:
                        command2 = command2 + [f"--payload_block_index",
                                             f"{row['payload_block_index']}"]
                        log.debug(f"payload_block_index: {row['payload_block_index']}")
                    else:
                        command2 = command2 + [f"--payload_block_index",
                                             f"1"]

                    #run command2
                    p = Popen(command2, universal_newlines=True)
                    log.info(f"...parsing allele frequency table and re-calculating allele frequencies")
                    p.communicate()  # wait for the commands to process
                else:
                    log.error(f"...cannot find CRISPResso output file: Alleles_frequency_table.zip with the following path:")
                    path_not_found = os.path.join(current_CRISPResso_out_dir,"Alleles_frequency_table.zip")
                    log.error(f"{path_not_found}")
                    log.error(f"...{row['Sample_ID']} was not processed")
                if os.path.isfile(os.path.join(current_CRISPResso_out_dir, "genotype_frequency.csv")):
                    with open(os.path.join(current_CRISPResso_out_dir, "genotype_frequency.csv"), "r") as handle:
                        next(handle)
                        writehandle.write(f"{row['Sample_ID']},")
                        writehandle.write(handle.readline())
                    log.info(f"...done")
                else:
                    current_result_file = os.path.join(current_CRISPResso_out_dir, "genotype_frequency.csv")
                    log.error(f"cannot find intermediate file {current_result_file}")
                    log.error(f"...{row['Sample_ID']} was not processed")

                # move the genotype zip file
                zip_path = os.path.join(current_CRISPResso_out_dir, "Alleles_frequency_table.zip")
                old_zip_path = "_".join([zip_path.rstrip(r'.zip'), "genotype.zip"])
                moved_zip_path = os.path.join(path2_allelsFreqTabs,f"{row['Sample_ID']}_alFreqRecal.zip")
                shutil.move(old_zip_path, moved_zip_path)

            # write genotype explanation
            writehandle.write(f"\n------------------------below are the abbreviations and genotype explanations---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
            if edit_type == "INS":
                    writehandle.write(f"\nAbbreviations:\n"
                                    f"    wt = wildtype\n"
                                    f"    PL = payload\n"
                                    f"    Prot = protein sequence without the payload (hereafter referred to as \"protein sequence\")\n"
                                    f"    mut = mutated: not matching neither the wt nor HDR amplicon at the protein level\n"
                                    f"    HDR = protein sequence matching the HDR amplicon\n"
                                    f"    okPL = payload whose sequence is correct at the protein level\n"
                                    f"Genotype explanation:\n"
                                    f"    wt_allele: wildtype allele - the *whole* amplicon sequence perfectly matched the wt amplicon at the *nucleotide* level\n"
                                    f"    HDR_perfect: the *whole* amplicon sequence perfectly matched the HDR amplicon at the *nucleotide* level. Please note that the 'HDR perfect' in CRISPResso results only looks at the payload sequence.\n"
                                    f"    wtProt_noPL: protein sequence matched wt protein sequence. no payload detected\n"
                                    f"    wtProt_okPL: protein sequence matched wt protein sequence. payload is correct at the protein level\n"
                                    f"    mutProt_noPL: protein sequence is mutated (see mut definition above). no payload detected\n"
                                    f"    mutProt_okPL: protein sequence is mutated (see mut definition above). payload is correct at the protein level\n"
                                    f"    mutProt_mutPL: protein sequence is mutated. payload is mutated (see mut definition above)\n"
                                    f"    wtProt_mutPL: protein sequence matched wt protein sequence. payload is mutated (see mut definition above)\n"
                                    f"Notes:\n"
                                    f"    num_post_filter_reads: number of reads that passed the quality filter\n"
                                    f"    weighted_avg_percent_identity: weighted average of the percent identity of the reads (that aligned to the HDR amplicon)\n"
                                    f"    weighted_avg_num_of_mismatches: weighted average of the number of mismatches of the reads (that aligned to the HDR amplicon)\n"
                                        )
            elif edit_type == "SNP":
                    writehandle.write(f"\nAbbreviations:\n"
                                    f"    wt = wildtype\n"
                                    f"    SNP = SNP-payload\n"
                                    f"    Prot = protein sequence without the payload (hereafter referred to as \"protein sequence\")\n"
                                    f"    mut = mutated: not matching neither the wt nor HDR amplicon at the protein level\n"
                                    f"    HDR = protein sequence matching the HDR amplicon\n"
                                    f"Genotype explanation:\n"
                                    f"    wt_allele: wildtype allele - the *whole* amplicon sequence perfectly matched the wt amplicon at the *nucleotide* level\n"
                                    f"    HDR_perfect: the *whole* amplicon sequence perfectly matched the HDR amplicon at the *nucleotide* level. Please note that the 'HDR perfect' in CRISPResso results only looks at the payload sequence.\n"
                                    f"    wtProt_hdrSNP: protein sequence matched wt protein sequence. SNP-payload-associated amino acid matched the intended change\n"
                                    f"    mutProt_hdrSNP: protein sequence is mutated (see mut definition above). SNP-payload-associated amino acid matched the intended change\n"
                                    f"    wtProt_mutSNP: protein sequence matched wt protein sequence. SNP-payload-associated amino acid is mutated\n"
                                    f"    mutProt_wtSNP: protein sequence is mutated (see mut definition above). SNP-payload-associated amino acid matched wt\n"
                                    f"    mutProt_mutSNP: protein sequence is mutated. SNP-payload-associated amino acid is mutated (see mut definition above)\n"
                                    f"    wtProt_wtSNP: protein sequence matched wt protein sequence. SNP-payload-associated amino acid matched wt\n"
                                    f"Notes:\n"
                                    f"    num_post_filter_reads: number of reads that passed the quality filter\n"
                                    f"    weighted_avg_percent_identity: weighted average of the percent identity of the reads (that aligned to the HDR amplicon)\n"
                                    f"    weighted_avg_num_of_mismatches: weighted average of the number of mismatches of the reads (that aligned to the HDR amplicon)\n"
                                        )

        # convert csv to xlsx
        csv_path = os.path.join(path2workDir,f"{out_basename}_genotype_freq.csv")
        command3 = [f"python", f"{os.path.join(wd,'csv2xlsx.py')}","--path2csv", csv_path , "--mode", f"{edit_type}"]
        p = Popen(command3, universal_newlines=True)
        p.communicate()  # wait for the commands to process
        
        # all done
        log.info(f"Done processing all samples in the input csv file: {path2csv}")
        log.info(f"CRISPResso results are in folder: {path2_CRISPResso_out}")
        excel_path = csv_path.replace(".csv", ".xlsx")
        log.info(f"DeepGenotype output files are {csv_path} and {os.path.basename(excel_path)}")



        os.chdir(wd)  # change to the saved working dir

    except Exception as e:
        print("Unexpected error:", str(sys.exc_info()))
        print("additional information:", e)
        PrintException()

##########################
## function definitions ##
##########################
def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

if __name__ == "__main__": main()
