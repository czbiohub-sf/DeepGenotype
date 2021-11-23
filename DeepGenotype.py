import argparse
import sys
import linecache
import pandas as pd
import os
import subprocess
from subprocess import Popen

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
    parser.add_argument('--fastq_R1_suffix', default="R1_001.fastq.gz", type=str, help='(optional) suffix to add to sample ID to map to fastq files, e.g. R1_001.fastq.gz', metavar='')
    parser.add_argument('--fastq_R2_suffix', default="R2_001.fastq.gz", type=str, help='(optional) suffix to add to sample ID to map to fastq files, e.g. R2_001.fastq.gz', metavar='')
    parser.add_argument('--sample_name_addon', default="R2_001.fastq.gz", type=str, help='(optional) name of a column in csv file containing additional suffix to help map sample name to fastq file name', metavar='')
    config = parser.parse_args()
    if len(sys.argv) == 1:  # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

config = vars(parse_args())
path2workDir = config["path2workDir"]
path2fastqDir = config['path2fastqDir']
path2csv= config['path2csv']
quantification_win_size= config['quantification_win_size']
fastq_R1_suffix= config['fastq_R1_suffix']
fastq_R2_suffix= config['fastq_R2_suffix']
sample_name_addon = config["sample_name_addon"]

path2_stdout = os.path.join(path2workDir, "CRISPResso_run_logs")
path2_CRISPResso_out = os.path.join(path2workDir, "CRISPResso_outputs")

#####################
##      main       ##
#####################
def main():
    try:
        # read input csv file
        df = pd.read_csv(os.path.join(path2csv))

        # run CRISPResso for each sample_ID
        for index, row in df.iterrows():
            fastq_r1 = f"{path2fastqDir}/{row['Sample_ID']}{row[sample_name_addon]}{fastq_R1_suffix}"
            fastq_r2 = f"{path2fastqDir}/{row['Sample_ID']}{row[sample_name_addon]}{fastq_R2_suffix}"
            command = [f"CRISPResso",
                       f"--fastq_r1", f"{fastq_r1}",
                       f"--fastq_r2", f"{fastq_r2}",
                       f"--amplicon_seq", f"{row['WT_amplicon_sequence']}",
                       f"--expected_hdr_amplicon_seq", f"{row['HDR_amplicon_sequence']}",
                       f"--amplicon_name", f"{row['gene_name']}",
                       f"--guide_seq", f"{row['gRNA_sequence']}",
                       f"--name", f"{row['Sample_ID']}.out",
                       f"--quantification_window_size", f"{quantification_win_size}",
                       f"--output_folder", f"{path2_CRISPResso_out}",
                       ]

            # print(command)
            path_to_stderr_file = os.path.join(path2_stdout, f"{row['Sample_ID']}.stderr.txt")
            path_to_stdout_file = os.path.join(path2_stdout, f"{row['Sample_ID']}.stdout.txt")
            mystdput = open(path_to_stdout_file, 'w+')
            mystderr = open(path_to_stderr_file, 'w+')
            p = Popen(command, stdout=mystdput, stderr=mystderr, universal_newlines=True)
            print(f"Processing sample: {row['Sample_ID']} with CRISPResso ...")
            p.communicate()  # now wait plus that you can send commands to process
        print("done")

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
