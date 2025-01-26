#load packages
import argparse, sys, os, re, zipfile, logging, warnings
from Bio.pairwise2 import align
from Bio.pairwise2 import format_alignment
from utils import *
warnings.filterwarnings('ignore')
log = logging.getLogger("Process_alleles_freq_table.py")
log.setLevel(logging.ERROR) #set the level of warning displayed

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='This script processes the allele frequency table output from CRISPResso2')
    parser.add_argument('--path', default="", type=str, help='path to the allele freq. table output zip file')
    parser.add_argument('--allele_freq_file', default="", type=str, help='the name of the freq. table output zip file')
    parser.add_argument('--wt_amp', default="", type=str, help='sequence of the wt amplicon')
    parser.add_argument('--HDR_amp', default="", type=str, help='sequence of the HDR amplicon')
    parser.add_argument('--ENST_ID', default="", type=str, help='ENST ID')
    parser.add_argument('--payload_block_index', default=1, type=int, help='if the SNP is downstream of re-cut-preventing SNP, set this to 2')
    parser.add_argument('--min_reads_for_genotype', default=3, type=int, help='minimum number of reads for genotype to be considered successful')

    config = parser.parse_args()
    if len(sys.argv)==1: # print help message if arguments are not valid
        parser.print_help()
        sys.exit(1)
    return config

config = vars(parse_args())

#####################
##      main       ##
#####################    
def main():
    try: 
        zip_dir = config['path']
        zip_name = config['allele_freq_file']
        zip_path = os.path.join(zip_dir,zip_name)
        wt_amp = config['wt_amp']
        HDR_amp = config['HDR_amp']
        ENST_ID = config['ENST_ID']
        SNP_block_of_interest = config['payload_block_index']
        SNP_idx = SNP_block_of_interest - 1

        ##################################################
        #prepration work before processing the alignments#
        ##################################################

        #fetch transcript
        try:
            mytranscript = fetch_ensembl_transcript(ensembl_transcript_id = ENST_ID, exon_annot = True)
        except Exception as e:
            log.error(f"Failed to fetch transcript {ENST_ID} from Ensembl. This may be due to connection issues, Ensembl server issues, or the ENST_ID is incorrect. the original error message is: {e}")
            sys.exit(1)

        #get cds seqs
        cdsSeqs = get_cds_seq_in_transcript(mytranscript)

        #concatenate cds segments
        whole_cds = ''
        if mytranscript.annotations['transcript_strand'] == 1:
            for part in cdsSeqs:
                whole_cds = whole_cds + part
        else:
            for part in reversed(cdsSeqs):
                whole_cds = whole_cds + part
        #print(f"whole CDS:\n{whole_cds}")

        #annotate coding window in wt_amplicon
        wt_amp_cds_coords = get_cds_coord_in_amplicon(whole_cds, wt_amp) # if strand = "-", the start and end cooresponds to the revcom sequence
        #print(f"\nwt amp cds coords:\n{wt_amp_cds_coords}")

        #get wt coding window in HDR amp (not including payloads)
        HDR_amp_cds_coords = get_cds_coord_in_amplicon(whole_cds, HDR_amp) # if strand = "-", the start and end cooresponds to the revcom sequence
        #print(f"\nHDR amp cds coords:\n{HDR_amp_cds_coords}")

        #get payload coordinates in HDR amp#
        #alignment to compute differences between wt and HDR amplicon
        aln = align.localms(wt_amp,HDR_amp,2, -0.1, -10, -5) #changed scoring to favor mismatch over gaps
        aln_wt = format_alignment(*aln[0]).split("\n")[0]
        aln_sym = format_alignment(*aln[0]).split("\n")[1]
        aln_HDR = format_alignment(*aln[0]).split("\n")[2]
        #check gaps in the alignment (there shouldn't be gaps)
        gap_pat = re.compile("-+")
        if re.search(gap_pat, aln_sym) is not None:
            sys.exit("Error: gaps found in the alignment between wt amplicon and HDR amplicon sequence")
        #find mismatches
        mm_pat = re.compile("\.+")
        #find gaps
        mms=[]
        for m_obj in mm_pat.finditer(aln_sym):
            mms.append(m_obj.span())
        SNP_coord_in_HDR_amp = mms # *regardless of strand*, the start and end cooresponds to the HDR amplicon
        #retrive the sequence at the mismatch sites
        SNP_wt = []
        SNP_HDR = []
        for item in mms:
            SNP_wt.append(aln_wt[item[0]:item[1]])
            SNP_HDR.append(aln_HDR[item[0]:item[1]])
        #get the SNP block of interest, and convert the SNP name to payload
        payload_coord_in_HDR_amp = SNP_coord_in_HDR_amp[SNP_idx]
        #get payload strand
        payload_strand = "+"
        if mytranscript.annotations['transcript_strand'] == -1:
            payload_strand = "-"
            

        #translate wt amplicon
        wt_translation=[]
        for coord_set in wt_amp_cds_coords:
            translation = translate_amp(seq=wt_amp, start = coord_set[0], end = coord_set[1], strand = coord_set[2], frame = coord_set[3])
            wt_translation.append(translation)
            #print(translation)

        #translate HDR amplicaon
        HDR_translation=[]
        for coord_set in HDR_amp_cds_coords:
            translation = translate_amp(seq=HDR_amp, start = coord_set[0], end = coord_set[1], strand = coord_set[2], frame = coord_set[3])
            HDR_translation.append(translation)
            #print(translation)

        payload_coord = payload_coord_in_HDR_amp # map legacy name to current name, #no need to convert coordinates, because they are already in reference to the + strand HDR amp
        HDR_amp_cds_coords_PosRef = get_amp_abs_coords(HDR_amp, HDR_amp_cds_coords) # convert coordinate so that they are all in reference to the + strand HDR amp, HDR_amp_cds_coords_PosRef canbe only used to check indels locations and *not translation*
        wt_amp_cds_coords_PosRef = get_amp_abs_coords(wt_amp, wt_amp_cds_coords) # convert coordinate so that they are all in reference to the + strand wt amp, wt_amp_cds_coords_PosRef canbe only used to check indels locations and *not translation*

        #translate payload
        payload_aa_pos = ''
        payload_HDR_translation = ""
        payload_wt_translation =  ""
        #check if payload is in cds
        payload_in_cds_flag = False
        for idx, coordset in enumerate(HDR_amp_cds_coords_PosRef):
            if payload_coord[0]>=coordset[0] and payload_coord[1]<=coordset[1]:
                payload_in_cds_flag = True
                #translate payload
                payload_aa_pos = find_payload_aa_pos(cds_st=coordset[0],cds_en=coordset[1],payload_st=payload_coord_in_HDR_amp[0],payload_en=payload_coord_in_HDR_amp[1],strand = coordset[2],frame = coordset[3])
                payload_HDR_translation = HDR_translation[idx][payload_aa_pos[0]: payload_aa_pos[1]]
                payload_wt_translation =  wt_translation[idx][payload_aa_pos[0]: payload_aa_pos[1]]        
        if (payload_in_cds_flag == False):
            sys.exit("ERROR: the SNP block of interest is not in the coding region")

                                

        #payload_translation = translate_payload(seq=HDR_amp, start = payload_coord_in_HDR_amp[0], end = payload_coord_in_HDR_amp[1], strand = payload_strand, frame = 1)

        #print(f"\nwt amp translation {wt_translation}")
        #print(f"HDR_translation {HDR_translation}")
        #print(f"payload_translation {payload_translation}")
 
        ########################
        #process the alignments#
        ########################
        with zipfile.ZipFile(zip_path, 'r') as myzip, open(os.path.join(zip_dir,'Alleles_frequency_table_genotype.txt'),'w') as writehandle:
            with myzip.open('Alleles_frequency_table.txt') as filehandle:
                next(filehandle) #skip header
                #print(f"Reference_Name\tRead_Status\tn_Reads\tperc_Reads\tn_deleted\tn_inserted\tn_mutated")
                writehandle.write(f"Aligned_Sequence	Reference_Sequence\tReference_Name\tRead_Status\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\tGenotype\t%Identity\t#Mismatches\n")
                for line in filehandle:
                    line_deco = line.decode()

                    fields = line_deco.rstrip().split("\t")
                    read = fields[0]
                    ref = fields[1]
                    Reference_Name = fields[2]
                    Read_Status = fields[3]
                    n_deleted = fields[4]
                    n_inserted = fields[5]
                    n_mutated = fields[6]
                    n_Reads = fields[7]
                    perc_Reads = fields[8]

                    if int(n_Reads) < config['min_reads_for_genotype']:
                        continue #skip if the number of reads is less than the minimum number of reads for genotype

                    writehandle.write(line_deco.rstrip()) # write the original line

                    ###########################################################################################################
                    #TRIM the gaps on both ends of the ref, these those gaps represent extra seq in the reads (not insertions)#
                    ###########################################################################################################
                    leading_gaps = []
                    trailing_gaps = []
                    if ref.startswith("-"):
                        leading_gaps = cal_leading_gaps(ref)
                        ref = ''.join(ref[idx] for idx in range(len(ref)) if not any([idx in range(st,en) for st,en in leading_gaps]))
                        read = ''.join(read[idx] for idx in range(len(read)) if not any([idx in range(st,en) for st,en in leading_gaps]))
                        #print(f"{read}\n{ref}", end='\n')
                    if ref.endswith("-"):
                        trailing_gaps = cal_trailing_gaps(ref)
                        ref = ''.join(ref[idx] for idx in range(len(ref)) if not any([idx in range(st,en) for st,en in trailing_gaps]))
                        read = ''.join(read[idx] for idx in range(len(read)) if not any([idx in range(st,en) for st,en in trailing_gaps]))
                        #print(f"{read}\n{ref}", end='\n')
                    # compute percentidentity and num_of_mismatches
                    num_of_mismatches = compute_num_of_mismatches(read, ref)
                    percent_identity = compute_percent_identity(ref, num_of_mismatches)
                    #############################
                    #Start calculating genotypes#
                    #############################            
                    #############################
                    #read mapped to HDR amplicon#
                    #############################
                    if Reference_Name == "HDR": # 
                        # compute percentidentity and nu
                        if read == ref:
                            #print("\tHDR allele (perfect HDR edit)", end="\n")
                            writehandle.write("\tperfect HDR edit\t")
                        ##########################
                        #mutations only, no indel#
                        ##########################
                        elif int(n_deleted)==0 and int(n_inserted)==0: 
                            #print(f"\n{read}\n{ref}\t{n_deleted}\t{n_inserted}\t{n_mutated}\t{n_Reads}\t{perc_Reads}", end='\n')
                            #check wt protein and payload
                            protein_correct = check_protein_and_payload(read, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["protein_correct"]
                            payload_genotype = check_protein_and_payload(read, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["payload_genotype"]

                            #produce output
                            if protein_correct==True:
                                #print(f"\twt protein + {payload_genotype} payload)")
                                writehandle.write(f"\twt protein + {payload_genotype} SNP\t")
                            elif protein_correct==False:
                                #print(f"wt protein + {payload_genotype} payload)")
                                writehandle.write(f"\tmutant protein + {payload_genotype} SNP\t")    

                        ########################
                        #indels in HDR allele  #
                        ########################
                        elif (int(n_deleted)!=0 and int(n_inserted)==0):  #deletion in the read，no insertion
                            #check if deletions are in wt cds
                            deletion_in_HDR_amp_cds_noPL_Flag = any([check_deletion_in_cds(read = read, amp_cds_coords = HDR_amp_cds_coords_PosRef)["deletion_in_cds_flag"],
                                                                    check_deletion_overlap_cds(read = read, amp_cds_coords = HDR_amp_cds_coords_PosRef)["deletion_overlap_cds_flag"]])        
                        
                            #check translation 
                            protein_correct = check_protein_and_payload(read, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["protein_correct"]
                            payload_genotype = check_protein_and_payload(read, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["payload_genotype"]
                            
                            if protein_correct == False:
                                deletion_in_HDR_amp_cds_noPL_Flag = True # hitch a ride
                            
                            #no need to check deletion in SNP, because the SNP translation will take care it this
                            #produce output
                            if deletion_in_HDR_amp_cds_noPL_Flag==True:
                                #print("\twt allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write(f"\tmutant protein + {payload_genotype} SNP\t")
                            elif deletion_in_HDR_amp_cds_noPL_Flag==False:
                                #print("\twt allele (wt protein + HDR payload)", end="\n")   
                                writehandle.write(f"\twt protein + {payload_genotype} SNP\t")       
                        elif (int(n_deleted)==0 and int(n_inserted)!=0):  #insertion in the read, no deletion
                            #check if insertion are in wt cds
                            insertion_in_HDR_amp_cds_noPL_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = HDR_amp_cds_coords_PosRef)["insertion_in_cds_flag"]
                            #check if insertion are in payload cds
                            insertion_in_payload_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = [payload_coord])["insertion_in_cds_flag"]
                            
                            #trim off inserted seq (in both reads and ref)
                            reg_gap_windows = cal_gap_win(seq = ref)
                            ref_trimmed = ''.join(ref[idx] for idx in range(len(ref)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))
                            read_trimmed = ''.join(read[idx] for idx in range(len(read)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))    
                            
                            #check translation 
                            protein_correct = check_protein_and_payload(read_trimmed, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["protein_correct"]
                            payload_genotype = check_protein_and_payload(read_trimmed, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["payload_genotype"]
                            
                            #summarize flag
                            mut_protein_summary_flag = any([insertion_in_HDR_amp_cds_noPL_Flag,
                                                        (not protein_correct)])        
                            
                            if insertion_in_payload_Flag:
                                payload_genotype = "mutant"
                            
                            #produce output
                            if mut_protein_summary_flag==True:
                                #print("\twt allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write(f"\tmutant protein + {payload_genotype} SNP\t")
                            elif mut_protein_summary_flag==False:
                                #print("\twt allele (wt protein + HDR payload)", end="\n")
                                writehandle.write(f"\twt protein + {payload_genotype} SNP\t") 
                                
                        elif(int(n_deleted)!=0 and int(n_inserted)!=0): #both insertion and deletion in the read
                            #check if insertion are in wt cds
                            insertion_in_HDR_amp_cds_noPL_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = HDR_amp_cds_coords_PosRef)["insertion_in_cds_flag"]
                            #check if insertion are in payload cds
                            insertion_in_payload_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = [payload_coord])["insertion_in_cds_flag"]            

                            #trim off inserted seq (in both reads and ref)
                            reg_gap_windows = cal_gap_win(seq = ref)
                            ref_trimmed = ''.join(ref[idx] for idx in range(len(ref)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))
                            read_trimmed = ''.join(read[idx] for idx in range(len(read)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))                     

                            #check if deletions are in wt cds
                            deletion_in_HDR_amp_cds_noPL_Flag = any([check_deletion_in_cds(read = read_trimmed, amp_cds_coords = HDR_amp_cds_coords_PosRef)["deletion_in_cds_flag"],
                                                                    check_deletion_overlap_cds(read = read_trimmed, amp_cds_coords = HDR_amp_cds_coords_PosRef)["deletion_overlap_cds_flag"]])
                            #check translation 
                            protein_correct = check_protein_and_payload(read_trimmed, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["protein_correct"]
                            payload_genotype = check_protein_and_payload(read_trimmed, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["payload_genotype"]
                            
                            #summarize flag
                            indel_in_HDR_amp_cds_noPL_Flag = any([insertion_in_HDR_amp_cds_noPL_Flag,
                                                                deletion_in_HDR_amp_cds_noPL_Flag, 
                                                                (not protein_correct)])
                    
                            #produce output
                            if (indel_in_HDR_amp_cds_noPL_Flag==True):
                                #print("\twt allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write(f"\tmutant protein + {payload_genotype} SNP\t")
                            elif (indel_in_HDR_amp_cds_noPL_Flag==False):
                                #print("\twt allele (wt protein + HDR payload)", end="\n")
                                writehandle.write(f"\twt protein + {payload_genotype} SNP\t")  
                    #############################
                    #read mapped to wt amplicon #
                    #############################        
                    else: 
                        if read == ref:
                            #print("\twt allele", end="\n")
                            writehandle.write("\twt allele\t")
                        ##########################
                        #mutations only, no indel#
                        ##########################
                        elif int(n_deleted)==0 and int(n_inserted)==0: 
                            #print(f"\n{read}\n{ref}\t{n_deleted}\t{n_inserted}\t{n_mutated}\t{n_Reads}\t{perc_Reads}", end='\n')
                            #check translation 
                            protein_correct = check_protein_and_payload(read, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["protein_correct"]
                            payload_genotype = check_protein_and_payload(read, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["payload_genotype"]

                            #produce output
                            if protein_correct==True:
                                #print(f"\twt protein + {payload_genotype} payload)")
                                writehandle.write(f"\twt protein + {payload_genotype} SNP\t")
                            elif protein_correct==False:
                                #print(f"\tmutant protein + {payload_genotype} payload)")
                                writehandle.write(f"\tmutant protein + {payload_genotype} SNP\t")
                        ########################
                        #indels in wt allele  #
                        ########################
                        elif (int(n_deleted)!=0 and int(n_inserted)==0):  #deletion in the read，no insertion
                            #check if deletions are in wt cds
                            deletion_in_wt_amp_cds_noPL_Flag = any([check_deletion_in_cds(read = read, amp_cds_coords = wt_amp_cds_coords_PosRef)["deletion_in_cds_flag"],
                                                                    check_deletion_overlap_cds(read = read, amp_cds_coords = wt_amp_cds_coords_PosRef)["deletion_overlap_cds_flag"]])
                            #check translation 
                            protein_correct = check_protein_and_payload(read, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["protein_correct"]
                            payload_genotype = check_protein_and_payload(read, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["payload_genotype"]
                            
                            if protein_correct == False:
                                deletion_in_wt_amp_cds_noPL_Flag = True #hitch a ride
                            
                            #no need to check deletion in SNP, because the SNP translation will take care it this
                            #produce output
                            if deletion_in_wt_amp_cds_noPL_Flag==True:
                                #print("\twt allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write(f"\tmutant protein + {payload_genotype} SNP\t")
                            elif deletion_in_wt_amp_cds_noPL_Flag==False:
                                #print("\twt allele (wt protein + HDR payload)", end="\n")   
                                writehandle.write(f"\twt protein + {payload_genotype} SNP\t")      
                        elif (int(n_deleted)==0 and int(n_inserted)!=0):  #insertion in the read, no deletion
                            #check if insertion are in wt cds
                            insertion_in_wt_amp_cds_noPL_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = wt_amp_cds_coords_PosRef)["insertion_in_cds_flag"]
                            #check if insertion are in payload cds
                            insertion_in_payload_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = [payload_coord])["insertion_in_cds_flag"]
                            
                            #trim off inserted seq (in both reads and ref)
                            reg_gap_windows = cal_gap_win(seq = ref)
                            ref_trimmed = ''.join(ref[idx] for idx in range(len(ref)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))
                            read_trimmed = ''.join(read[idx] for idx in range(len(read)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))    
                            
                            #check translation 
                            protein_correct = check_protein_and_payload(read_trimmed, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["protein_correct"]
                            payload_genotype = check_protein_and_payload(read_trimmed, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["payload_genotype"]
                            
                            #summarize flag
                            mut_protein_summary_flag = any([insertion_in_wt_amp_cds_noPL_Flag,
                                                        (not protein_correct)])        
                            
                            if insertion_in_payload_Flag:
                                payload_genotype = "mutant"
                            
                            #produce output
                            if mut_protein_summary_flag==True:
                                #print("\twt allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write(f"\tmutant protein + {payload_genotype} SNP\t")
                            elif mut_protein_summary_flag==False:
                                #print("\twt allele (wt protein + HDR payload)", end="\n")
                                writehandle.write(f"\twt protein + {payload_genotype} SNP\t") 
                        elif(int(n_deleted)!=0 and int(n_inserted)!=0): #both insertion and deletion in the read
                            #check if insertion are in wt cds
                            insertion_in_wt_amp_cds_noPL_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = wt_amp_cds_coords_PosRef)["insertion_in_cds_flag"]
                            #check if insertion are in payload cds
                            insertion_in_payload_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = [payload_coord])["insertion_in_cds_flag"]            

                            #trim off inserted seq (in both reads and ref)
                            reg_gap_windows = cal_gap_win(seq = ref)
                            ref_trimmed = ''.join(ref[idx] for idx in range(len(ref)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))
                            read_trimmed = ''.join(read[idx] for idx in range(len(read)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))                     
                        
                            #check if deletions are in wt cds
                            deletion_in_wt_amp_cds_noPL_Flag = any([check_deletion_in_cds(read = read_trimmed, amp_cds_coords = wt_amp_cds_coords_PosRef)["deletion_in_cds_flag"],
                                                                    check_deletion_overlap_cds(read = read_trimmed, amp_cds_coords = wt_amp_cds_coords_PosRef)["deletion_overlap_cds_flag"]])
                            #check translation 
                            protein_correct = check_protein_and_payload(read_trimmed, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["protein_correct"]
                            payload_genotype = check_protein_and_payload(read_trimmed, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation)["payload_genotype"]
                            
                            #summarize flag
                            indel_in_wt_amp_cds_noPL_Flag = any([insertion_in_wt_amp_cds_noPL_Flag,
                                                                deletion_in_wt_amp_cds_noPL_Flag, 
                                                                (not protein_correct)])
                    
                            #produce output
                            if indel_in_wt_amp_cds_noPL_Flag==True:
                                #print("\twt allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write(f"\tmutant protein + {payload_genotype} SNP\t")
                            elif indel_in_wt_amp_cds_noPL_Flag==False:
                                #print("\twt allele (wt protein + HDR payload)", end="\n")
                                writehandle.write(f"\twt protein + {payload_genotype} SNP\t") 
                    #write percent identity and num_of_mismatches
                    writehandle.write(f"{percent_identity}\t{num_of_mismatches}\n")

        #write the allele frequency table to a new zip file
        new_zip_path = "_".join([zip_path.rstrip(r'.zip'),"genotype.zip"])
        if os.path.exists(new_zip_path):
            os.remove(new_zip_path)
        with zipfile.ZipFile(new_zip_path,compression=zipfile.ZIP_BZIP2, compresslevel=1, mode = 'x') as myzip2:
            wd = os.getcwd() #save current working dir
            os.chdir(zip_dir) #change to the the folder containg the file to be zipped
            myzip2.write(os.path.join('Alleles_frequency_table_genotype.txt'))
            os.chdir(wd) #change to the saved working dir
        if os.path.exists(os.path.join(zip_dir,'Alleles_frequency_table_genotype.txt')): #remove unzipped raw output
            os.remove(os.path.join(zip_dir,'Alleles_frequency_table_genotype.txt'))    
        log.info(f"Done processing alignments, saved the results in zip file: {new_zip_path}")

        ############################
        #calculate allele frequency#
        ############################
        #calculate allele frequency
        log.info(f"Calculating allele frequency")
        genotype_freq = {"wt_allele":0,
                        "HDR_perfect":0,
                        "wtProt_hdrSNP":0,
                        "mutProt_hdrSNP":0,
                        "mutProt_mutSNP":0,
                        "wtProt_mutSNP":0,
                        "mutProt_wtSNP":0,
                        "wtProt_wtSNP":0}

        with zipfile.ZipFile(new_zip_path, 'r') as myzip:
            d_total_num_reads = 0
            d_num_reads_post_read_group_filter = 0
            d_num_alleles = 0
            perc_identity_list = []
            num_of_mismatches_list = []
            with myzip.open('Alleles_frequency_table_genotype.txt') as filehandle:
                next(filehandle) #skip header
                #print(f"Reference_Name\tRead_Status\tn_Reads\tperc_Reads\tn_deleted\tn_inserted\tn_mutated")
                for line in filehandle:
                    line_deco = line.decode()

                    fields = line_deco.rstrip().split("\t")
                    read = fields[0]
                    ref = fields[1]
                    Reference_Name = fields[2]
                    Read_Status = fields[3]
                    n_deleted = fields[4]
                    n_inserted = fields[5]
                    n_mutated = fields[6]
                    n_Reads = fields[7]
                    perc_Reads = fields[8]
                    genotype = fields[9]
                    perc_identity = fields[10]
                    num_of_mismatches = fields[11]

                    d_total_num_reads += int(n_Reads)
                    d_num_alleles += 1
                    if int(n_Reads) >= config['min_reads_for_genotype']:
                        d_num_reads_post_read_group_filter += int(n_Reads)
                        

                    if Reference_Name == "HDR": #only consider HDR for weighted average of percent identity and num_of_mismatches
                        perc_identity_list.append(float(perc_identity))
                        num_of_mismatches_list.append(int(num_of_mismatches))

                    #print(f"{Reference_Name}\t{n_Reads}\t{perc_Reads}\t{n_deleted}\t{n_inserted}\t{n_mutated}\t{genotype}", end='\n')
                    if genotype == "perfect HDR edit":
                        genotype_freq["HDR_perfect"] += int(n_Reads)
                    elif genotype == "wt protein + HDR SNP":
                        genotype_freq["wtProt_hdrSNP"] += int(n_Reads)
                    elif genotype == "wt protein + mutant SNP":
                        genotype_freq["wtProt_mutSNP"] += int(n_Reads)  
                    elif genotype == "wt protein + wt SNP":
                        genotype_freq["wtProt_wtSNP"] += int(n_Reads)                  
                    elif genotype == "mutant protein + HDR SNP":
                        genotype_freq["mutProt_hdrSNP"] += int(n_Reads)                  
                    elif genotype == "mutant protein + mutant SNP":
                        genotype_freq["mutProt_mutSNP"] += int(n_Reads)                  
                    elif genotype == "mutant protein + wt SNP":
                        genotype_freq["mutProt_wtSNP"] += int(n_Reads)      
                    elif genotype == "wt allele":
                        genotype_freq["wt_allele"] += int(n_Reads)  

            # compute the percentage of each genotype
            for key in genotype_freq:
                if d_total_num_reads == 0:
                    genotype_freq[key] = 0
                else:   
                    genotype_freq[key] = genotype_freq[key] / d_total_num_reads * 100

        # get reads statistics from CRISPResso
        with open(os.path.join(zip_dir,'CRISPResso_mapping_statistics.txt'),'r') as readstat_file:
            header_line = next(readstat_file)
            second_line = next(readstat_file)
            fields = second_line.rstrip().split("\t")
            num_reads_in_fq = fields[0]
            num_reads_post_QC = fields[1]
            num_reads_aligned = fields[2]
            num_computed_aln = fields[3]


        #write genotype to file
        with open(os.path.join(zip_dir,'genotype_frequency.csv'),'w') as writehandle:
            #THIS IS THE OLD OUTPUTR\
            #writehandle.write("Sample,wt_allele,HDR_perfect,wtProt_wtSNP,wtProt_hdrSNP,mutProt_wtSNP,mutProt_hdrSNP,mutProt_mutSNP,wtProt_mutSNP\n")
            # writehandle.write(f"{float(genotype_freq['wt_allele']):.4f}%,"
            #                 f"{float(genotype_freq['HDR_perfect']):.4f}%,"
            #                 f"{float(genotype_freq['wtProt_wtSNP']):.4f}%,"  
            #                 f"{float(genotype_freq['wtProt_hdrSNP']):.4f}%,"
            #                 f"{float(genotype_freq['mutProt_wtSNP']):.4f}%,"               
            #                 f"{float(genotype_freq['mutProt_hdrSNP']):.4f}%,"
            #                 f"{float(genotype_freq['mutProt_mutSNP']):.4f}%,"
            #                 f"{float(genotype_freq['wtProt_mutSNP']):.4f}%\n")
            writehandle.write("Sample,wt_allele,HDR_perfect,wtProt_hdrSNP,mutProt_hdrSNP,wtProt_mutSNP,mutProt_wtSNP,mutProt_mutSNP,wtProt_wtSNP,num_clean_reads,weighted_avg_perc_identity,weighted_avg_num_of_mismatches,num_reads_in_fq,num_reads_post_QC,num_reads_aligned,[dg]num_reads_post_read_group_filter,num_alleles,[dg]num_alleles\n")
            writehandle.write(f"{float(genotype_freq['wt_allele']):.4f}%,"
                            f"{float(genotype_freq['HDR_perfect']):.4f}%,"
                            f"{float(genotype_freq['wtProt_hdrSNP']):.4f}%,"
                            f"{float(genotype_freq['mutProt_hdrSNP']):.4f}%,"
                            f"{float(genotype_freq['wtProt_mutSNP']):.4f}%,"
                            f"{float(genotype_freq['mutProt_wtSNP']):.4f}%,"
                            f"{float(genotype_freq['mutProt_mutSNP']):.4f}%,"
                            f"{float(genotype_freq['wtProt_wtSNP']):.4f}%,"
                            f"{compute_weighted_average(perc_identity_list, num_of_mismatches_list):.4f},"
                            f"{compute_weighted_average(num_of_mismatches_list, num_of_mismatches_list):.4f},"
                            f"{num_reads_in_fq},"
                            f"{num_reads_post_QC},"
                            f"{num_reads_aligned},"
                            f"{d_num_reads_post_read_group_filter},"
                            f"{num_computed_aln},"
                            f"{d_num_alleles}\n"
                            )
                            

    except Exception  as e:
        print("Unexpected error:", str(sys.exc_info()))
        print("additional information:", e)
        PrintException()











    
if __name__ == "__main__": main()    
