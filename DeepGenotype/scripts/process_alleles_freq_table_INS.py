#load packages
import argparse, sys, os, re, zipfile,  logging, warnings
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
    parser.add_argument('--payload_block_index', default=1, type=int, help='1-based index of the insertion block to treat as the payload if multiple are found')
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
        payload_block_index = config['payload_block_index']


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


        #get payload coordinates in HDR amp 
        payload_seq = ""
        payload_len = 0
        payload_coord_in_HDR_amp = [] # *regardless of strand*, the start and end cooresponds to the HDR amplicon
        payload_strand = "+"
        if mytranscript.annotations['transcript_strand'] == -1:
            payload_strand = "-"
            
        # Align WT to HDR and look for insertion blocks (gaps in wt_aln)
        aln = align.localms(wt_amp,HDR_amp,2, -1, -.5, -.1)
        wt_aln = format_alignment(*aln[0]).split("\n")[0]

        # Gather all gap blocks in the WT portion
        insertion_blocks = list(re.finditer(r"-+", wt_aln))

        if len(insertion_blocks) == 0:
            # No insertion blocks found
            log.warning("No insertion blocks found between WT and HDR. Payload information will be empty.")
        elif len(insertion_blocks) == 1:
            # Only one block found — proceed as before
            m = insertion_blocks[0]
            payload_len = m.span()[1] - m.span()[0]
            payload_seq = HDR_amp[m.span()[0] : m.span()[1]]
            payload_coord_in_HDR_amp = [m.span()[0], m.span()[1]]
        else:
            # Multiple insertion blocks found — user must choose
            log.warning(f"Multiple insertion blocks found ({len(insertion_blocks)}). Using block "
                        f"#{payload_block_index} as the designated payload.")
            # Check that the user-specified index is valid
            if payload_block_index < 1 or payload_block_index > len(insertion_blocks):
                raise ValueError(f"Invalid payload_block_index={payload_block_index}. Must be between "
                                f"1 and {len(insertion_blocks)}.")
            # Select the block corresponding to payload_block_index (1-based)
            chosen_block = insertion_blocks[payload_block_index - 1]
            payload_len = chosen_block.span()[1] - chosen_block.span()[0]
            payload_seq = HDR_amp[chosen_block.span()[0] : chosen_block.span()[1]]
            payload_coord_in_HDR_amp = [chosen_block.span()[0], chosen_block.span()[1]]

        ## old code, replaced by the new code above
        # for m in re.finditer("-+",wt_aln):
        #     if m:
        #         payload_len = m.span()[1]-m.span()[0]
        #         payload_seq = HDR_amp[m.span()[0]:m.span()[1]]
        #         payload_coord_in_HDR_amp = [m.span()[0],m.span()[1]]

        #print(f"\npayload coord in HDR amp: {payload_coord_in_HDR_amp}")
        #print(f"payload length: {payload_len}") 
        #print(f"payload strand: {payload_strand}") 
        #print(f"payload seq: {payload_seq}") 

            
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
            
        #translate payload
        payload_translation = translate_payload(seq=HDR_amp, start = payload_coord_in_HDR_amp[0], end = payload_coord_in_HDR_amp[1], strand = payload_strand, frame = 1)

        #print(f"\nwt amp translation {wt_translation}")
        #print(f"HDR_translation {HDR_translation}")
        #print(f"payload_translation {payload_translation}")


        payload_coord = payload_coord_in_HDR_amp # map legacy name to current name, #no need to convert coordinates, because they are already in reference to the + strand HDR amp
        HDR_amp_cds_coords_PosRef = get_amp_abs_coords(HDR_amp, HDR_amp_cds_coords) # convert coordinate so that they are all in reference to the + strand HDR amp, HDR_amp_cds_coords_PosRef canbe only used to check indels locations and *not translation*
        wt_amp_cds_coords_PosRef = get_amp_abs_coords(wt_amp, wt_amp_cds_coords) # convert coordinate so that they are all in reference to the + strand wt amp, wt_amp_cds_coords_PosRef canbe only used to check indels locations and *not translation*

        ########################
        #process the alignments#
        ########################
        # read the allele frequency table (from CRISPResso) and write an updated table to Alleles_frequency_table_genotype.txt
        with zipfile.ZipFile(zip_path, 'r') as myzip, open(os.path.join(zip_dir,'Alleles_frequency_table_genotype.txt'),'w') as writehandle:
            with myzip.open('Alleles_frequency_table.txt') as filehandle:
                next(filehandle) #skip header
                #print(f"Reference_Name\tRead_Status\tn_Reads\tperc_Reads\tn_deleted\tn_inserted\tn_mutated")
                writehandle.write(f"Aligned_Sequence	Reference_Sequence\tReference_Name\tRead_Status\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\tGenotype\n")
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

                    writehandle.write(line_deco.rstrip())
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
                        if read == ref:
                            #print("\tHDR allele (perfect HDR edit)", end="\n")
                            writehandle.write("\tHDR allele (perfect HDR edit)\t")
                        ##########################
                        #mutations only, no indel#
                        ##########################
                        elif int(n_deleted)==0 and int(n_inserted)==0: 
                            #print(f"\n{read}\n{ref}", end='\n')
                            payload_correct = check_protein_correctness(read,HDR_amp_cds_coords,HDR_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation)["payload_correct"]
                            protein_correct = check_protein_correctness(read,HDR_amp_cds_coords,HDR_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation)["protein_correct"]

                            #produce output
                            if (protein_correct==True and payload_correct==True):
                                #print("\tHDR allele (wt protein + correct payload)", end="\n")
                                writehandle.write("\tHDR allele (wt protein + correct payload)\t")
                            elif (protein_correct==True and payload_correct==False):
                                #print("\tHDR allele (wt protein + mutant payload)", end="\n")
                                writehandle.write("\tHDR allele (wt protein + mutant payload)\t")
                            elif (protein_correct==False and payload_correct==True):
                                #print("\tHDR allele (mutant protein + correct payload)", end="\n")
                                writehandle.write("\tHDR allele (mutant protein + correct payload)\t")    
                            elif (protein_correct==False and payload_correct==False):
                                #print("\tHDR allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write("\tHDR allele (mutant protein + mutant payload)\t")    
                        ########################
                        #indels in HDR allele  #
                        ########################
                        elif (int(n_deleted)!=0 and int(n_inserted)==0):  #deletion in the read，no insertion
                            #check if deletions are in wt cds
                            deletion_in_HDR_amp_cds_noPL_Flag = any([check_deletion_in_cds(read = read, amp_cds_coords = HDR_amp_cds_coords_PosRef)["deletion_in_cds_flag"],
                                                                    check_deletion_overlap_cds(read = read, amp_cds_coords = HDR_amp_cds_coords_PosRef)["deletion_overlap_cds_flag"]])
                            #check if deletions are in payload cds
                            deletion_in_payload_Flag = any([check_deletion_in_cds(read = read, amp_cds_coords = [payload_coord])["deletion_in_cds_flag"],
                                                        check_deletion_overlap_cds(read = read, amp_cds_coords = [payload_coord])["deletion_overlap_cds_flag"]])
                            
                            payload_correct = check_protein_correctness(read,HDR_amp_cds_coords,HDR_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation)["payload_correct"]
                            protein_correct = check_protein_correctness(read,HDR_amp_cds_coords,HDR_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation)["protein_correct"]
                            
                            if payload_correct == False:
                                deletion_in_payload_Flag = True #hitch a ride with the deletion flag
                            if protein_correct == False:
                                deletion_in_HDR_amp_cds_noPL_Flag = True #hitch a ride with the deletion flag

                            #produce output
                            if (deletion_in_HDR_amp_cds_noPL_Flag==True and deletion_in_payload_Flag==True):
                                #print("\tHDR allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write("\tHDR allele (mutant protein + mutant payload)\t")
                            elif (deletion_in_HDR_amp_cds_noPL_Flag==True and deletion_in_payload_Flag==False):
                                #print("\tHDR allele (mutant protein + correct payload)", end="\n")
                                writehandle.write("\tHDR allele (mutant protein + correct payload)\t")
                            elif (deletion_in_HDR_amp_cds_noPL_Flag==False and deletion_in_payload_Flag==True):
                                #print("\tHDR allele (wt protein + mutant payload)", end="\n")
                                writehandle.write("\tHDR allele (wt protein + mutant payload)\t")    
                            elif (deletion_in_HDR_amp_cds_noPL_Flag==False and deletion_in_payload_Flag==False):
                                #print("\tHDR allele (wt protein + correct payload1)", end="\n")   
                                writehandle.write("\tHDR allele (wt protein + correct payload)\t")      
                        elif (int(n_deleted)==0 and int(n_inserted)!=0):  #insertion in the read, no deletion
                            #check if insertion are in wt cds
                            insertion_in_HDR_amp_cds_noPL_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = HDR_amp_cds_coords_PosRef)["insertion_in_cds_flag"]
                            #check if insertion are in payload cds
                            insertion_in_payload_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = [payload_coord])["insertion_in_cds_flag"]
                            #produce output
                            
                            #trim off inserted seq (in both reads and ref)
                            reg_gap_windows = cal_gap_win(seq = ref)
                            ref_trimmed = ''.join(ref[idx] for idx in range(len(ref)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))
                            read_trimmed = ''.join(read[idx] for idx in range(len(read)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))   
                            
                            payload_correct = check_protein_correctness(read_trimmed,HDR_amp_cds_coords,HDR_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation)["payload_correct"]
                            protein_correct = check_protein_correctness(read_trimmed,HDR_amp_cds_coords,HDR_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation)["protein_correct"]
                            
                            if payload_correct == False:
                                insertion_in_payload_Flag = True #hitch a ride with the deletion flag
                            if protein_correct == False:
                                insertion_in_HDR_amp_cds_noPL_Flag = True #hitch a ride with the deletion flag
                                
                            if (insertion_in_HDR_amp_cds_noPL_Flag==True and insertion_in_payload_Flag==True):
                                #print("\tHDR allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write("\tHDR allele (mutant protein + mutant payload)\t")
                            elif (insertion_in_HDR_amp_cds_noPL_Flag==True and insertion_in_payload_Flag==False):
                                #print("\tHDR allele (mutant protein + correct payload)", end="\n")
                                writehandle.write("\tHDR allele (mutant protein + correct payload)\t")
                            elif (insertion_in_HDR_amp_cds_noPL_Flag==False and insertion_in_payload_Flag==True):
                                #print("\tHDR allele (wt protein + mutant payload)", end="\n")
                                writehandle.write("\tHDR allele (wt protein + mutant payload)\t")
                            elif (insertion_in_HDR_amp_cds_noPL_Flag==False and insertion_in_payload_Flag==False):
                                #print("\tHDR allele (wt protein + correct payload)", end="\n") 
                                writehandle.write("\tHDR allele (wt protein + correct payload)\t") 
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
                                                                    check_deletion_overlap_cds(read = read_trimmed, amp_cds_coords = HDR_amp_cds_coords_PosRef)["deletion_overlap_cds_flag"]
                                                                    ])
                            
                            #check if deletions are in payload cds
                            deletion_in_payload_Flag = any([check_deletion_in_cds(read = read_trimmed, amp_cds_coords = [payload_coord])["deletion_in_cds_flag"],
                                                        check_deletion_overlap_cds(read = read_trimmed, amp_cds_coords = [payload_coord])["deletion_overlap_cds_flag"]
                                                        ])

                            payload_correct = check_protein_correctness(read_trimmed,HDR_amp_cds_coords,HDR_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation)["payload_correct"]
                            protein_correct = check_protein_correctness(read_trimmed,HDR_amp_cds_coords,HDR_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation)["protein_correct"]
                            
                            #summarize flag
                            indel_in_HDR_amp_cds_noPL_Flag = any([insertion_in_HDR_amp_cds_noPL_Flag,
                                                                deletion_in_HDR_amp_cds_noPL_Flag,
                                                                (not protein_correct)])  #hitch a ride with the indel flag
                            indel_in_payload_Flag = any([insertion_in_payload_Flag,
                                                        deletion_in_payload_Flag,
                                                        not payload_correct])   #hitch a ride with the indel flag   
                            #produce output
                            if (indel_in_HDR_amp_cds_noPL_Flag==True and indel_in_payload_Flag==True):
                                #print("\tHDR allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write("\tHDR allele (mutant protein + mutant payload)\t")
                            elif (indel_in_HDR_amp_cds_noPL_Flag==True and indel_in_payload_Flag==False):
                                #print("\tHDR allele (mutant protein + correct payload)", end="\n")
                                writehandle.write("\tHDR allele (mutant protein + correct payload)\t")    
                            elif (indel_in_HDR_amp_cds_noPL_Flag==False and indel_in_payload_Flag==True):
                                #print("\tHDR allele (wt protein + mutant payload)", end="\n")
                                writehandle.write("\tHDR allele (wt protein + mutant payload)\t")    
                            elif (indel_in_HDR_amp_cds_noPL_Flag==False and indel_in_payload_Flag==False):
                                #print("\tHDR allele (wt protein + correct payload)", end="\n")
                                writehandle.write("\tHDR allele (wt protein + correct payload)\t")    

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
                            #print(f"\n{read}\n{ref}", end='\n')
                            protein_correct = check_protein_correctness(read,wt_amp_cds_coords,wt_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation)["protein_correct"]
                            #produce output
                            if (protein_correct==True):
                                #print("\twt allele (wt protein + no payload)", end="\n")
                                writehandle.write("\twt allele (wt protein)\t")
                            elif (protein_correct==False):
                                #print("\twt allele (mutant protein + no payload)", end="\n")
                                writehandle.write("\twt allele (mutant protein)\t")
                        #######################
                        #indels in wt allele  #
                        #######################
                        elif (int(n_deleted)!=0 and int(n_inserted)==0):  #deletion in the read，no insertion
                            #check if deletions are in cds
                            deletion_in_cds_flag = any([check_deletion_in_cds(read = read, amp_cds_coords = wt_amp_cds_coords_PosRef)["deletion_in_cds_flag"],
                                                        check_deletion_overlap_cds(read = read, amp_cds_coords = wt_amp_cds_coords_PosRef)["deletion_overlap_cds_flag"]
                                                    ])
                            
                            protein_correct = check_protein_correctness(read,wt_amp_cds_coords,wt_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation)["protein_correct"]
                            if protein_correct == False:
                                deletion_in_cds_flag= True #hitch a ride with the del flag
                        
                            if deletion_in_cds_flag == True:
                                #print("\twt allele (mutant protein)", end="\n")
                                writehandle.write("\twt allele (mutant protein)\t")
                            else:
                                #print("\twt allele (wt protein)", end="\n")
                                writehandle.write("\twt allele (wt protein)\t")
                        elif (int(n_deleted)==0 and int(n_inserted)!=0):  #insertion in the read, no deletion
                            #check if insertion are in cds
                            insertion_in_cds_flag = check_insertion_in_cds(ref = ref, amp_cds_coords = wt_amp_cds_coords_PosRef)["insertion_in_cds_flag"]
                            
                            #trim off inserted seq (in both reads and ref)
                            reg_gap_windows = cal_gap_win(seq = ref)
                            ref_trimmed = ''.join(ref[idx] for idx in range(len(ref)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))
                            read_trimmed = ''.join(read[idx] for idx in range(len(read)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))    
                            
                            protein_correct = check_protein_correctness(read_trimmed,wt_amp_cds_coords,wt_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation)["protein_correct"]
                            
                            if protein_correct == False:
                                insertion_in_cds_flag= True #hitch a ride with the del flag
                                
                            if insertion_in_cds_flag == True:
                                #print("\twt allele (mutant protein)", end="\n")
                                writehandle.write("\twt allele (mutant protein)\t")
                            else:
                                #print("\twt allele (wt protein)", end="\n")
                                writehandle.write("\twt allele (wt protein)\t")
                        elif(int(n_deleted)!=0 and int(n_inserted)!=0): #both insertion and deletion in the read
                            #check if insertion are in cds
                            insertion_in_cds_flag = check_insertion_in_cds(ref = ref, amp_cds_coords = wt_amp_cds_coords_PosRef)["insertion_in_cds_flag"]
                            #trim off inserted seq (in both reads and ref)
                            reg_gap_windows = cal_gap_win(seq = ref)
                            ref_trimmed = ''.join(ref[idx] for idx in range(len(ref)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))
                            read_trimmed = ''.join(read[idx] for idx in range(len(read)) if not any([idx in range(st,en) for st,en in reg_gap_windows]))    
                            #check if deletions are in cds (using the trimmed sequences, this is IMPORTANT!)
                            deletion_in_cds_flag = any([check_deletion_in_cds(read = read_trimmed, amp_cds_coords = wt_amp_cds_coords_PosRef)["deletion_in_cds_flag"],
                                                        check_deletion_overlap_cds(read = read_trimmed, amp_cds_coords = wt_amp_cds_coords_PosRef)["deletion_overlap_cds_flag"]
                                                    ])
                            #read_gap_loc = check_deletion_in_cds(read = read_trimmed, amp_cds_coords = wt_amp_cds_coords)["read_gap_loc"]
                            protein_correct = check_protein_correctness(read_trimmed,wt_amp_cds_coords,wt_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation)["protein_correct"]
                            
                            if protein_correct == False:
                                deletion_in_cds_flag= True #hitch a ride with the del flag
                            
                            if any([insertion_in_cds_flag, deletion_in_cds_flag]):
                                #print("\twt allele (mutant protein)", end="\n")
                                writehandle.write("\twt allele (mutant protein)\t")
                            else:
                                #print("\twt allele (wt protein)", end="\n")
                                writehandle.write("\twt allele (wt protein)\t")

                    #write percent identity and num_of_mismatches
                    writehandle.write(f"{percent_identity}\t{num_of_mismatches}\n")

        #write the allele frequency table to a new zip file (genotype.zip)
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
        #calculate allele frequency and write to genotype_frequency.csv
        log.info(f"Calculating allele frequency")
        genotype_freq = {"wt_allele":0,
                        "HDR_perfect":0,
                        "wtProt_noPL":0,
                        "wtProt_OKPL":0,
                        "mutProt_noPL":0,
                        "mutProt_OKPL":0,
                        "mutProt_mutPL":0,
                        "wtProt_mutPL":0}
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
                    if genotype == "HDR allele (perfect HDR edit)":
                        genotype_freq["HDR_perfect"] += int(n_Reads)
                    elif genotype == "HDR allele (wt protein + correct payload)":
                        genotype_freq["wtProt_OKPL"] += int(n_Reads)
                    elif genotype == "HDR allele (wt protein + mutant payload)":
                        genotype_freq["wtProt_mutPL"] += int(n_Reads)        
                    elif genotype == "HDR allele (mutant protein + correct payload)":
                        genotype_freq["mutProt_OKPL"] += int(n_Reads)                  
                    elif genotype == "HDR allele (mutant protein + mutant payload)":
                        genotype_freq["mutProt_mutPL"] += int(n_Reads)                  
                    elif genotype == "wt allele":
                        genotype_freq["wt_allele"] += int(n_Reads)  
                    elif genotype == "wt allele (wt protein)":
                        genotype_freq["wtProt_noPL"] += int(n_Reads)
                    elif genotype == "wt allele (mutant protein)":
                        genotype_freq["mutProt_noPL"] += int(n_Reads)    
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
            writehandle.write(f"wt_allele,"
                            f"HDR_perfect,"
                            f"wtProt_noPL,"
                            f"wtProt_OKPL,"
                            f"mutProt_noPL,"
                            f"mutProt_OKPL,"
                            f"mutProt_mutPL,"
                            f"wtProt_mutPL,"
                            f",num_reads_in_fq,num_reads_post_QC,num_reads_aligned,[dg]num_reads_post_read_group_filter,num_alleles,[dg]num_alleles\n")
            writehandle.write(f"{float(genotype_freq['wt_allele']):.4f}%,"
                            f"{float(genotype_freq['HDR_perfect']):.4f}%,"
                            f"{float(genotype_freq['wtProt_noPL']):.4f}%,"
                            f"{float(genotype_freq['wtProt_OKPL']):.4f}%,"
                            f"{float(genotype_freq['mutProt_noPL']):.4f}%,"
                            f"{float(genotype_freq['mutProt_OKPL']):.4f}%,"
                            f"{float(genotype_freq['mutProt_mutPL']):.4f}%,"
                            f"{float(genotype_freq['wtProt_mutPL']):.4f}%,"
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
