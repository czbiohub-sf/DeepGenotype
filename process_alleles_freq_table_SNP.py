#load packages
import argparse
import sys
import linecache
import pandas as pd
import os
import re
import zipfile
import requests
import logging
from Bio.pairwise2 import align
from Bio.pairwise2 import format_alignment
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import warnings
warnings.filterwarnings('ignore')
log = logging.getLogger("Process_alleles_freq_table.py")
log.setLevel(logging.CRITICAL) #set the level of warning displayed

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():
    parser= MyParser(description='This script does processes the allele freq. table output from CRISPResso2')
    parser.add_argument('--path', default="", type=str, help='path to the allele freq. table output zip file')
    parser.add_argument('--allele_freq_file', default="", type=str, help='the name of the freq. table output zip file')
    parser.add_argument('--wt_amp', default="", type=str, help='sequence of the wt amplicon')
    parser.add_argument('--HDR_amp', default="", type=str, help='sequence of the HDR amplicon')
    parser.add_argument('--ENST_ID', default="", type=str, help='ENST ID')
    parser.add_argument('--SNP_block_of_interest', default=1, type=int, help='if the SNP is downstream of re-cut-preventing SNP, set this to 2')

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
        SNP_block_of_interest = config['SNP_block_of_interest']
        SNP_idx = SNP_block_of_interest - 1

        #supporting functions
        def revcom(seq):
            return str(Seq(seq).reverse_complement())

        def adjust_frame(seq, frame): #adjusts the frame 
            if frame == 1:
                return(seq)
            elif frame == 2:
                return(seq[1:])
            elif frame == 3:
                return(seq[2:])
            else:
                return("invalid frame!")
            
        def translate_amp(seq, start, end, strand, frame): # for amplicons, the start, end refers to the revcom if strand = "-"
            #seq = seq[start:end].replace("-", "")
            if strand == "-":
                seq = revcom(seq) #revcom
                seq = seq[start:end]
                coding_dna = adjust_frame(seq,frame) #adjust frame
                coding_dna=coding_dna.replace("-", "")
                trans = str(Seq(coding_dna).translate())    #translate
                return(trans)      
            else:
                seq = seq[start:end]
                coding_dna = adjust_frame(seq,frame) #adjust frame
                coding_dna = coding_dna.replace("-", "")
                trans = str(Seq(coding_dna).translate())    #translate
                return(trans)      

        def translate_payload(seq, start, end, strand, frame): # for payloads, the start, end refers to the seq (HDR amplicon) , not the revcom 
            #seq = seq[start:end].replace("-", "")
            if strand == "-":
                seq = seq[start:end]
                seq = revcom(seq) #revcom
                coding_dna = adjust_frame(seq,frame) #adjust frame
                coding_dna=coding_dna.replace("-", "")
                trans = str(Seq(coding_dna).translate())    #translate
                return(trans)      
            else:
                seq = seq[start:end]
                coding_dna = adjust_frame(seq,frame) #adjust frame
                coding_dna = coding_dna.replace("-", "")
                trans = str(Seq(coding_dna).translate())    #translate
                return(trans)  
            
        def infer_strand(tag_for,tag_rev,HDR_amp): #use the provided tag sequence to infer the strand and frame of the coding sequence
            mFor = re.finditer(tag_for, HDR_amp)
            mRev = re.finditer(tag_rev, HDR_amp)  
            if sum(1 for _ in mFor) == 1: #found tag in sense strand
                for m in re.finditer(tag_for, HDR_amp):
                    st,en = m.span()
                    frame = st%3+1
                    #print(f"{row['gene_name']} {st} {en} + {frame}")
                    return(["+",frame, st, en])
            else:
                for m in re.finditer(tag_rev, HDR_amp):
                    st,en = m.span()
                    frame = (len(HDR_amp)-en)%3+1
                    #print(f"{row['gene_name']} {st} {en} - {frame}")
                    return(["-",frame, st, en])    
                
        #get cds coordinate in amplicons
        #align each cds to amplicon and mark cds segments
        def get_cds_coord_in_amplicon(whole_cds, amplicon):
            """
            Aligns CDS to amplicon, and get the coordinates of coding part of amplicon
        
            Notes:
            whole_cds is always +1 frame
            
            """
            log.info(f"Searching for cds segments in amplicon...")

            #alignment
            aln = align.localms(amplicon, whole_cds, 0.1, -10, -0.2, -.1)
            aln_rc = align.localms(revcom(amplicon), whole_cds, 0.1, -10, -0.2, -.1)

            #get coding sequence coordinates
            cds_coord_in_amp = []
            if aln[0].score >= aln_rc[0].score:  # amplicon matched to sense CDS
                #print("amp matched to CDS")
                #print(format_alignment(*aln[0]))
                coords = get_st_en_from_aln(aln) #get start and end of amplicon match, and cds match
                #print(f"start and end of amplicon match, and cds match\n{coords['seq1_st']}-{coords['seq1_en']} {coords['seq2_st']}-{coords['seq2_en']}")

                #get frame
                frame = int(coords['seq2_st'])%3
                if frame ==0:
                    frame =1
                elif frame == 1: 
                    frame = 3
                elif frame == 2: 
                    frame = 2

                #print(f"strand= '-'\n{frame}")
                #return coordinates
                cds_coord_in_amp.append((coords['seq1_st'],coords['seq1_en'],"+",frame))

            else: # revcom amp matched to CDS
                #print("revcom amp matched to CDS")
                #print(format_alignment(*aln_rc[0]))
                coords = get_st_en_from_aln(aln_rc)
                #print(f"start and end of amplicon match, and cds match\n{coords['seq1_st']}-{coords['seq1_en']} {coords['seq2_st']}-{coords['seq2_en']}")

                #get frame
                frame = int(coords['seq2_st'])%3
                if frame ==0:
                    frame =1
                elif frame == 1: 
                    frame = 3
                elif frame == 2: 
                    frame = 2

                #print(f"strand= '-'\nframe={frame}")
                #return coordinates
                cds_coord_in_amp.append((coords['seq1_st'],coords['seq1_en'],"-",frame))
                    
            log.info(f"Done searching cds segments in the amplicon")
            return(cds_coord_in_amp)

        def get_st_en_from_aln(aln):
            seq1_start = int(format_alignment(*aln[0]).split("\n")[0].split()[0])
            seq1_match = format_alignment(*aln[0]).split("\n")[0].split()[1]
            seq1_match = re.sub(r'-+','',seq1_match)
            seq1_st = seq1_start - 1
            seq1_en = seq1_start + len(seq1_match) -1
            
            seq2_start = int(format_alignment(*aln[0]).split("\n")[2].split()[0])
            seq2_match = format_alignment(*aln[0]).split("\n")[2].split()[1]
            seq2_match = re.sub(r'-+','',seq2_match)
            seq2_st = seq2_start - 1
            seq2_en = seq2_start + len(seq2_match) -1    

            return({"seq1_st":seq1_st,
                    "seq1_en":seq1_en,
                    "seq2_st":seq2_st,
                    "seq2_en":seq2_en})

        #get cds seq
        def get_cds_seq_in_transcript(mytranscript):
            '''
            return a list of tuples, each tuple is a set of start,end of a segment of cds
            '''
            wholeSeq = str(mytranscript.seq)
            wholeSeq_rc = str(Seq(wholeSeq).reverse_complement())
            cds_seqs = []
            for feat in mytranscript.features:
                if feat.type == "cds":
                    st = feat.location.start.position
                    en = feat.location.end.position
                    if (feat.strand==-1): # neg strand
                        cds_seqs.append(str(Seq(wholeSeq_rc[st:en]).reverse_complement())) # the coord are respective to the revcom of the retrieved seq (weird)
                    else:  # pos strand
                        cds_seqs.append(wholeSeq[st:en]) 
            return(cds_seqs)  

        def cal_gap_positions(seq):
            gap = re.compile("-+")
            #find gaps
            gaps=[]
            for m in gap.finditer(seq):
                gaps.append(m.span())
            #go through all gaps and project each gap to a point in the ungapped seq 
            gap_positions = []
            gap_accu_len = 0
            for gap_coord in gaps:
                gap_st = gap_coord[0]
                gap_en = gap_coord[1]
                adjusted_gap_st = gap_st - gap_accu_len
                gap_positions.append(adjusted_gap_st)
                gap_accu_len+=gap_en-gap_st
            return(gap_positions)

        def cal_gap_win(seq):
            gap = re.compile("-+")
            #find gaps
            gaps=[]
            for m in gap.finditer(seq):
                gaps.append(m.span())
            return(gaps)

        def cal_leading_gaps(seq):
            gap = re.compile("^-+")
            #find gaps
            gaps=[]
            for m in gap.finditer(seq):
                gaps.append(m.span())
            return(gaps)

        def cal_trailing_gaps(seq):
            gap = re.compile("-+$")
            #find gaps
            gaps=[]
            for m in gap.finditer(seq):
                gaps.append(m.span())
            return(gaps)

        def get_amp_abs_coords(seq, coords_list): # get the coords in reference to the amp, not the revcom in cases of '-'
            new_coord_list = []
            for coords in coords_list:
                if coords[2]=="+":
                    new_coord_list.append(coords)
                else:
                    st = coords[0]
                    en = coords[1]
                    new_st = len(seq) - en
                    new_en = len(seq) - st
                    new_coord_list.append((new_st, new_en, coords[2], coords[3]))
            return(new_coord_list)


        def check_deletion_in_cds(read,amp_cds_coords):
            '''
            input: read, amp_cds_coords
            this function will first map the gap in the reads to the amplicon
            then check if any of the gap(deletions) are in the cds
            returns a flag indicating if deletions are found in the cds
            '''
            read_gap_loc = cal_gap_positions(read)
            #check gap in the read (deletion in the reads)
            deletion_in_cds_flag = False
            for coord_set in amp_cds_coords:
                for read_gap in read_gap_loc:
                    if read_gap>=coord_set[0] and read_gap<=coord_set[1]:
                        deletion_in_cds_flag = True
            return({"deletion_in_cds_flag":deletion_in_cds_flag,
                    "read_gap_loc":read_gap_loc})


        def check_insertion_in_cds(ref,amp_cds_coords):
            '''
            input: ref, amp_cds_coords
            this function will first map the gap(insertions) in the ref to the amplicon
            then check if any of the gap(insertions) are in the cds
            returns a flag indicating if insertions are found in the cds
            '''
            #input: ref, wt_amp_cds_coords
            ref_gap_loc = cal_gap_positions(ref)
            #check gap in the ref (insertion in the reads)
            insertion_in_cds_flag = False
            for coord_set in amp_cds_coords:
                for ref_gap in ref_gap_loc:
                    if ref_gap>=coord_set[0] and ref_gap<=coord_set[1]:
                        insertion_in_cds_flag = True 
            return({"insertion_in_cds_flag":insertion_in_cds_flag,
                    "ref_gap_loc":ref_gap_loc,})

        def find_next_num_div3(num):
            if num%3==0:
                return num
            elif (num+1)%3==0:
                return num+1
            elif (num+2)%3==0:
                return num+2
            
        def find_payload_aa_pos(cds_st, cds_en, payload_st, payload_en, strand, frame):
            frame_adjust = 0
            if frame == 2:
                frame_adjust = 1
            if frame == 3:
                frame_adjust = 2

            if strand == "-":
                en2en = cds_en - payload_en
                #print(en2en)
                aa_st = find_next_num_div3(en2en - frame_adjust)/3
                en2st = cds_en - payload_st
                #print(en2st)
                aa_en = find_next_num_div3(en2st- frame_adjust)/3
                return((int(aa_st))-1, int(aa_en))
            if strand == "+":
                st2st = payload_st - cds_st 
                #print(st2st)
                aa_st = find_next_num_div3(st2st - frame_adjust)/3
                en2st = payload_en - cds_st
                #print(en2st)
                aa_en = find_next_num_div3(en2st- frame_adjust)/3
                return((int(aa_st))-1, int(aa_en))
    
        #fetch ensembl transcript
        def fetch_ensembl_transcript(
            ensembl_transcript_id, 
            exon_annot = False):

            '''
            Fetch the requested Ensembl transcript.
            Get the requested Ensembl transcript, together with exon and
            coding region (CDS) boundaries.
            Parameters
            ----------
            ensembl_transcript_id : str
            the ensembl transcript id, of the form ENST...
            
            Returns
            -------
            `Bio.SeqRecord`
            The requested transcript sequence, in 5' -> 3' order, together
            with exon and CDS features. The coordinates of exons and CDS
            features are relative to the sequence fragment.
            '''
            
            base_url = "http://rest.ensembl.org"

            # First, fetch the transcript sequence
            url = base_url + f"/sequence/id/{ensembl_transcript_id}"

            log.info(f"Querying Ensembl for sequence of {ensembl_transcript_id}")
            response = requests.get(url, { "type": "genomic",
                                        "content-type": "application/json" })

            try:
                response.raise_for_status()
            except requests.HTTPError:
                log.error("Ensembl sequence REST query returned error "
                        "{}".format(response.text))
                raise ValueError(reponse.text)

            response_data = response.json()
            
            try:
                description = response_data['desc'].split(':')
                species = description[1]
                try:
                    chromosome_number = int(description[2])
                except:
                    chromosome_number = str(description[2])
                sequence_left = int(description[3])
                sequence_right = int(description[4])
                transcript_strand = int(description[5])

                if sequence_left > sequence_right:
                    raise ValueError(f"Expected left sequence boundary {sequence_left} "
                                    f"<= right sequence boundary {sequence_right}: did "
                                    "the format of the Ensembl REST response change?")
                
                sequence_id = response_data['id']
                
                seq_str = response_data['seq']

                log.info(f"Retrieved sequence {response_data['desc']} of length "
                        f"{sequence_right - sequence_left} for species {species} on "
                        f"strand {transcript_strand}")
            except (KeyError, ValueError) as e:
                log.error(e)
                log.error('Error parsing sequence metadata from Ensembl REST response - '
                        'did the format of the response change?')
                raise ValueError(e)
            
            seq = Seq(seq_str)
            
            record = SeqRecord(seq, id=sequence_id,
                            description=":".join(description))
            if exon_annot:

                url = base_url + f"/overlap/id/{ensembl_transcript_id}"

                log.info(f"Querying Ensembl for overlaps of {ensembl_transcript_id}")
                response = requests.get(url, { "feature": ["cds", "exon"],
                                            "content-type": "application/json" })
                try:
                    response.raise_for_status()
                except requests.HTTPError:
                    log.error("Ensembl sequence REST query returned error "
                            "{}".format(response.text))
                    raise ValueError(reponse.text)

                response_data = response.json()

                try:
                    # Handle the unlikely event of a single piece of information overlapping a lonely transcript

                    if not hasattr(response_data, '__iter__'):
                        response_data = [response_data]

                    for response_datum in response_data:
                        if response_datum['Parent'] != ensembl_transcript_id:
                            continue

                        if response_datum['assembly_name'] != species:
                            continue

                        # Store feature locations 0-indexed from the left-most sequence boundary

                        record.features.append(SeqFeature(
                            location=FeatureLocation(
                                int(response_datum['start']) - sequence_left,
                                int(response_datum['end']) - sequence_left + 1,
                                strand=int(response_datum['strand'])),
                            type=response_datum['feature_type']))
                    num_exon_boundaries = len([f for f in record.features
                                            if f.type == 'exon'])

                    num_cds_boundaries = len([f for f in record.features
                                            if f.type == 'cds'])

                    log.info(f"Retrieved {num_exon_boundaries} exons and "
                            f"{num_cds_boundaries} coding regions for transcript "
                            f"{ensembl_transcript_id}")
                except (KeyError, ValueError) as e:
                    log.error(e)
                    log.error('Error parsing overlap metadata from Ensembl REST response - '
                            'did the format of the response change?')
                    raise ValueError(e)

            record.annotations['reference_species'] = species
            record.annotations['reference_chromosome_number'] = chromosome_number
            record.annotations['reference_left_index'] = sequence_left
            record.annotations['reference_right_index'] = sequence_right
            record.annotations['transcript_strand'] = transcript_strand

            # Finally, sort features by their start locations
            record.features.sort(key=lambda f: f.location.start)
                
            return record


        ##################################################
        #prepration work before processing the alignments#
        ##################################################

        #fetch transcript
        mytranscript = fetch_ensembl_transcript(ensembl_transcript_id = ENST_ID, exon_annot = True)

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
        #alignment
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
                writehandle.write(f"Aligned_Sequence	Reference_Sequence\tReference_Name\tRead_Status\tn_deleted\tn_inserted\tn_mutated\t#Reads\t%Reads\tGenotype\n")
                for line in filehandle:
                    line_deco = line.decode()
                    writehandle.write(line_deco.rstrip())

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
                    #############################
                    #Start calculating genotypes#
                    #############################            
                    #############################
                    #read mapped to HDR amplicon#
                    #############################
                    if Reference_Name == "HDR": # 
                        if read == ref:
                            #print("\tHDR allele (perfect HDR edit)", end="\n")
                            writehandle.write("\tperfect HDR edit\n")
                        ##########################
                        #mutations only, no indel#
                        ##########################
                        elif int(n_deleted)==0 and int(n_inserted)==0: 
                            #print(f"\n{read}\n{ref}\t{n_deleted}\t{n_inserted}\t{n_mutated}\t{n_Reads}\t{perc_Reads}", end='\n')
                            payload_genotype = "missing or nc"
                            protein_correct = False

                            #check wt protein and payload
                            counter=0
                            translations=[]
                            for idx, coord_set in enumerate(HDR_amp_cds_coords):                      
                                translation = translate_amp(seq=read, start = coord_set[0], end = coord_set[1], strand = coord_set[2], frame = coord_set[3])
                                #print("\n")
                                #print(wt_translation[idx].strip('*'))
                                
                                coord_set2 = HDR_amp_cds_coords_PosRef[idx]
                                if payload_coord[0]>=coord_set2[0] and payload_coord[1]<=coord_set2[1]: # payload is in this cds, check payload translation
                                    this_payload_aa_pos = find_payload_aa_pos(cds_st=coord_set2[0],cds_en=coord_set2[1],payload_st=payload_coord_in_HDR_amp[0],payload_en=payload_coord_in_HDR_amp[1],strand = coord_set2[2],frame = coord_set2[3])  
                                    #print(this_payload_aa_pos)
                                    payloadtranslation =  translation[this_payload_aa_pos[0]: this_payload_aa_pos[1]]
                                    #print(payloadtranslation)
                                    if payloadtranslation == payload_wt_translation:
                                        payload_genotype = "wt"
                                    elif payloadtranslation == payload_HDR_translation:
                                        payload_genotype = "HDR"
                                    else:
                                        payload_genotype = "mutant"
                                        
                                    #remove payload and check translation
                                    this_translation_noPL = translation[0:this_payload_aa_pos[0]] + translation[this_payload_aa_pos[1]:]
                                    translation_noPL = wt_translation[idx][0:this_payload_aa_pos[0]] + wt_translation[idx][this_payload_aa_pos[1]:]
                                    #print(f"this trans (no PL): {this_translation_noPL} wt trans: {translation_noPL}")
                                    if this_translation_noPL == translation_noPL: # mark translatin correct
                                        counter+=1
                                else: #no payload, check the direct translation
                                    if HDR_translation[idx] == translation: # mark translatin correct
                                        counter+=1

                            if counter==len(HDR_amp_cds_coords):
                                protein_correct=True

                            #produce output
                            if protein_correct==True:
                                #print(f"\twt protein + {payload_genotype} payload)")
                                writehandle.write(f"\twt protein + {payload_genotype} SNP\n")
                            elif protein_correct==False:
                                #print(f"wt protein + {payload_genotype} payload)")
                                writehandle.write(f"\tmutant protein + {payload_genotype} SNP\n")    

                        ########################
                        #indels in HDR allele  #
                        ########################
                        elif (int(n_deleted)!=0 and int(n_inserted)==0):  #deletion in the read，no insertion
                            #check if deletions are in wt cds
                            deletion_in_HDR_amp_cds_noPL_Flag = check_deletion_in_cds(read = read, amp_cds_coords = HDR_amp_cds_coords_PosRef)["deletion_in_cds_flag"]
                            #check if deletions are in payload cds
                            deletion_in_payload_Flag = check_deletion_in_cds(read = read, amp_cds_coords = [payload_coord])["deletion_in_cds_flag"]
                            #produce output
                            if (deletion_in_HDR_amp_cds_noPL_Flag==True and deletion_in_payload_Flag==True):
                                #print("\tHDR allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write("\tmutant protein + mutant SNP\n")
                            elif (deletion_in_HDR_amp_cds_noPL_Flag==True and deletion_in_payload_Flag==False):
                                #print("\tHDR allele (mutant protein + HDR payload)", end="\n")
                                writehandle.write("\tmutant protein + HDR SNP\n")
                            elif (deletion_in_HDR_amp_cds_noPL_Flag==False and deletion_in_payload_Flag==True):
                                #print("\tHDR allele (wt protein + mutant payload)", end="\n")
                                writehandle.write("\twt protein + mutant SNP\n")    
                            elif (deletion_in_HDR_amp_cds_noPL_Flag==False and deletion_in_payload_Flag==False):
                                #print("\tHDR allele (wt protein + HDR payload)", end="\n")   
                                writehandle.write("\twt protein + HDR SNP\n")      
                        elif (int(n_deleted)==0 and int(n_inserted)!=0):  #insertion in the read, no deletion
                            #check if insertion are in wt cds
                            insertion_in_HDR_amp_cds_noPL_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = HDR_amp_cds_coords_PosRef)["insertion_in_cds_flag"]
                            #check if insertion are in payload cds
                            insertion_in_payload_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = [payload_coord])["insertion_in_cds_flag"]
                            #produce output
                            if (insertion_in_HDR_amp_cds_noPL_Flag==True and insertion_in_payload_Flag==True):
                                #print("\tHDR allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write("\tmutant protein + mutant SNP\n")
                            elif (insertion_in_HDR_amp_cds_noPL_Flag==True and insertion_in_payload_Flag==False):
                                #print("\tHDR allele (mutant protein + HDR payload)", end="\n")
                                writehandle.write("\tmutant protein + HDR SNP\n")
                            elif (insertion_in_HDR_amp_cds_noPL_Flag==False and insertion_in_payload_Flag==True):
                                #print("\tHDR allele (wt protein + mutant payload)", end="\n")
                                writehandle.write("\twt protein + mutant SNP\n")
                            elif (insertion_in_HDR_amp_cds_noPL_Flag==False and insertion_in_payload_Flag==False):
                                #print("\tHDR allele (wt protein + HDR payload)", end="\n") 
                                writehandle.write("\twt protein + HDR SNP\n") 
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
                            deletion_in_HDR_amp_cds_noPL_Flag = check_deletion_in_cds(read = read_trimmed, amp_cds_coords = HDR_amp_cds_coords_PosRef)["deletion_in_cds_flag"]
                            #check if deletions are in payload cds
                            deletion_in_payload_Flag = check_deletion_in_cds(read = read_trimmed, amp_cds_coords = [payload_coord])["deletion_in_cds_flag"]

                            #summarize flag
                            indel_in_HDR_amp_cds_noPL_Flag = any([insertion_in_HDR_amp_cds_noPL_Flag,
                                                                deletion_in_HDR_amp_cds_noPL_Flag])
                            indel_in_payload_Flag = any([insertion_in_payload_Flag,
                                                        deletion_in_payload_Flag])     
                            #produce output
                            if (indel_in_HDR_amp_cds_noPL_Flag==True and indel_in_payload_Flag==True):
                                #print("\tHDR allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write("\tmutant protein + mutant SNP\n")
                            elif (indel_in_HDR_amp_cds_noPL_Flag==True and indel_in_payload_Flag==False):
                                #print("\tHDR allele (mutant protein + HDR payload)", end="\n")
                                writehandle.write("\tmutant protein + HDR SNP\n")    
                            elif (indel_in_HDR_amp_cds_noPL_Flag==False and indel_in_payload_Flag==True):
                                #print("\tHDR allele (wt protein + mutant payload)", end="\n")
                                writehandle.write("\twt protein + mutant SNP\n")    
                            elif (indel_in_HDR_amp_cds_noPL_Flag==False and indel_in_payload_Flag==False):
                                #print("\tHDR allele (wt protein + HDR payload)", end="\n")
                                writehandle.write("\twt protein + HDR SNP\n")    
                    #############################
                    #read mapped to wt amplicon #
                    #############################        
                    else: 
                        if read == ref:
                            #print("\twt allele", end="\n")
                            writehandle.write("\twt allele\n")
                        ##########################
                        #mutations only, no indel#
                        ##########################
                        elif int(n_deleted)==0 and int(n_inserted)==0: 
                            #print(f"\n{read}\n{ref}\t{n_deleted}\t{n_inserted}\t{n_mutated}\t{n_Reads}\t{perc_Reads}", end='\n')
                            payload_genotype = "missing or nc"
                            protein_correct = False

                            #check wt protein and payload
                            counter=0
                            translations=[]
                            for idx, coord_set in enumerate(wt_amp_cds_coords):                      
                                translation = translate_amp(seq=read, start = coord_set[0], end = coord_set[1], strand = coord_set[2], frame = coord_set[3])
                                #print("\n")
                                #print(wt_translation[idx].strip('*'))
                                
                                coord_set2 = HDR_amp_cds_coords_PosRef[idx]
                                if payload_coord[0]>=coord_set2[0] and payload_coord[1]<=coord_set2[1]: # payload is in this cds, check payload translation
                                    this_payload_aa_pos = find_payload_aa_pos(cds_st=coord_set2[0],cds_en=coord_set2[1],payload_st=payload_coord_in_HDR_amp[0],payload_en=payload_coord_in_HDR_amp[1],strand = coord_set2[2],frame = coord_set2[3])  
                                    #print(this_payload_aa_pos)
                                    payloadtranslation =  translation[this_payload_aa_pos[0]: this_payload_aa_pos[1]]
                                    #print(payloadtranslation)
                                    if payloadtranslation == payload_wt_translation:
                                        payload_genotype = "wt"
                                    elif payloadtranslation == payload_HDR_translation:
                                        payload_genotype = "HDR"
                                    else:
                                        payload_genotype = "mutant"
                                        
                                    #remove payload and check translation
                                    this_translation_noPL = translation[0:this_payload_aa_pos[0]] + translation[this_payload_aa_pos[1]:]
                                    translation_noPL = wt_translation[idx][0:this_payload_aa_pos[0]] + wt_translation[idx][this_payload_aa_pos[1]:]
                                    #print(f"this trans (no PL): {this_translation_noPL} wt trans: {translation_noPL}")
                                    if this_translation_noPL == translation_noPL: # mark translatin correct
                                        counter+=1
                                else: #no payload, check the direct translation
                                    if HDR_translation[idx] == translation: # mark translatin correct
                                        counter+=1

                            if counter==len(HDR_amp_cds_coords):
                                protein_correct=True

                            #produce output
                            if protein_correct==True:
                                #print(f"\twt protein + {payload_genotype} payload)")
                                writehandle.write(f"\twt protein + {payload_genotype} SNP\n")
                            elif protein_correct==False:
                                #print(f"\tmutant protein + {payload_genotype} payload)")
                                writehandle.write(f"\tmutant protein + {payload_genotype} SNP\n")
                        ########################
                        #indels in wt allele  #
                        ########################
                        elif (int(n_deleted)!=0 and int(n_inserted)==0):  #deletion in the read，no insertion
                            #check if deletions are in wt cds
                            deletion_in_wt_amp_cds_noPL_Flag = check_deletion_in_cds(read = read, amp_cds_coords = wt_amp_cds_coords_PosRef)["deletion_in_cds_flag"]
                            #check if deletions are in payload cds
                            deletion_in_payload_Flag = check_deletion_in_cds(read = read, amp_cds_coords = [payload_coord])["deletion_in_cds_flag"]
                            #produce output
                            if (deletion_in_wt_amp_cds_noPL_Flag==True and deletion_in_payload_Flag==True):
                                #print("\twt allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write("\tmutant protein + mutant SNP\n")
                            elif (deletion_in_wt_amp_cds_noPL_Flag==True and deletion_in_payload_Flag==False):
                                #print("\twt allele (mutant protein + HDR payload)", end="\n")
                                writehandle.write("\tmutant protein + HDR SNP\n")
                            elif (deletion_in_wt_amp_cds_noPL_Flag==False and deletion_in_payload_Flag==True):
                                #print("\twt allele (wt protein + mutant payload)", end="\n")
                                writehandle.write("\twt protein + mutant SNP\n")    
                            elif (deletion_in_wt_amp_cds_noPL_Flag==False and deletion_in_payload_Flag==False):
                                #print("\twt allele (wt protein + HDR payload)", end="\n")   
                                writehandle.write("\twt protein + HDR SNP\n")      
                        elif (int(n_deleted)==0 and int(n_inserted)!=0):  #insertion in the read, no deletion
                            #check if insertion are in wt cds
                            insertion_in_wt_amp_cds_noPL_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = wt_amp_cds_coords_PosRef)["insertion_in_cds_flag"]
                            #check if insertion are in payload cds
                            insertion_in_payload_Flag = check_insertion_in_cds(ref = ref, amp_cds_coords = [payload_coord])["insertion_in_cds_flag"]
                            #produce output
                            if (insertion_in_wt_amp_cds_noPL_Flag==True and insertion_in_payload_Flag==True):
                                #print("\twt allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write("\tmutant protein + mutant SNP\n")
                            elif (insertion_in_wt_amp_cds_noPL_Flag==True and insertion_in_payload_Flag==False):
                                #print("\twt allele (mutant protein + HDR payload)", end="\n")
                                writehandle.write("\tmutant protein + HDR SNP\n")
                            elif (insertion_in_wt_amp_cds_noPL_Flag==False and insertion_in_payload_Flag==True):
                                #print("\twt allele (wt protein + mutant payload)", end="\n")
                                writehandle.write("\twt protein + mutant SNP\n")
                            elif (insertion_in_wt_amp_cds_noPL_Flag==False and insertion_in_payload_Flag==False):
                                #print("\twt allele (wt protein + HDR payload)", end="\n") 
                                writehandle.write("\twt protein + HDR SNP\n") 
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
                            deletion_in_wt_amp_cds_noPL_Flag = check_deletion_in_cds(read = read_trimmed, amp_cds_coords = wt_amp_cds_coords_PosRef)["deletion_in_cds_flag"]
                            #check if deletions are in payload cds
                            deletion_in_payload_Flag = check_deletion_in_cds(read = read_trimmed, amp_cds_coords = [payload_coord])["deletion_in_cds_flag"]

                            #summarize flag
                            indel_in_wt_amp_cds_noPL_Flag = any([insertion_in_wt_amp_cds_noPL_Flag,
                                                                deletion_in_wt_amp_cds_noPL_Flag])
                            indel_in_payload_Flag = any([insertion_in_payload_Flag,
                                                        deletion_in_payload_Flag])     
                            #produce output
                            if (indel_in_wt_amp_cds_noPL_Flag==True and indel_in_payload_Flag==True):
                                #print("\twt allele (mutant protein + mutant payload)", end="\n")
                                writehandle.write("\tmutant protein + mutant SNP\n")
                            elif (indel_in_wt_amp_cds_noPL_Flag==True and indel_in_payload_Flag==False):
                                #print("\twt allele (mutant protein + HDR payload)", end="\n")
                                writehandle.write("\tmutant protein + HDR SNP\n")    
                            elif (indel_in_wt_amp_cds_noPL_Flag==False and indel_in_payload_Flag==True):
                                #print("\twt allele (wt protein + mutant payload)", end="\n")
                                writehandle.write("\twt protein + mutant SNP\n")    
                            elif (indel_in_wt_amp_cds_noPL_Flag==False and indel_in_payload_Flag==False):
                                #print("\twt allele (wt protein + HDR payload)", end="\n")
                                writehandle.write("\twt protein + HDR SNP\n") 

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

                    #print(f"{Reference_Name}\t{n_Reads}\t{perc_Reads}\t{n_deleted}\t{n_inserted}\t{n_mutated}\t{genotype}", end='\n')
                    if genotype == "perfect HDR edit":
                        genotype_freq["HDR_perfect"] += float(perc_Reads)
                    elif genotype == "wt protein + HDR SNP":
                        genotype_freq["wtProt_hdrSNP"] += float(perc_Reads)
                    elif genotype == "wt protein + mutant SNP":
                        genotype_freq["wtProt_mutSNP"] += float(perc_Reads)  
                    elif genotype == "wt protein + wt SNP":
                        genotype_freq["wtProt_wtSNP"] += float(perc_Reads)                  
                    elif genotype == "mutant protein + HDR SNP":
                        genotype_freq["mutProt_hdrSNP"] += float(perc_Reads)                  
                    elif genotype == "mutant protein + mutant SNP":
                        genotype_freq["mutProt_mutSNP"] += float(perc_Reads)                  
                    elif genotype == "mutant protein + wt SNP":
                        genotype_freq["mutProt_wtSNP"] += float(perc_Reads)      
                    elif genotype == "wt allele":
                        genotype_freq["wt_allele"] += float(perc_Reads)  


        #write genotype to file
        with open(os.path.join(zip_dir,'genotype_frequency.csv'),'w') as writehandle:
            writehandle.write(f"wt_allele,"
                            f"HDR_perfect,"
                            f"wtProt_hdrSNP,"
                            f"mutProt_hdrSNP,"
                            f"mutProt_mutSNP,"
                            f"wtProt_mutSNP,"
                            f"mutProt_wtSNP,"
                            f"wtProt_wtSNP\n")
            writehandle.write(f"{float(genotype_freq['wt_allele']):.4f}%,"
                            f"{float(genotype_freq['HDR_perfect']):.4f}%,"
                            f"{float(genotype_freq['wtProt_hdrSNP']):.4f}%,"
                            f"{float(genotype_freq['mutProt_hdrSNP']):.4f}%,"
                            f"{float(genotype_freq['mutProt_mutSNP']):.4f}%,"
                            f"{float(genotype_freq['wtProt_mutSNP']):.4f}%,"
                            f"{float(genotype_freq['mutProt_wtSNP']):.4f}%,"
                            f"{float(genotype_freq['wtProt_wtSNP']):.4f}%\n")


    except Exception  as e:
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
