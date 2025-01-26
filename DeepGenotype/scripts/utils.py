#load packages
import sys, linecache, re, requests, logging
from Bio.pairwise2 import align
from Bio.pairwise2 import format_alignment
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

log = logging.getLogger("Process_alleles_freq_table.py")
log.setLevel(logging.ERROR) #set the level of warning displayed

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))

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
        coding_dna=coding_dna.replace("-", "N")
        trans = str(Seq(coding_dna).translate())    #translate
        return(trans)      
    else:
        seq = seq[start:end]
        coding_dna = adjust_frame(seq,frame) #adjust frame
        coding_dna = coding_dna.replace("-", "N")
        trans = str(Seq(coding_dna).translate())    #translate
        return(trans)      

def translate_payload(seq, start, end, strand, frame): # for payloads, the start, end refers to the seq (HDR amplicon) , not the revcom 
    #seq = seq[start:end].replace("-", "")
    if strand == "-":
        seq = seq[start:end]
        seq = revcom(seq) #revcom
        coding_dna = adjust_frame(seq,frame) #adjust frame
        coding_dna=coding_dna.replace("-", "N")
        trans = str(Seq(coding_dna).translate())    #translate
        return(trans)      
    else:
        seq = seq[start:end]
        coding_dna = adjust_frame(seq,frame) #adjust frame
        coding_dna = coding_dna.replace("-", "N")
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
            if (read_gap-1)>=coord_set[0] and (read_gap-1)<=coord_set[1]: # read gap is 1-based index, to match 0-based index of coord, use read_gap-1
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
            if (ref_gap-1)>=coord_set[0] and (ref_gap-1)<=coord_set[1]: # ref gap is 1-based index, to match 0-based index of coord, use ref_gap-1
                insertion_in_cds_flag = True 
    return({"insertion_in_cds_flag":insertion_in_cds_flag,
            "ref_gap_loc":ref_gap_loc,})

def check_deletion_overlap_cds(read,amp_cds_coords):
    read_gap_win = cal_gap_win(read)
    #check for overlap between gap and cds
    deletion_overlap_cds_flag = False
    for coord_set in amp_cds_coords:
        for gap_win in read_gap_win:
            if ( gap_win[0] <= coord_set[0] <= gap_win[1]) or (gap_win[0] <= coord_set[1] <= gap_win[1]):
                deletion_overlap_cds_flag = True    
    return({"deletion_overlap_cds_flag":deletion_overlap_cds_flag,
            "read_gap_win":read_gap_win})

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

def compute_percent_identity(ref, num_of_mismatches):
    return 100 - (num_of_mismatches / len(ref)) * 100

def compute_num_of_mismatches(read, ref):
    return sum(1 for r, t in zip(read, ref) if r != t)

def compute_weighted_average(values, weights):
    if len(values) != len(weights):
        raise ValueError("The number of values and weights must be the same.")
    if sum(weights) == 0:
        if len(values) == 0:
            return 0
        else:
            # if the sum of weights is 0, return the average of the values
            return sum(values) / len(values)
    weighted_sum = sum(v * w for v, w in zip(values, weights))
    total_weight = sum(weights)
    return weighted_sum / total_weight

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
    
def check_protein_and_payload(read, wt_amp_cds_coords, HDR_amp_cds_coords_PosRef, payload_wt_translation, payload_HDR_translation, payload_coord, payload_coord_in_HDR_amp, HDR_translation, wt_translation): #must remove insertions in read before calling this function    
    # used in process_alleles_freq_table_SNP.py
    payload_genotype = "missing or nc"
    protein_correct = False

    #check wt protein and payload
    counter=0
    for idx, coord_set in enumerate(wt_amp_cds_coords): # doesn't matter if wt or HDR
        
        #convert gaps into "N"s
        read = re.sub(r'-','N',read)
        
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

    if counter==len(wt_amp_cds_coords):
        protein_correct=True

    return({"payload_genotype":payload_genotype,
            "protein_correct":protein_correct})
            
def check_protein_correctness(read,wt_or_HDR_amp_cds_coords,wt_or_HDR_translation, payload_coord_in_HDR_amp, payload_strand, payload_translation): #must remove insertions in read before calling this function    
    # used in process_alleles_freq_table_INS.py
    payload_correct = False
    protein_correct = False
    #check payload
    translation = translate_payload(seq=read, start =  payload_coord_in_HDR_amp[0], end = payload_coord_in_HDR_amp[1], strand = payload_strand, frame = 1)
    if translation == payload_translation:
        payload_correct = True
    #check wt protein
    counter=0
    for idx, coord_set in enumerate(wt_or_HDR_amp_cds_coords):
        translation = translate_amp(seq=read, start = coord_set[0], end = coord_set[1], strand = coord_set[2], frame = coord_set[3])
        #print("\n")
        #print(wt_translation[idx].strip('*'))
        #print(translation)
        if wt_or_HDR_translation[idx] == translation:
            counter+=1
    if counter==len(wt_or_HDR_amp_cds_coords):
        protein_correct=True
    return({"payload_correct":payload_correct,
            "protein_correct":protein_correct})