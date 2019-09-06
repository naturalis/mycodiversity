"""**********************************************************"""
"""
SCRIPT 2: This is the post-processing script of the PROFUNGIS pipeline
It allows to update a reference ZOTU file with a new FASTA file
It also updates two mappings: the tracker and the ZOTU FASTA file mapping 
Please read the README file for details
"""

import os
import csv
import pandas
import re
from Bio import SeqIO
from random import randint


"""
check last existing PK
"""
def fetchdigital(magicnumber):
        count_digits = 0
        digit_value = 10
        while(magicnumber>0):
            magicnumber = magicnumber/digit_value
            count_digits = count_digits + 1
        return count_digits

"""
Retrieve FASTA file
"""
def fetchsrrname(fastafile):
    srr_filename = os.path.basename(fastafile)
    srr_label = os.path.splitext(srr_filename)[0]
    srr_output_parse = srr_label.split("_")
    srr_idname = srr_output_parse[0]
    return srr_idname

"""
Generation of outputs from updated dataframes
"""
def createdf_nohead(csv_df):
    generate_df_nohead = pandas.read_csv(csv_df, header=0)
    return generate_df_nohead
    
def createdf_head(csv_df):
    generate_df_head = pandas.read_csv(csv_df)
    return generate_df_head
    
def create_csv_from_mapping(srrfilename, nameheaders, fastamapping):
    with open(srrfilename, "w") as fp:
        writer = csv.writer(fp)
        writer.writerow(nameheaders)
        writer.writerows(fastamapping)
    fp.close()
    
def drop_column_fromdf(df_name, column_name):
    drop_col = df_name.drop(column_name, axis=1)
    return drop_col

"""
Generation of new PK, keep track of last PK of ZOTU reference file
"""
def generate_newpk_refset(ref_keys, update_keys):
    pk_suffix = "MDDBOTU" 
    new_pk_values = []
    new_pk_keys_list = []
    increment_pk = 1 
    max_value_pk =  max(ref_keys) #list_store_keys, #the next new pk will be incremented based on the max pk of refset
    pos_suffix_pk = max_value_pk.find(pk_suffix) + 7
    id_suffix_pk = max_value_pk[:pos_suffix_pk]
    max_value_key = max_value_pk[pos_suffix_pk:]
    max_number_pk = int(max_value_key)
    number_keys_generate = len(update_keys) #list_update_keys, #how many new keys do you need? save this var on file
    while(increment_pk < len(update_keys) + 1):
        new_pk_values.append(max_number_pk + increment_pk)
        increment_pk += 1
    max_count_keyvalue = max_number_pk + number_keys_generate
    for newseq_record in range(max_number_pk + 1, max_number_pk + number_keys_generate + 1): #for appending 0s, the amount depends on the digit, lim 10^6
        generate_new_primary_key = "0" * (6 - len(str(newseq_record))) + str(newseq_record)
        new_pk_keys_list.append(pk_suffix + str(generate_new_primary_key))
    new_key_zotu_mapping = dict(zip(update_keys, new_pk_keys_list))
    return new_key_zotu_mapping
    
def save_df_file(df_name, file_name):
    encoding_form = "utf-8"
    index_no = False
    df_name.to_csv(file_name, encoding = encoding_form, index = index_no) 
    
def change_column_order(dfname, colname, indexpos):
    cols = dfname.columns.tolist()
    cols.remove(colname)
    cols.insert(indexpos, colname)
    return dfname[cols]
"""
save tracker file
#keep track of how many new seq are in the fasta file, and how many already exist
"""    
def save_record_track(trackfile, tracker): 
    with open(trackfile, 'a') as f:
        tracker.to_csv(f, header = False, index = False)
    f.close()
    
    
    """
    provide and generate output file names
    """   
def Main(fastafiletest, referencefile, mappingfile):
    otu_id, otu_sequence = [], []
    list_updated_keys, list_update_keys, list_store_keys = [], [], [] #holders for keeping keys,ids
    header_columns_fasta_def = ["refsequence_pk","sequence"]    #note in what you are fetching
    generate_new_refseq_list = "appendnewseqinref1.csv"         #name for new refset generated, thus we keep also previous version
    marker_label = "marker_type"
    length_label = "refseq_length"
    file_extension = ".csv"
    oldkey_suffix = "Zotu"
    marker_type = "ITS2" 
    updated_refseq_table = "refseq_table_up2.csv"
    updated_mapping_pk_srr_zotu = "updated_mapping_pk_srr_zotu.csv"
    record_track = "record_track.csv"
    
    oldkey_suffix_length = len(oldkey_suffix)
    
    for record in SeqIO.parse(fastafiletest, 'fasta'):  #create df from the new fasta list
        otu_id.append(record.id)
        otu_sequence.append(record.seq)
        
    srr_name = fetchsrrname(fastafiletest)
    new_sequence_table = srr_name + file_extension
    record_list = [list(item) for item in list(zip(otu_id, otu_sequence))]
    create_csv_from_mapping(new_sequence_table, header_columns_fasta_def, record_list)
    sequences_to_upload = createdf_nohead(new_sequence_table)
    reference_sequences = createdf_nohead(referencefile)
    count_refseq = len(reference_sequences.index)
    print 'THE LENGTH OF REFERENCE FILE:'
    print count_refseq #keep header here and drop what's not necessary
    seq_in_srr = len(otu_sequence)
    print 'THE LENGTH OF FASTA FILE:'
    print seq_in_srr
    """
    launch the functions needed to update a reference file and reference mapping
    """
    reference_srr_mapping = pandas.read_csv(mappingfile)
    newversion_ref_seq_no_marker = drop_column_fromdf(reference_sequences, marker_label)
    newversion_ref_seq = drop_column_fromdf(newversion_ref_seq_no_marker, length_label)  #newversion_ref_seq df with pk, seq only
    append_new_zotus_onrefset = newversion_ref_seq.merge(sequences_to_upload, on=list(newversion_ref_seq), how = "outer") #appending new zotulist
    append_new_zotus_onrefset.columns = ["refsequence_pk", "sequence"] #new header for merging set
    append_new_zotus_onrefset.drop_duplicates(subset=["sequence"], inplace=True, keep="first") #keep first duplicates only! This means that the original PK is kept
    
    #keep mapping of new zotus to which pk of ref set they map to (need to create a new df for this because the previous df dropped duplicates
    existing_zotus_onrefset = newversion_ref_seq.merge(sequences_to_upload, on = list(newversion_ref_seq), how = 'outer')
    fetch_existing_zotus_onrefset = existing_zotus_onrefset[existing_zotus_onrefset.duplicated(subset = ["sequence"], keep = 'first')][["refsequence_pk","sequence"]] #keep pk of existing sequences
    fetch_new_zotus_onrefset = existing_zotus_onrefset[existing_zotus_onrefset.duplicated(subset = ["sequence"], keep = 'last')][["refsequence_pk","sequence"]] #keep zotus of existing sequences
    #create the mapping, while keeping matching sequences
    remove_refseq_duplicate = fetch_existing_zotus_onrefset.merge(fetch_new_zotus_onrefset, on="sequence") #df which contains the mapping of existing sequences
    mapping_existing_zotu_seq = remove_refseq_duplicate.drop("sequence", axis = 1) #this contains the mapping new_zotus | refset pk_existing seq
    #GENERATE NEW KEYS, KEEP TRACK HOW MANY KEYS TO GENERATE
    list_update_keysreplace = append_new_zotus_onrefset["refsequence_pk"].values
    for previous_key in list_update_keysreplace:
        if previous_key[0:oldkey_suffix_length] == oldkey_suffix:
            list_update_keys.append(previous_key)
        else:
            list_store_keys.append(previous_key)
            
    updated_key_zotu_mapping = generate_newpk_refset(list_store_keys, list_update_keys)
    #now replace new pks with previous keys (zotus)
    append_new_zotus_onrefset['refsequence_pk'] = append_new_zotus_onrefset['refsequence_pk'].replace(updated_key_zotu_mapping) #new ref set, append attributes
    refseq_length = append_new_zotus_onrefset.sequence.str.len()
    append_new_zotus_onrefset.insert(loc=2, column = "marker_type", value = marker_type) #new refseq detected
    append_new_zotus_onrefset.insert(loc=3, column = "refseq_length", value = refseq_length) #append length of each seq
    save_df_file(append_new_zotus_onrefset, updated_refseq_table) #save   
    print 'new keys mapping:'
    for key, value in updated_key_zotu_mapping.iteritems(): #NOW APPEND MARKER AND LENGTH
        print key, value
    
    srr_previous = reference_srr_mapping["srr_name"].tolist() #1. old existing keys mapping
    mapping_existing_zotu_seq.insert(loc=2, column = 'srr_name', value = srr_name) #list which needs the new SRR value
    updated_key_zotu_mapping = pandas.DataFrame(updated_key_zotu_mapping.items(), columns = ["zotu_id", "refsequence_pk"]) #convert dictionary in dataframe
    updated_key_zotu_mapping.insert(loc=2, column = 'srr_name', value = srr_name)
    mapping_existing_zotu_seq.columns = ['zotu_id','refsequence_pk','srr_name']
    df_mapping_zotu_pk = pandas.concat([updated_key_zotu_mapping, mapping_existing_zotu_seq]) #MERGE PREVIOUS KEYS TO NEW KEYS FOR MAPPING TABLE
    col_to_change = 'zotu_id'
    column_to_change1 = 'refsequence_pk'
    df_mapping_pk_srr_zotu = change_column_order(df_mapping_zotu_pk, col_to_change, 2)
    df_mapping_pk_srr_zotu = df_mapping_pk_srr_zotu.merge(append_new_zotus_onrefset, on = "refsequence_pk", how = "left")
    df_mapping_pk_srr_zotu = drop_column_fromdf(df_mapping_pk_srr_zotu, "refseq_length")
    df_mapping_pk_srr_zotu = drop_column_fromdf(df_mapping_pk_srr_zotu, "marker_type")
    previous_mapping_pk_srr = createdf_head(mappingfile)
    previous_mapping_pk_srr = change_column_order(previous_mapping_pk_srr, column_to_change1, 0)
    prev_mapping_pk_srr = previous_mapping_pk_srr.rename(columns = {"otuseq_id":"zotu_id"})
    df_new_mapping = prev_mapping_pk_srr.append(df_mapping_pk_srr_zotu) #NOW MERGE THE PREVIOUS MAPPING WITH NEW MAPPING
    save_df_file(df_new_mapping, updated_mapping_pk_srr_zotu)
    updatecol_headers = ["update_srr", "seq_fasta", "previous_seq", "update_seq"]
    amount_seq_first = prev_mapping_pk_srr.shape[0]
    amount_seq_update = append_new_zotus_onrefset.shape[0]
    first_srr_name = prev_mapping_pk_srr['srr_name']
    first_srr = [[first_srr_name[0], amount_seq_first, 0, amount_seq_first]] 
    for i in first_srr:
        print i
    next_srr = [[srr_name, seq_in_srr, amount_seq_first, amount_seq_update]]
    update_record = pandas.DataFrame(next_srr) #UPDATE_RECORD IS NOW SAVED ON A FILE   
    save_record_track(record_track,update_record)
            
#Ask the user the fasta file ref
if __name__ == '__main__':	
	import sys
	try:
		Fastafile = sys.argv[1]	
		ReferenceSet = sys.argv[2]
		MappingSet = sys.argv[3]
	except:
        #give your csv table contains the fasta format of all sequences
		#print "give file name, needs to be in a fasta format"
		#print "For example:'srr0000.fa' such as PROFUNGIS output" 
		#print len(sys.argv)
		print "Usage: CreateCsvFromFasta.py <filename>.fa <refset>.csv <mapping>.csv"
		print "Missing files, please provide all 3 files: fasta file to upload, a reference table and a mapping table"
		raise SystemExit
	else:
	    if len(sys.argv) == 4:
	        Main(Fastafile, ReferenceSet, MappingSet)
	    else:   
	        print "Missing files, please provide all 3 files"
	        raise SystemExit