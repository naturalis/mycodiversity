# encoding=utf8
import requests
import sys
reload(sys)
sys.setdefaultencoding('utf8')
import os
import re
from bs4 import BeautifulSoup
from Bio import Entrez
import csv
import pandas as pd
import numpy as np
import conf_lit as na
from itertools import islice



__author__ = "airini"

"""
v25k: this version contains: retrieve and fetch data from PubMed, SRA, SRP, SRS, SRX and BioSample
consider to process SRX data -> fetch specific attributes -> exp name, primers and design, retrieve only unique values found in these att
"""


"""
***WHAT the script does***
A DOI is provided in input. Initially the script will try the NCBI DOI PMID converter to retrieve the PMID
If the above will not succeed, the EBI API will be used instead
The secondary source (sequence data repository) will be searched in the PubMed xml record
If the data element above will not be found, it will be requested to the user to provide it manually
Once the sequence data identifier (SRA, SRP, BioProject) is obtained OR provided, 
the SRA object mapping will be retrieved.

The following mappings will be saved:
DOI - PMID ids and metadata associated
PMID - SRA ids 
SRA - SRS - SRP - SRX - SRR - BioSample - BioProject ids
SRR list used as input for PROFUNGIS
SRS metadata will be saved -> only if sample metadata is been provided in SRA. This is given for metadata curation 
SRX metadata will be saved -> used for retrieving attributes used for the input of PROFUNGIS
"""

"""
NOTE
please make sure to add your ncbi user information in the conf_lit file for accessing NCBI data before launching
"""



def addRow(df,ls):
    """
    Adding a row to df - used for updating requested data
    """
    numEl = len(ls)
    newRow = pd.DataFrame(np.array(ls).reshape(1,numEl), columns = list(df.columns))
    df = df.append(newRow, ignore_index=True)
    return df

def retrieveElementValue(ele, xml_find):
	"""
	Searchers for a specific element in a xml
	"""
	item_list = []
	search_element = xml_find.find_all(ele)
	for item in search_element:
		element_value = item.get_text()
		item_list.append(element_value)
	return item_list

def cross_reference_sra(ssv):
	"""
	this checks if the name type falls within a given list of possible accessions
	"""
	list_sra_possible = ["SRA","SRP","SRS","ERP"]
	list_search_possible = []
	list_search_possible.append(ssv)
	save_sra = []
	for i in range(len(list_sra_possible)):
		for ii in range(len(list_search_possible)):
			if re.findall(list_sra_possible[i], list_search_possible[ii]):
				print(str(ii)+': '+list_sra_possible[i])
				save_sra.append(list_search_possible[ii])
	return save_sra
	

def send_request(doi_url, temp_xml):
	"""
	xml dump
	"""
	doi_xml_request = requests.get(doi_url)
	with open(temp_xml, 'wb') as f:
		f.write(doi_xml_request.content)
	f.close()
	return temp_xml

def save_mapping(save_map, output_name, output_path):
	completeName = os.path.join(output_path, output_name)
	save_map.to_csv(completeName, index=False)

def save_df_no_head(save_map, output_name, output_path):
	completeName = os.path.join(output_path, output_name)
	save_map.to_csv(completeName, index=False, header = None)

def handle_article_info(access_account, pubmed_id):
	"""
	handle the PubMed article record, fetch and check title
	"""
	medline_pubmed_db_name = "pubmed"
	Entrez.email = access_account
	handle = Entrez.efetch(medline_pubmed_db_name, id=pubmed_id, retmode = "xml")
	records = Entrez.read(handle) 
	records = records['PubmedArticle']
	for record in records:
		article_pubmed_title = record['MedlineCitation']['Article']['ArticleTitle']
	return article_pubmed_title

def handle_article_data_element(access_account, pubmed_id):
	"""
	E-utilities request (efetch)
	"""
	medline_pubmed_db_name = "pubmed"
	Entrez.email = access_account
	handle = Entrez.efetch(medline_pubmed_db_name, id=pubmed_id, retmode = "xml")
	records = Entrez.read(handle) 
	records = records['PubmedArticle']
	return records

def create_doi_pmid_header():
	doi_pmid_map = "doi_mappings.csv"
	mapping_header = ["DOI","DOI_link","PMID","Pubmed_link", "DOI_article_title"]
	with open(doi_pmid_map, 'wb') as mf:
		wr = csv.writer(mf)
		wr.writerow(mapping_header)
	mf.close()
	return doi_pmid_map

def create_df_doi_pmid(df_map, art_doi, article_link, pmid_id, pmid_link, pmid_name):
	mapping_df = pd.read_csv(df_map)
	mapping_row = []
	mapping_row.extend((art_doi, article_link, pmid_id, pmid_link, pmid_name))
	update_mapping = addRow(mapping_df, mapping_row)
	return update_mapping

def doi_ebi_pmid_converter(a_doi):
	"""
	Request EBI call and fetch items 'pmid' and article 'title'
	"""
	ebi_path = na.ebi_pmid_request["ebi_api_path"]
	ebi_format = na.ebi_pmid_request["ebi_xml_format"]
	ebi_xml_file = "ebi_xml"
	doi_pmid_path_EBI = ebi_path + a_doi + ebi_format
	get_ebi_xml = send_request(doi_pmid_path_EBI, ebi_xml_file)
	contents = open(get_ebi_xml, "r")
	soup = BeautifulSoup(contents, "xml") 
	print "PMID detected..."
	pmid_values = retrieveElementValue('pmid', soup)
	article_titles = retrieveElementValue('title', soup)
	article_title1 = article_titles[0]
	pmid_value1 = pmid_values[0]
	print "Article associated to: " + str(pmid_value1)
	return (pmid_value1, article_title1)
	

def doi_ncbi_pmid_converter_connector(art_doi):
	"""
	Request NCBI call and fetch items 'pmid' and article 'title' (faster)
	"""
	ncbi_account = na.ncbi_access1["account"]
	ncbi_sic_path =  na.ncbi_pmid_request["ncbi_pmid_api_converter"]
	ncbi_access_path = na.ncbi_pmid_request["ncbi_tool_path"]
	doi_link = "http://dx.doi.org/"
	ncbi_pmid_link = na.ncbi_pmid_request["ncbi_pmid_path"]
	temp_pars_idconv = "idconv"	#this contains the api request from utils
	"""create OUTPUT directory"""
	output_dir = "OutputResults"
	try:
		os.mkdir(output_dir)
	except OSError as exc_dir:
		##consider testing: if any(f for f in os.listdir(output_dir) if not f.startswith('.')):
		if any(f for f in os.listdir(output_dir) if not f.startswith('.')):
			print "Warning: path output name <OutputResults> already exists with contents"
			print "Please take a look at the files in this folder, consider saving them somewhere safe as all previous files in this folder will be overwritten."
			print "Check and launch programme again"
			raise SystemExit
	
	article_link = doi_link + art_doi
	"""
	generate variables for requests
	I handle the error
	Here is where I create the first mapping DOI - Pubmed
	"""
	doi_pmid_request_path = ncbi_sic_path + art_doi + ncbi_access_path + ncbi_account
	get_xml = send_request(doi_pmid_request_path, temp_pars_idconv)
	contents = open(get_xml, "r")
	soup = BeautifulSoup(contents, 'xml')

	check_doi = soup.find_all('error')
	if (check_doi):
		print "Not able to recognise DOI: <" + art_doi + ">"
		print "Please check the DOI or replace it with another one"
		raise SystemExit
	else:
		contents = open(get_xml, "r")
		soup = BeautifulSoup(contents, 'xml')
		mappings = soup.find_all('record')
		for mapping in mappings:
				pmid_id = mapping.get('pmid')
				pmc_id = mapping.get('pmcid')
				if (pmid_id) and (pmc_id):
					print 'found mappings..'
					print("PubMed id = %s" % pmid_id)
					print("PMC id = %s" % pmc_id)
					article_name = handle_article_info(ncbi_account, pmid_id)
					print article_name
					doi_pmid_mapping = create_doi_pmid_header()
					ncbi_pmid_record = ncbi_pmid_link + str(pmid_id)
					update_mapping = create_df_doi_pmid(doi_pmid_mapping, art_doi, article_link, pmid_id, ncbi_pmid_record, article_name)
					os.remove(doi_pmid_mapping)
					doi_pmid_mapping2 = "doi_pmid_map.csv"
					save_mapping(update_mapping, doi_pmid_mapping2, output_dir)
				else:
					print "********************************"
					print "*  STEP 1: DOI - PMID mapping  *"
					print "********************************"
					print "\nUnable to find a Pubmed article in association to '" + art_doi + "'"
					print "by using the NCBI converter. EBI API will be tried next"
					pmid_id, article_name = doi_ebi_pmid_converter(art_doi)
					doi_pmid_mapping = create_doi_pmid_header()
					ncbi_pmid_record = ncbi_pmid_link + str(pmid_id)
					checkmap = create_df_doi_pmid(doi_pmid_mapping, art_doi, article_link, pmid_id, ncbi_pmid_record, article_name)
					os.remove(doi_pmid_mapping)
					doi_pmid_mapping2 = "doi_pmid_map.csv"
					save_mapping(checkmap, doi_pmid_mapping2, output_dir)
	return pmid_id

def select_df_columns(data_frame, column_names):
	"""fetch columns from df"""
	new_frame = data_frame.loc[:, column_names]
	return new_frame

def save_sra_objects(sra_meta):
	"""just save everything"""
	path_sra_meta = os.path.dirname(sra_meta)
	df_sra_meta = pd.read_csv(sra_meta, header=0)
	sra_object_list = []
	sra_object_ids = ["Submission","SRAStudy","BioProject","BioSample","ProjectID","Sample","Run","Experiment"]
	run_ids = ["Run"]
	sra_object_df = select_df_columns(df_sra_meta, sra_object_ids)
	run_df = select_df_columns(df_sra_meta, run_ids)
	Profungis_input = "Profungis_accession_list.txt"
	sra_object_mapping = "SRA_object_mapping.csv"
	save_mapping(sra_object_df, sra_object_mapping, path_sra_meta)
	save_df_no_head(run_df, Profungis_input, path_sra_meta)
	#append here path else you can't recognize it
	sra_object_mapping = os.path.join(path_sra_meta, sra_object_mapping)
	return sra_object_mapping

def biosample_retrieve(sra_object_mapping):
	"""
	In order to retrieve sample metadata, I will need the BioSample ids
	"""
	df_object_mapping = pd.read_csv(sra_object_mapping, header = 0)
	#print df_object_mapping
	biosample_id_list = df_object_mapping["BioSample"].tolist()

	#for pos in biosample_id_list:
	#	print pos
	return biosample_id_list

def sra_object_retrieve(sra_id):
	"""
	I 'get' all the SRA info I need
	"""
	#call the following request as did for ncbi_account: ncbi_account = os.getenv("NCBI_ACCOUNT", na.ncbi_access1["account"])
	url_sra_request = na.sra_request["url_sra"]
	url_sra_id_request = url_sra_request + sra_id
	send_sra_request = requests.get(url_sra_id_request)
	sra_content = send_sra_request.content
	sra_object_dump = "SRA_submission_meta.csv"
	output_dir = "OutputResults"
	complete_sra_name = os.path.join(output_dir, sra_object_dump)
	with open(complete_sra_name, "wb") as f:
		f.write(sra_content)
	f.close()
	try:
		request_check = pd.read_csv(complete_sra_name, header = 0)
		first_col_check = request_check["Run"].tolist()
		item_check = first_col_check[0]
		if not is_nan(item_check):
			print "I fetched something.."
		elif is_nan(item_check):
			print "there is something"
	except:
		print "Sorry, SRA not accessible"
	else:
		sra_map_file = save_sra_objects(complete_sra_name)
		biosample_id_list = biosample_retrieve(sra_map_file)
	return biosample_id_list
	

def msg_no_secondary_source(pub):
	print "\n********************************"
	print "*    DATA RETRIEVAL MESSAGE    *"
	print "********************************\n"
	print "****************************************************************"
	print "I am so sorry, unfortunately the Pubmed record: <" + pub + "> "
	print "did not provide the SRA sequence source. "
	print "Would you please provided it for me?"
	print "**HINT**: you must dig in the text of the article, "
	print "In most cases the accession number is found under 'Materials'"
	print "or 'Acknowledgments'"
	print "I am not trained yet to retrieve it from free text\n"
	print "Type <Yes>/<y> if you want to provide the SRA accession or "
	print "<No>/<n> if you would like to come another time"
	print "No worries, all previous mappings will be found in the output folder"
	print "*****************************************************************\n"

def fetch_sra_manually():
	sra_input = raw_input("Please provide SRA submission/study id (ex: SRP043706): ")
	print "Checking for: " + sra_input
	biosamples = []
	try:
		biosamples = sra_object_retrieve(sra_input)
	except:
		print "Sorry I could not recognize <" + sra_input + "> and no data could be found"
		raise SystemExit
	return biosamples

def pmid_sra_connector(pubmed_id):
	"""
	here I try to slurp the pubmed article and see if I can manage to find the sequence id
	"""
	ncbi_account = na.ncbi_access1["account"]
	sample_metadata_list = []
	pubmed_record = handle_article_data_element(ncbi_account, pubmed_id)
	print "********************************"
	print "*  STEP 2: PMID - SRA mapping  *"
	print "********************************"
	for element_value in pubmed_record:
		try:
			databaseinfo = element_value['MedlineCitation']['Article']['DataBankList']
			for d in databaseinfo:
				databank = d["DataBankName"]
				if databank:
					secondarysource = "YES"
					print secondarysource
					accessionrepository = d["AccessionNumberList"]
            		accessionrepositoryname = accessionrepository[0]
            		seq_id = str(accessionrepositoryname)
            		sra_list_check = cross_reference_sra(seq_id) #check the repository value. If SRA type retrieve, else provide name and ids of other rep
            		if not sra_list_check:
            			print "Sequence DataBase detected, but NOT SRA, it is: " + databank 
            			print "This DataBase provided the following accession number(s):" 
            			for record_id in accessionrepository:
            				print record_id
            			print "It is not possible to automatically retrieve anymore"
            			print "********************************"
            			print "*       STEP 2 completed       *"
            			print "********************************"
            			raise SystemExit
            		if sra_list_check:
            				print "SRA detected in PubMed"
            				sample_metadata_list = sra_object_retrieve(seq_id)
            		else:
            			"Could not manage sequence source automatic detection"
		except:
			msg_no_secondary_source(pubmed_id)
			control_input = True
			while control_input:
				sra_check = raw_input("Would you like me to try to continue the mapping and SRA data retrieval? ").lower()
				print sra_check
				if sra_check in ("yes", "y"):
					sample_metadata_list = fetch_sra_manually()
					control_input = False
				elif sra_check in ("no", "n"):
					print "See you next time!"
					raise SystemExit
				else:
					print "Sorry, I did not get your input"
	return sample_metadata_list
			
			
def sra_biosample_connector(biosample_items):
	output_biosample_meta = "SRA_sample_metadata.csv"
	dump_sample_meta = os.path.join("OutputResults", output_biosample_meta)
	output_dir = "OutputResults"
	print "********************************"
	print "* STEP 3: SRA BioSample mapping*"
	print "********************************"
	print "Retrieving sample metadata..."
	size_meta = len(biosample_items)
	print "There are " + str(size_meta) + " Biosamples associated to this study"
	small_test_biosample_items = []
	sampling_test = 10
	sample_iterator = islice(biosample_items, sampling_test)
	#TESTING STEP: use sampling_test, else for all data use the whole list biosample_items
	print "For testing fetch request, I will only fetch and retrieve the first " + str(sampling_test) + " samples associated to this study"
	for sample_item in sample_iterator:
		small_test_biosample_items.append(sample_item)
	for biosample_item in small_test_biosample_items:
		handle = Entrez.efetch(db="biosample", id = biosample_item, format = "xml")
		out_handle = open("OutputResults/sample_meta","wb")
		out_handle.write(handle.read())
		out_handle.close()
		handle.close()
		contents = open("OutputResults/sample_meta","r")
		soup = BeautifulSoup(contents, "xml")
		try:
			fetch_attribute = retrieveElementValue('Attribute',soup)
			with open(dump_sample_meta, 'a') as mf:
				sample_meta_w = csv.writer(mf)
				sample_meta_w.writerow(fetch_attribute)
			mf.close()
		except:
			print "No attributes were found for the samples"
	contents2 = open("OutputResults/sample_meta", "r")
	soup2 = BeautifulSoup(contents2, "xml")
	try:
		list_attributes = soup2.find_all('Attribute')
		attribute_names = []
		for pos in list_attributes:
			get_attribute_name = pos["attribute_name"]
			attribute_names.append(get_attribute_name)
	except:
		print "Sorry, but no sample metadata could be retrieved from the Biosample DB"
	else:
		if attribute_names[0]:
			df_biosample_meta = pd.read_csv(dump_sample_meta, names=attribute_names)
			Biosample_column_values = pd.Series(small_test_biosample_items)
			df_biosample_meta.insert(loc=0, column='Biosample_id', value=Biosample_column_values)
			save_mapping(df_biosample_meta, "Biosample_meta_h.csv", output_dir)
		else:
			print "Sorry but no sample metadata could be saved because there is no metadata attached to this study"


def retrieve_bioproject_title(bioproject_id):
	handle = Entrez.efetch(db="bioproject", id = bioproject_id, format = "xml")
	out_handle = open("bioproject2","wb")
	out_handle.write(handle.read())
	out_handle.close()
	handle.close()
	contents = open("bioproject2", "r")
	soup = BeautifulSoup(contents, "xml") 
	project_name = soup.find('Title')
	project_value = project_name.get_text()
	os.remove("bioproject2")
	return project_value

def is_nan(x):
	return (x is np.nan or x != x)

def retrieve_exp_title(experiment_id):
	handle = Entrez.efetch(db="sra", id = experiment_id, format = "xml")
	out_handle = open("sra_exp1","wb")
	out_handle.write(handle.read())
	out_handle.close()
	handle.close()
	#retrieve element 
	contents = open("sra_exp1", "r")
	soup = BeautifulSoup(contents, "xml") 
	study_title_ele = soup.find('STUDY_TITLE')
	study_title = study_title_ele.get_text()
	return study_title


def update_lipubsra_mapping():
	output_dir = "OutputResults"
	sra_obj_map = output_dir + "/SRA_object_mapping.csv"
	doi_pub_map = output_dir + "/doi_pmid_map.csv"
	try:
		df_sra_obj_map = pd.read_csv(sra_obj_map, header = 0)
		study_list = df_sra_obj_map["SRAStudy"].tolist()
		fetch_study_id = study_list[0]
		bioproject_list = df_sra_obj_map["BioProject"].tolist()
		fetch_bioproject_id = bioproject_list[0]
		experiment_list = df_sra_obj_map["Experiment"].tolist()
		fetch_experiment_id = experiment_list[0]
		if not is_nan(fetch_bioproject_id):
			bioproject_study_title = retrieve_bioproject_title(fetch_bioproject_id)
		elif is_nan(fetch_bioproject_id):
			bioproject_study_title = retrieve_exp_title(fetch_experiment_id)
	except:
		pass
	df_lit_pub_map = pd.read_csv(doi_pub_map, header = 0)
	try:
		df_lit_pub_map["SRA"] = fetch_study_id
		sra_path = "https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=" + fetch_study_id
		df_lit_pub_map["SRA_link"] = sra_path
		df_lit_pub_map["SRA_study_title"] = bioproject_study_title
	except:
		print "Could not retrieve SRA title"
		pass
	save_mapping(df_lit_pub_map, "doi_pmid_sra.csv", output_dir)
	os.remove(doi_pub_map)

def fetch_exp_attribute(xml_content, exp_id_item):
	#EXP META1: experiment and design metadata
	attribute_sra_exp_id = ["EXP_ID"]
	attribute_name_list = ["DESIGN_DESCRIPTION","INSTRUMENT_MODEL"]
	attribute_header_list = ["TITLE","INSTRUMENT_MODEL","DESIGN_DESCRIPTION","LIBRARY_NAME","LIBRARY_SOURCE","LIBRARY_STRATEGY","LIBRARY_SELECTION"]
	experiment_attribute_prime_head = []
	experiment_sample_prime_head = []
	experiment_primer_list = ["EXPERIMENT_ATTRIBUTE_TAG1", "EXPERIMENT_ATTRIBUTE_VALUE1","EXPERIMENT_ATTRIBUTE_TAG2","EXPERIMENT_ATTRIBUTE_VALUE2","EXPERIMENT_ATTRIBUTE_TAG3","EXPERIMENT_ATTRIBUTE_VALUE3"]
	attribute_value_list = []
	#this is the header of SRX 
	header_attributes = attribute_sra_exp_id + attribute_header_list + experiment_primer_list
	attribute_value_list.append(exp_id_item)
	experiment_primer_list_values = []
	for att in attribute_header_list:
		try:
			attribute_name = xml_content.find(att)
			attribute_value = attribute_name.get_text()
		except:
			print "Experiment attribute " + str(att) + " not found"
			attribute_value_list.append("NOT FOUND")
		else:
			attribute_value_list.append(attribute_value)

	#META2: now handle the primer list from experiment attribute
	exp_attribute_retrieve = xml_content.find_all("EXPERIMENT_ATTRIBUTE")
	if exp_attribute_retrieve:
		for tag_exp_item in exp_attribute_retrieve:
			try:
				nametag = tag_exp_item.find("TAG")
				valuetag = tag_exp_item.find("VALUE")
				nametag_value = nametag.get_text()
				valuetag_value = valuetag.get_text()
			except:
				experiment_attribute_prime_head.append("Exp_Attribute")
				experiment_primer_list_values.append("NOT FOUND")
			else:
				primer_attributes_amount = len(exp_attribute_retrieve)
				experiment_attribute_prime_head.append(nametag_value)
				experiment_primer_list_values.append(valuetag_value)
	else:
		experiment_attribute_prime_head = ["Sequencing Direction", "Primer1", "Primer2"]
		primer_value_amount = 3
		experiment_primer_list_values = ["NOT FOUND"] * primer_value_amount	

	#META3: now handle the primer list from sample attribute		
	full_head = attribute_sra_exp_id + attribute_header_list + experiment_attribute_prime_head
	exp_sample_attribute_retrieve = xml_content.find_all("SAMPLE_ATTRIBUTE")
	if exp_sample_attribute_retrieve:
		for tag_exp_item in exp_sample_attribute_retrieve:
			try:
				nametag = tag_exp_item.find("TAG")
				valuetag = tag_exp_item.find("VALUE")
				nametag_value = nametag.get_text()
				valuetag_value = valuetag.get_text()
			except:
				print "Unable to fetch"
				experiment_attribute_prime_head.append("Exp_Sample_Primer_Item")
				experiment_primer_list_values.append("NOT FOUND")
			else:
				if nametag_value == "sample_type":
					full_head.append(nametag_value)
					experiment_primer_list_values.append(valuetag_value)
				if nametag_value == "primers":
					full_head.append(nametag_value)
					experiment_primer_list_values.append(valuetag_value)
	else:
		if len(full_head) < 13:
			full_head.append("Sample_Primer_Item")
			experiment_primer_list_values.append("NOT FOUND")
	if len(full_head) < 12:
		full_head.append("Sample_Primer_Item")
		experiment_primer_list_values.append("NOT FOUND")
	concatenate_attribute_values = attribute_value_list + experiment_primer_list_values
	return concatenate_attribute_values, full_head


def search_fetch_exp_attribute():
	sra_object_file = "OutputResults/SRA_object_mapping.csv"
	output_exp_meta = "SRA_exp_metadata.csv"
	dump_exp_meta = os.path.join("OutputResults", output_exp_meta)
	df_sra_meta = pd.read_csv(sra_object_file, header=0)
	srx_object_list = []
	srx_ids = ["Experiment"]
	srx_list = df_sra_meta["Experiment"].tolist()
	print "********************************"
	print "*STEP 4: SRA Experiment mapping*"
	print "********************************"
	print "Retrieving experiment metadata..."
	print "There are " + str(len(srx_list)) + " experiments associated to this study"
	small_test_exp_items = []
	sampling_test = 10
	exp_iterator = islice(srx_list, sampling_test)
	print "For testing fetch request, I will only fetch and retrieve the first " + str(sampling_test) + " experiments associated to this study"
	for sample_item in exp_iterator:
		small_test_exp_items.append(sample_item)
	#for biosample_item in small_test_biosample_items:
	for experiment_item in small_test_exp_items:
		handle = Entrez.efetch(db="sra", id = experiment_item, format = "xml")
		out_handle = open("sra_exp_t1","wb")
		out_handle.write(handle.read())
		out_handle.close()
		handle.close()
		handle_exp_attributes = []
#retrieve element
		contents = open("sra_exp_t1", "r")
		soup = BeautifulSoup(contents, "xml") 
		exp_attributes, exp_header = fetch_exp_attribute(soup, experiment_item)

		with open(dump_exp_meta, 'a') as mf:
			exp_meta_w = csv.writer(mf)
			exp_meta_w.writerow(exp_attributes)
		mf.close()
	df_exp_meta = pd.read_csv(dump_exp_meta, names=exp_header)
	save_mapping(df_exp_meta, "Exp_meta_h.csv", "OutputResults")
	

def Main(article_doi):

	#STEP1: first mapping (doi - pmid)
	var_pmid = doi_ncbi_pmid_converter_connector(article_doi)
	print "mapping DOI - PMID saved\n"
	print "********************************"
	print "*      STEP 1 completed        *"
	print "********************************\n"
	#STEP2: second mapping(sec source mapping)
	get_sample_metadata = pmid_sra_connector(var_pmid)
	print "mapping PUBMED - SRA saved"
	print "SRA - SRA objects mapping saved"
	print "********************************"
	print "*       STEP 2 completed       *"
	print "********************************\n"
	#STEP3: third mapping (SRA BioSample)
	sra_biosample_connector(get_sample_metadata)
	print "SRA SRS BioSample mapping and meta saved"
	print "********************************"
	print "*       STEP 3 completed       *"
	print "********************************\n"
	update_lipubsra_mapping()
	#STEP4: obtain experimental metadata 
	search_fetch_exp_attribute()
	print "SRA SRX mapping and meta saved"
	print "********************************"
	print "*       STEP 4 completed       *"
	print "********************************\n"


if __name__ == '__main__':
	import sys
	try:
		article = sys.argv[1]		#doi test
	except:
		print "Usage: test_data_access1.py <DOI>"
		print "Please provide the DOI of the article you would like to retrieve information"
		print "you can try DOI: 10.1126/science.1256688 for testing"
		raise SystemExit
	else:
		if len(sys.argv) == 2:	#check if value is provided	
			Main(article)	
		else:
			print "Please follow specifications."
			raise SystemExit