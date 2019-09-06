import requests
import sys
import os
import re
from bs4 import BeautifulSoup
from Bio import Entrez
import csv
import pandas as pd
import numpy as np
import conf_lit as na


__author__ = "airini"

"""
***WHAT the script does***
A DOI is provided in input. Initially the script will try the NCBI DOI PMID converter to retrieve the PMID
If the above will not succeed, the EBI API will be used instead
The secondary source (sequence data repository) will be searched in the PubMed xml record
If the data element above will not be found, it will be requested to the user to provide it manually
Once the sequence data identifier (SRA, SRP, BioProject) is obtained OR provided, t
he SRA object mapping will be retrieved.

The following mappings will be saved:
DOI - PMID ids and metadata associated
PMID - SRA ids 
SRA - SRS - SRP - SRX - SRR - BioSample - BioProject ids
SRR list used as input for PROFUNGIS
SRS metadata will be saved -> given for metadata curation 
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
	print "Article associated to: "
	pmid_values = retrieveElementValue('pmid', soup)
	article_titles = retrieveElementValue('title', soup)
	article_title1 = article_titles[0]
	pmid_value1 = pmid_values[0]
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
		print "Warning: path output name <OutputResults> already exists"
		print "Please change it and launch again"
		raise SystemExit
	
	article_link = doi_link + art_doi
	"""
	generate variables for requests
	I handle the error
	Here is where I create the first mapping DOI - Pubmed
	The Pubmed ID is provided, latter used for 
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
					print "First mapping 'DOI - PMID' will be now saved to output file..."
					doi_pmid_mapping2 = "doi_pmid_map.csv"
					save_mapping(update_mapping, doi_pmid_mapping2, output_dir)

				else:
					print "unable to find a Pubmed article in association to '" + art_doi + "' by using the NCBI converter"
					print "EBI API will be tried next"
					pmid_id, article_name = doi_ebi_pmid_converter(art_doi)
					doi_pmid_mapping = create_doi_pmid_header()
					ncbi_pmid_record = ncbi_pmid_link + str(pmid_id)
					checkmap = create_df_doi_pmid(doi_pmid_mapping, art_doi, article_link, pmid_id, ncbi_pmid_record, article_name)
					print "First mapping 'DOI - PMID' will be now saved to output file..."
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
	return sra_object_mapping

def biosample_retrieve(sra_object_mapping):
	"""
	In order to retrieve sample metadata, I will need the BioSample ids
	"""
	df_object_mapping = pd.read_csv(sra_object_mapping, header = 0)
	df_biosample = df_object_mapping["BioSample"]
	df_biosample.to_csv("biosample_id.csv")
	df_fetch_biosample_id = pd.read_csv("biosample_id.csv", names=["BioSample"])
	biosample_id_list = df_fetch_biosample_id["BioSample"].tolist()
	return biosample_id_list

def sra_object_retrieve(sra_id):
	"""
	I 'get' all the SRA info I need
	"""
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
	sra_map_file = save_sra_objects(complete_sra_name)
	search_sra_map_file = output_dir + '/' + sra_map_file
	biosample_id_list = biosample_retrieve(search_sra_map_file)
	return biosample_id_list
	

def msg_no_secondary_source(pub):
	print "****************************************************************************************"
	print "Key Error, I am so sorry, unfortunately the Pubmed record: <" + pub + "> did not provide a sequence source"
	print "Would you please provided it for me?"
	print "HINT: you must dig in the text of the article, I am not trained yet to retrieve it from free text"
	print "Type YES if you want to provide the SRA accession or No if you would like to come another time"
	print "No worries, all previous mappings will be found in the output folder"
	print "****************************************************************************************"

def fetch_sra_manually():
	sra_input = raw_input("Please provide SRA submission/study id (ex: SRP043706): ")
	print "checking for: " + sra_input
	try:
		biosamples = sra_object_retrieve(sra_input)
	except:
		print "Sorry I could not recognize <" + sra_input + "> and no data could be found"

def pmid_sra_connector(pubmed_id):
	"""
	here I try to slurp the pubmed article and see if I can manage to find the sequence id
	"""
	ncbi_account = na.ncbi_access1["account"]
	sample_metadata_list = ""
	pubmed_record = handle_article_data_element(ncbi_account, pubmed_id)
	for element_value in pubmed_record:
		try:
			databaseinfo = element_value['MedlineCitation']['Article']['DataBankList']
			for d in databaseinfo:
				databank = d["DataBankName"]
				if databank:
					secondarysource = "YES"
					accessionrepository = d["AccessionNumberList"]
            		accessionrepositoryname = accessionrepository[0]
            		seq_id = str(accessionrepositoryname)
        			#Here consider only SRA detection
            		find_sra = re.match(r'[SRA,SRP,SRS]', seq_id)
            		if find_sra:
            			print "it is an SRA, I can now retrieve from the NCBI SRA automatically"
            			sample_metadata_list = sra_object_retrieve(seq_id)
            		else:
            			print "Database detected. Repository name: " + databank
            			print "You can visualize data by using these accession number(s)" + seq_id

		except:
			msg_no_secondary_source(pubmed_id)
			control_input = True
			while control_input:
				sra_check = raw_input("Would you like me to try to continue the mapping and SRA data retrieval? ").lower()
				print sra_check

				if sra_check in ("yes", "y"):
					fetch_sra_manually()
					control_input = False
				elif sra_check in ("no", "n"):
					print "See you next time!"
					raise SystemExit
				else:
					print "Sorry, I did not get your input"
	return sample_metadata_list
			
			
def sra_biosample_connector(biosample_items):
	dump_sample_meta = os.path.join("OutputResults", "SRA_sample_metadata.csv")
	#PLEASE NOTE THAT some studies contain many samples to retrieve, many requests 

	########TEST 
	"""The below I use for testing a small set"""
	##handle_list_biosamples = ["SAMN02864703","SAMN02864583"]
	##################################
	#for biosample_item in handle_list_samples: #for small test
	for biosample_item in biosample_items:
		handle = Entrez.efetch(db="biosample", id = biosample_item, format = "xml")
		out_handle = open("sample_meta","wb")
		out_handle.write(handle.read())
		out_handle.close()
		handle.close()

		contents = open("sample_meta","r")
		soup = BeautifulSoup(contents, "xml")
		fetch_attribute = retrieveElementValue('Attribute',soup)
		with open(dump_sample_meta, 'a') as mf:
			sample_meta_w = csv.writer(mf)
			sample_meta_w.writerow(fetch_attribute)

def update_lipubsra_mapping():
	output_dir = "OutputResults"
	sra_obj_map = output_dir + "/SRA_object_mapping.csv"
	doi_pub_map = output_dir + "/doi_pmid_map.csv"
	try:
		df_sra_obj_map = pd.read_csv(sra_obj_map, header = 0)
		df_study = df_sra_obj_map["SRAStudy"]
		df_study.to_csv("confirm_study_id.csv")
		df_fetch_study = pd.read_csv("confirm_study_id.csv", names=["SRAStudy"])
		fetch_study_list = df_fetch_study["SRAStudy"].tolist()
		fetch_study_id = fetch_study_list[0]
	except:
		pass
	df_lit_pub_map = pd.read_csv(doi_pub_map, header = 0)
	try:
		df_lit_pub_map["SRA"] = fetch_study_id
		sra_path = "https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=" + fetch_study_id
		df_lit_pub_map["SRA_link"] = sra_path
	except:
		pass
	save_mapping(df_lit_pub_map, "doi_pmid_sra.csv", output_dir)
	os.remove(doi_pub_map)
	

def Main(article_doi):

	#first step
	var_pmid = doi_ncbi_pmid_converter_connector(article_doi)
	print "first mapping DOI - PUBMED completed"
	#second step (sec source mapping)
	get_sample_metadata = pmid_sra_connector(var_pmid)
	print "second mapping PUBMED - SRA completed"
	#third mapping, metadata retrieval
	sra_biosample_connector(get_sample_metadata)
	print "Retrieving all metadata associated to SRA"
	update_lipubsra_mapping()


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