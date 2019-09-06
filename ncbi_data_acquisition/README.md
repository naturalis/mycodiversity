# ncbi_data_acquisition
retrieving NCBI PubMed, SRA mappings and sample metadata associated to a scientific publication

# lit_pub_sra_mapper.py
Literature PubMed SRA mapper 
The lit_pub_sra_mapper has been created to requesting data from the NCBI databases Pubmed and SRA.
It provides all the mappings associated to:
- a scientific article (DOI), 
- the reference NCBI PubMed record 
- cross reference mapping to NCBI SRA. 
It also provides a SRA object mapping to SRS (Samples), SRR (Sequence Run) and SRX (Sequence Experiment), NCBI BioProject and BioSample database mappings. 
It additionally downloads all the metadata associated to the SRA sample(s).

The lit_pub_sra_mapper.py retrieves data automatically as long as data is provided in the given repository requests.
lit_pub_sra_mapper.py will only 'request' help from the user when it is not able to localize any sequence source (SRA) in the PubMed record.
In this case, the user has the opporturnity to help, thus provide the SRA accession and the script can continue retrieving data.
The user is also able to terminate at this stage, thus not provide the information requested. 
In this case, all the previous outputs and mappings are still provided, except for the SRA data. 

## How to run lit_pub_sra_mapper.py
  launch the script followed by a DOI
  >python lit_pub_sra_mapper.py <DOI>
	
  ex: python lit_pub_sra_mapper.py 10.1126/science.1256688

### Configuration 

 - provide in the conf_lit.py file your personal NCBI account (ncbi_access["account"]) and save file
 - please follow the information 'API key management' of the NCBI link:"https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/"
 
### Requirements

Please check the requirements.txt file for the script's dependencies 

### Cases

The script takes into account four main scenarios:

CASE1:

	DOI is not able to retrieve data with the NCBI DOI - PMID converter, tries with EBI API and its PubMed record provides a Sequence source  (aka 'DataBankName' element)

CASE2:

	DOI is not able to use the NCBI DOI - PMID converter, uses EBI API and its PubMed record does not provide a sequence source

CASE3:

	DOI is able to retrieve data with the NCBI DOI - PMID converter, but PubMed record contains no sequence source

CASE4:

	DOI is able to retrieve data with the NCBI DOI - PMID converter, and PubMed record provides the sequence source

#Examples to use with the 4 Cases

CASE1:

DOI: 10.1126/science.1256688
SRP: SRP043706

CASE2:

DOI: 10.1111/j.1365-294X.2010.04898.x
SRP: SRP102378

CASE3:

DOI: 10.1038/ismej.2013.28
SRP: SRP015917


CASE4:
DOI: 10.1073/pnas.1508382112
SRP:  SRP052716


For CASE1 and CASE3, the script retrieves data automatically, without the need to wait for an input information from user
For CASE2 and CASE4, the user will need to provide the SRA accession manually in order to continue with all the mappings and metadata retrieval

###MAPPING PIPELINE METHODS

doi_ncbi_pmid_converter_connector : will request and retrieve PubMed elements by using the NCBI and EBI API

pmid_sra_connector : will retrieve data from sra by using the PubMed elements

sra_biosample_connector : the cross reference sra biosample will retrieve the sample metadata associated with the SRA study


###OUTPUTS

The following mapping files and outputs will be generated 

>doi_pmid_sra.csv -> contains the mappings and links from article, pubmed and sra

>Profungis_accession_list.txt -> the input file requested for PROFUNGIS pipeline is generated. It contains the sequence run accession list of the raw fastq files.

>SRA_object_mapping.csv -> very useful file which contains all possible sra object mappings associated with a study (SRA, SRP, SRS, SRR, SRX, BioSample and BioProject)

>SRA_sample_metadata.csv -> the sample metadata associated with the sample collected

>SRA_submission_meta.csv -> the metadata associated with the SRA submission record





