import os
from os import listdir
from os.path import join,isfile
import tempfile,pprint
import ntpath
ntpath.basename("a/b/c")
import string,math

#################
#Common locations
#################
where_i_am = os.path.dirname(os.path.abspath(__file__))

#Where the pickled datasets are stored
pickle_datasets_location = join(where_i_am,'../data/raw')
#Where we keep the jsonified datasets.
json_datasets_location = join(where_i_am,'../data/json')
#Where we keep the metadata for searching bulk.
meta_datasets_location = join(where_i_am,'../data/meta')
#Where the original nucleotide sequences are.
nucleotides_datasets_location = join(where_i_am,'../data/nucleotides')
#Where we keep the statistics files.
stats_datasets_location = join(where_i_am,'../data/stats')
#Where we keep the gene files.
gene_datasets_location = join(where_i_am,'../data/genes')
#Where we keep the gene statistics files.
gene_stats_location = join(where_i_am,'../data/genes_stats')
#Where we keep the gene files.
geneprofile_datasets_location = join(where_i_am,'../data/gene_profiles')
#Where we keep the fasta files of entire datasets
fasta_datasets_location = join(where_i_am,'../data/fasta')
#Where we keep the fastas as formatted for mmseqs input.
mmseqs_input_location = join(where_i_am,'../data/mmseqs_input')
#Where we keep the fasta files separated by genes
fastagenes_datasets_location = join(where_i_am,'../data/fasta_genes')
#Where we keep the results of deep-search clustering.
deepsearch_datasets_location = join(where_i_am,'../data/deep_search')
#Where we keep the one off analytics files
analytics_location = join(where_i_am,'../data/analytics')
#Where we keep the cdrs data
cdrs_datasets_location = join(where_i_am,'../data/cdrs')
#Where we keep the cdrs stats data
cdrs_stats_datasets_location = join(where_i_am,'../data/cdrs_stats')
#Where we keep the results of the shallow searches.
shallow_search_results = join(where_i_am,'../data/shallow_search')
#Fasta file for each dataset.
full_fasta = join(where_i_am,'../data/full_fasta')
#Directory where we keep data for the web
web_data = join(where_i_am,'../data/web')
#Where we keep the basic data for AbStudio.
base_data_location = join(where_i_am,'../data/base_data')
#Where we keep the germlines.
germlines_location = join(where_i_am,'../data/germlines')
#Where we keep the diamond databases
diamond_dbs = join(where_i_am,'../data/diamond_dbs')

#All datasets
datasets = [f for f in listdir(json_datasets_location) if not isfile(join(json_datasets_location, f))]

##################
#Common Functions#
##################

from anarci.anarci import run_anarci

#based on IMGT number, get the region
def get_region(imgt,chain):

	if 0<=imgt and imgt<=26:
		return "fw"+chain+"1"
	if 27<=imgt and imgt<=38:
		return "cdr"+chain+"1"
	if 39<=imgt and imgt<=55:
		return "fw"+chain+"2"
	if 56<=imgt and imgt<=65:
		return "cdr"+chain+"2"
	if 66<=imgt and imgt<=104:
		return "fw"+chain+"3"
	if 105<=imgt and imgt<=117:
		return "cdr"+chain+"3"
	if 118<=imgt and imgt<=129:
		return "fw"+chain+"4"
	
#Number and transform the sequence into the format used in the db region->imgt->aa
#If get_germline is set to  true, germline is returned as well.
def number_and_transform(sequence,get_germline=False):

		#Get the numbered sequence out...
		numbered_sequence = number_sequence(sequence)
		
		germlines = {'V':'','J':'','species':''}
		if get_germline == True:
			germlines['V'] = numbered_sequence[2][0][0]['germlines']['v_gene'][0][1]
			germlines['J'] = numbered_sequence[2][0][0]['germlines']['j_gene'][0][1]
			germlines['species'] = numbered_sequence[2][0][0]['germlines']['v_gene'][0][0]
		
		#Extract the chain.		
		chain = numbered_sequence[2][0][0]['chain_type'].lower()
		if chain == 'k':
			chain = 'l'
		#Make it sortable and identifiable.
		numbered_sequence = transform_anarci_output(numbered_sequence[1][0][0][0])
		
		#transform to output as seen in the db
		organized = {}
		for elem in numbered_sequence:
			#We do not care for gaps.
			if numbered_sequence[elem] == '-':
				continue
			region = get_region(int(elem[0]),chain)
			if region not in organized:
				organized[region] = {}
			 
			organized[region][str(elem[0])+elem[1].replace(' ','')] = numbered_sequence[elem]
		if get_germline == True:
			return {'numbered':organized,'germlines':germlines}
		else:
			return organized

#Number sequence using anarci
def number_sequence(query_seq,force_human = False,assign_germline=True):
	#Number the query sequence
	#
	if force_human == True:
		res = run_anarci([('q',query_seq)],scheme='imgt',assign_germline=True,allowed_species=["human"])
	else:
		res = run_anarci([('q',query_seq)],scheme='imgt',assign_germline=True)
	
	return res

#Transform anarci list into dictionary.
def transform_anarci_output(anarci_output):
		
	sequence = {}
	for elem in anarci_output:
		
		num,insertion = elem[0]
		insertion = insertion.replace(' ','')
		if elem[1] == '-':#We do not care for gaps.
			continue
		sequence[(num,insertion)] = elem[1]
		
	return sequence

#Get the leaf of a path.
def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


#Fetch all files in directory and subdirectories.
def list_folders(directory):
	folders = []
	for _dir in listdir(directory):
		folders.append( join(directory,_dir) )
	return folders

#Fetch all files in directory and subdirectories.
def list_file_paths(directory):
   for dirpath,_,filenames in os.walk(directory):
       for f in filenames:
           yield os.path.abspath(os.path.join(dirpath, f))

#Create a temporary folder.
def create_temp_folder():
	tempdir = tempfile.mkdtemp()
	return tempdir

#Read the sequence format and return python-formatted variety
def read_sequence(data,target_region=None):
	sequence = {}
	for region in data:
		if target_region != None and region!=target_region:
			continue
		for imgt_id in data[region]:
			#See if we have insertions?
			try:
				_id = int(imgt_id)
				_id = (_id,'','')
			except ValueError:
				insertion = imgt_id[-1]
				_id = int(imgt_id[0:(len(imgt_id)-1)])
				_id = (_id,'',insertion)
			sequence[_id] = data[region][imgt_id]

	return sequence

#Transrofm given json-formatted sequence to proper sequence file.
def get_primary_sequence(data,target_region=None,with_ids=False):

	#We use imgt numbering. this means that certain entries go after insertions
	#32 31A 31B ... 33A 33B 33
	#60 60A 60B ... 61A 61B 61
	#111 111A 111B ... 112A 112B 112C
	unusual_inserts = [33,61,112]
	
	sequence = get_sorted_sequence(data,target_region=target_region)
	p_sequence = ""
	for element in sorted(sequence):
		print element
		p_sequence+= sequence[element]
	if with_ids == True:
		return sequence
	return p_sequence

#Transrofm given json-formatted sequence to proper sequence file.
def get_sorted_sequence(data,target_region=None):

	#We use imgt numbering. this means that certain entries go after insertions
	#32 31A 31B ... 33A 33B 33
	#60 60A 60B ... 61A 61B 61
	#111 111A 111B ... 112A 112B 112C
	unusual_inserts = [33,61,112]
	
	sequence = read_sequence(data,target_region=target_region)
	
	p_sequence = ""
	#Dealing with IMGT ordering.
	for uid in unusual_inserts:
		
		imgt_id = (uid,'','')
		
		if imgt_id in sequence:
			aa = sequence[imgt_id]
			del sequence[imgt_id]
			new_id = (uid,1000,'')
			sequence[new_id] = aa	
		for alpha in string.ascii_uppercase:
			imgt_id = (uid,'',alpha)
			if imgt_id in sequence:
				aa = sequence[imgt_id]
				del sequence[imgt_id]
				new_id = (uid,1000-ord(alpha),alpha)
				sequence[new_id] = aa
	
	return sequence

#Calculate sequence identity when both sequences are region->imgt->aa
#The identity is normalized to the FIRST sequence being submitted (the query).
#If target_region == None, the alignment is calculated over the entire sequence.
#It returns a json object with the identity and the absolute length difference.
#Match length set to true will force the query and template to be the same length  - works only on a region.
def sequence_identity(query,template,target_region=None,match_length=True):
	
	#Total number of matches.	
	matches = 0
	#Total number of entries we took into consideration.
	total = 0

	template_length = 0

	for region in query:
		if target_region !=None and target_region!=region:
			continue
		if match_length==True and target_region!=None and len(template[region])!=len(query[region]):
			return {'identity': 0,'length_difference':10000}
		template_length+=len(template[region])
		for imgt in query[region]:
			total+=1
			try:
				if template[region][imgt] == query[region][imgt]:
					matches+=1
			except KeyError: #Given imgt entry is not in the template.
				pass
	#print "Matches O:",matches
	#print "Template OX",template_length
	#Calculate the absolute length difference
	len_diff = abs(template_length-total)

	#Calcualte the identity.
	try:
		identity = int(100*(float(matches)/float(total)))
	except ZeroDivisionError:
		#TODO - this means there are quite erroneous sequences in there.
		return {'identity': 0,'length_difference':100000}
	return {'identity': identity,'length_difference':len_diff}

	
#sequence identity on query being imgt->aa and template region->imgt->aa 
def sequence_identity_raw_format(query,template):
	#put all the imgt entries in one pot.
	template_data = {}	
	for region in template:
		
		for imgt in template[region]:
			
			template_data[imgt] = template[region][imgt]
	
	match = 0
	total = 0
	for imgt in query:
		
		total+=1
		try:
			if query[imgt] == template_data[imgt]:
				match+=1
			
		except KeyError:
			
			pass
	perc = int(100*(float(match)/float(total)))
	return perc
	
#Get sequence identity of two anarci-numbered sequences.
#Normalize length to query and treat entries which are not found as mismatches.
def sequence_identity_parsed_format(query,template):

	matches = 0
	tot_res = 0
	for entry in sorted(query):
		
		try:
			if query[entry] == '-':
				continue

			tot_res+=1
			if template[entry] == query[entry]:
				matches+=1
			
		except KeyError:#entry not found in template
			
			pass
	identity = int(100*(float(matches)/float(tot_res)))
	return identity

#Align two sequences and print the results onto the command line.
def align_sequences(s1,s2,show=False):
	
	#Put sequences in the format we need.
	s1 = number_and_transform(s1)
	s2 = number_and_transform(s2)

	#Caalculate the sequence identity - average between t-q and t-q	
	identity_1 = sequence_identity(s1,s2)['identity']
	identity_2 = sequence_identity(s2,s1)['identity']

	identity = int((identity_1+identity_2)/2.0)

	s1 = get_sorted_sequence(s1)
	s2 = get_sorted_sequence(s2)
	
	#Get the keys.
	ks = []
	for imgt in s1:
		ks.append(imgt)
	for imgt in s2:
		if imgt not in ks:
			ks.append(imgt)
	s1_out = ''
	s2_out = ''
	ali = ''
	region_data = ''
	for imgt in sorted(ks):
		#The alignment quality.
		if imgt in s1 and imgt in s2:
			if s1[imgt] == '-' and s2[imgt] == '-':
				continue
			if s1[imgt] == s2[imgt]:
				ali+='|'
			else:
				ali+='.'
		else:
			ali+='.'
		#Amino acid identities.
		if imgt in s1:
			s1_out+=s1[imgt]
		else:
			s1_out+='-'
		if imgt in s2:
			s2_out+=s2[imgt]
		else:
			s2_out+='-'
		#Region annotation#TODO hardcoded this bit.
		reg = get_region(imgt[0],'H')
		if 'cdr' in reg:
			region_data+='^'
		else:
			region_data+=' '

	if show ==True:
		print s1_out
		print ali		
		print s2_out
		print region_data
		print "Identity=",identity,'%\n\n'
	else:
		out = s1_out+'\n'+ali+'\n'+s2_out+'\n'+region_data
		return out

#Perform sequence alignment using SW.
def sw_alignment(s1,s2):
	from Bio import pairwise2
	alignments = pairwise2.align.globalxx(s1, s2)
	matches = alignments[0][2]
	
	id_1 = int(100*(float(matches)/float(len(s1))))
	id_2 = int(100*(float(matches)/float(len(s2))))

	return int((float(id_1+id_2)/2.0))
	

#Strictly for debugging purposes.
if __name__ == '__main__':
	
	import sys
	cmd = sys.argv[1]

	if cmd == 'test_temp_folder':
		print create_temp_folder()
	if cmd == 'test_identity_calculation':
		query = "DIQMTQSPSSLSASVGDRVTITCSASSSVSYMNWYQQKPGKAPKRLIYDTSKLASGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQWSSNPPTFGGGTKVEIK"
		template = "AIQLTQSPSSLSASVGDRVTITCRASQDISSALVWYQQKPGKAPKLLIYDASSLESGVPSRFSGSESGTDFTLTISSLQPEDFATYYCQQFNSYPLTFGGGTKVEIK"
		query = number_and_transform(query)
		template = number_and_transform(template)
		print sequence_identity(template,query)
	#Debug listing folders
	if cmd == 'list_folders':
		list_folders('../data/deep_search')

	if cmd == 'get_primary_sequence':

		sequence = "QSGAEVKKPGSSVKVSCKASGYTFTNYYIYWVRQAPGQGLEWIGGINPTSGGSNFNEKFKTRVTITADESSTTAYMELSSLRSEDTAFYFCTRQGLWFDSDGRGFDFWGQGTTVTVSS"
			    
		data = number_and_transform(sequence)
		
		cdr = get_primary_sequence(data,target_region = 'cdrh3')
		full_seq = get_primary_sequence(data)
		print "CDR",cdr
		print "full_seq",full_seq
		if cdr in sequence:
			print "cdr ok"
		if full_seq == sequence:
			print "full seq ok"
	if cmd == 'sw_align':
		s1 = 'CASGYTTWERTY'
		s2 = 'CASAYDDTTWERTY'
		print sw_alignment(s1,s2)

	if cmd == 'primary_sequence':
		sequence = "QSGAEVKKPGSSVKVSCKASGYTFTNYYIYWVRQAPGQGLEWIGGINPTSGGSNFNEKFKTRVTITADESSTTAYMELSSLRSEDTAFYFCTRQGLWFDSDGRGFDFWGQGTTVTVSS"
		numbered= get_sorted_sequence(number_and_transform(sequence))
		print len(numbered),len(sequence)
		for elem in sorted(numbered):
			print elem,numbered[elem]
		print get_primary_sequence(number_and_transform(sequence))
	

