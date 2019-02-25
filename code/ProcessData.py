#System imports.
import cPickle as pickle
import os
from os import listdir
from os.path import isfile, join
import json
import gzip

#Imports of our own code
from Common import json_datasets_location,gene_datasets_location



#Split files by genes.
def process_file(ds_path,ds_file,ds_name):

	#Where the gene outputs go.
	gene_out = join(gene_datasets_location,ds_name)
	
	data_name = ds_file.replace('.json.gz','')
	
	#Global genes.
	genes = []
	meta_line = True
	for line in gzip.open(join(ds_path,ds_file),'rb'):
		if meta_line == True:
			meta_line=False
			continue
		sequence = json.loads(line)
		
		metadata = {'redundancy':int(sequence['redundancy']),'v':sequence['v'],'j':sequence['j']}
			
		#######
		#GENES#
		#######
		#Split by genes and write out as we go along.
		gene_id = metadata['v']+'@'+metadata['j']

		#Create the folder.
		gene_folder = metadata['v'].split('*')[0]+metadata['j'].split('*')[0]
		gene_folder = join(gene_out,gene_folder)
		if not os.path.exists(join(gene_folder)):
			os.mkdir(gene_folder)

		#Split by file (isotype, organism)
		gene_file = join(gene_folder,ds_file.replace('.json.gz','')+'--'+gene_id)

		if gene_file not in genes:
			genes.append(gene_file)
		
		gene_file = open(gene_file,'a')
		gene_file.write(json.dumps(sequence)+'\n')
		gene_file.close()

		

	#Write out the genes.
	print "Writing out genes..."
	for gene_file in genes:
		print "Compressing ",gene_file
		os.system('gzip -9 '+gene_file)
	
	

#Load given dataset by name
def read_dataset(ds_name):

	ds_path = join(json_datasets_location,ds_name)

	ds_files = [f for f in listdir(ds_path) if isfile(join(ds_path, f))]
	#Load all the files from this dataset.
	i = 0
	for ds_file in ds_files:
		i+=1
		print "Reading",ds_file,'[',i,'/',len(ds_files),']'
		process_file(ds_path,ds_file,ds_name)
		
#Go through all available datasets and process them	
def read_available_datasets():
	
	datasets = [f for f in listdir(json_datasets_location) if not isfile(join(json_datasets_location, f))]
	
	for ds in datasets:
		print "Current ds",ds
		read_dataset(ds)

#Go through all available datasets and process them	
def parallel_read_available_datasets():

	datasets = [f for f in listdir(json_datasets_location) if not isfile(join(json_datasets_location, f))]
	
	for ds in datasets:
		
		#Create the gene directory if necessary
		gene_out = join(gene_datasets_location,ds)
		if not os.path.exists(gene_out):
			os.mkdir(gene_out)

		
		print "python ProcessData.py read_single_dataset ",ds,'&'

#Create parallel script for a single file.
def parallel_script_single_ds(ds):
	ds_path = join(json_datasets_location,ds)

	ds_files = [f for f in listdir(ds_path) if isfile(join(ds_path, f))]
	#Load all the files from this dataset.
	i = 0
	for ds_file in ds_files:
		i+=1
		print "python ProcessData.py process_single_file ",ds_path,ds_file,ds,'&'
		
if __name__ == '__main__':

	import sys
	cmd = sys.argv[1]
	
	#Debugging processing single files.
	if cmd == 'debug_file_processing':
		ds_path = "/data/colocolo1/not-backed-up/krawczyk/AbMap/data/json/22_Greiff_mouse/"
		ds_file = 'H_MOUSE_IGHM_Greiff_NUMBERED_part_42.json.gz'
		ds_name = '22_Greiff_mouse'
		process_file(ds_path,ds_file,ds_name)
	if cmd == 'process_single_file':
		ds_path = sys.argv[2]
		ds_file = sys.argv[3]
		ds_name = sys.argv[4]
		process_file(ds_path,ds_file,ds_name)
		pass
	#Create the parallel script for a single dataset
	if cmd == 'parallel_single_ds':
		#Create parallelized file.
		parallel_script_single_ds(sys.argv[2])
	#Create the parallel script.
	if cmd == 'parallel_datasets':
		#Create parallelized file.
		parallel_read_available_datasets()
	if cmd == 'read_single_dataset':
		#Read a single dataset
		ds = sys.argv[2]
		read_dataset(ds)
		
	if cmd == 'serial':
		
		#Read all available data, split up by gene and collect statistics.
		read_available_datasets()
	
	
