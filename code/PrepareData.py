#System imports.
import cPickle as pickle
import os,csv,sys
csv.field_size_limit(sys.maxsize)
from ast import literal_eval
from subprocess import check_call
from os import listdir
from os.path import isfile, join
import json,shutil,time,gzip

#Imports of our own code
from Common import pickle_datasets_location,json_datasets_location,list_file_paths,meta_datasets_location,path_leaf,nucleotides_datasets_location

#Load given dataset by name
def load_dataset(ds_name):

	print ds_name

	ds_path = join(pickle_datasets_location,ds_name)
	#Load all the files from this dataset.


	for ds_file in list_file_paths(ds_path):

		ds_file_name = ds_file.split('/raw/')[1].replace('/','_')

		#Ignore the entire big files.
		if 'Summary' in ds_file:
			continue
		if '.csv' not in ds_file:
			continue

		#Make sure that there are some sequences in there.
		f_size = int(os.stat(ds_file).st_size) 
		if f_size == 0 :
			continue
		
		#check that there is the corresponding .json file with metadata
		metadata = ds_file.replace('.csv','.json')
		#TODO -- some files do not have associated metadata -- run the Consistency.py step.
	    	if not os.path.exists(metadata):
			print 'Missing metadata at',metadata
            		continue
		#print "Examining",ds_file
		#Copy metadata where it should be.
		meta_path = join(meta_datasets_location,ds_name)
		if not os.path.exists(meta_path):
			os.mkdir(meta_path)
		shutil.copyfile(metadata,join(meta_path,ds_file_name.replace('.csv','.json')))
		
		#Read the metadata.
		meta_json = json.load(open(join(meta_path,ds_file_name.replace('.csv','.json'))))


		#check that there is the corresponding .fasta nucleotide file
		nucleotidesdata = ds_file.replace('.csv','.fasta')
		#TODO -- some files do not have associated metadata.
	    	if not os.path.exists(nucleotidesdata):
			print 'Missing nucleotides at',nucleotidesdata
            		continue
		
		#Copy nucleotides where it should be.
		nucleotides_path = join(nucleotides_datasets_location,ds_name)
		if not os.path.exists(nucleotides_path):
			os.mkdir(nucleotides_path)
		nucleotide_file = join(nucleotides_path,ds_file_name.replace('.csv','.fasta'))
		shutil.copyfile(nucleotidesdata,nucleotide_file)
		

		
		#New json file.
		#See if we have this folder?
		out_ds_path = join(json_datasets_location,ds_name)

		json_file = join(out_ds_path,ds_file_name.replace('.csv','.json'))
		#Do not overwrite if it already exists.
		if os.path.exists(json_file+'.gz'):
			#print "We already have",ds_file
			continue
		

		#Read in the file.
		if not os.path.exists(out_ds_path):
			os.mkdir(out_ds_path)
		
		#Write out the json file.
		#print "Jsonifying..."
		out = open(json_file,'w')

		#The first line is the header with metadata.
		out.write(json.dumps(meta_json)+'\n')
		
		#Naming sequences
		seq_ordinal = 0
		with open(ds_file, 'rb') as csvfile:
			
			lines = csv.reader(csvfile, delimiter=',', quotechar='"')
			try:
				for line in lines:
					seq_ordinal+=1
					#This is the current reading
				
					fasta,redundancy,cdr3,num_errors,errors,vgene,jgene,numbered,seq_name = line
			
	    				data_element = {'name':seq_ordinal,'data':numbered,'v':vgene,'j':jgene,'cdr3':cdr3,'redundancy':int(redundancy),'seq':fasta,'num_errors':num_errors,'errors':errors,'original_name':literal_eval(seq_name)[0]}
					json_line = json.dumps(data_element)
					out.write(json_line+'\n')
			except:
				print ds_file
				print "Problematic ds",ds_name
				raise
		out.close()
		#print "Compressing Json..."
		os.system('gzip '+join(out_ds_path,json_file))
		#print "Compressing Nucleotide..."
		os.system("gzip "+nucleotide_file)
		
		
#Go through all available datasets and make them in the correct representation.		
def read_available_datasets():
	
	datasets = [f for f in listdir(pickle_datasets_location) if not isfile(join(pickle_datasets_location, f))]
	
	for ds in datasets:
		print "Current ds",ds
		load_dataset(ds)
#Create a script to do everything in parallel.
def create_parallel_script():
	datasets = [f for f in listdir(pickle_datasets_location) if not isfile(join(pickle_datasets_location, f))]
	
	for ds in datasets:
		print "python PrepareData.py do_single ",ds,'&'
		
if __name__ == '__main__':
	
	import sys

	cmd = sys.argv[1]

	
	#Creating a parallel script to prepare all data.
	if cmd == 'create_parallel':
		create_parallel_script()
	
	#Preparing data from a single dataset.
	if cmd == 'do_single':
		ds = sys.argv[2]
		load_dataset(ds)

	if cmd == 'all':
		pass

		#Read the pickle files, translate them into gzip-json format.
		read_available_datasets()
	if cmd == 'debug':
		#Debugging doing single dataset.
		load_dataset('Galson_2015')
	
