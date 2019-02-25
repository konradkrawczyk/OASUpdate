from os import listdir
import json,os
from os.path import join,isfile
from Common import pickle_datasets_location,list_file_paths

random_elements = "abcdefghijklmn"

if __name__ == '__main__':
	
	import sys
	#The dataset name we are dealing with - should correspond to the dataset name in ../data/raw
	dsname = sys.argv[1]

	for element in listdir(pickle_datasets_location):

		if element != dsname:
			continue
		if not isfile(join(pickle_datasets_location,element)):
			print '->',element
			author = ""
			for f in list_file_paths(join(pickle_datasets_location,element)):
				if '.json' in f:
					
					json_data = json.load(open(f))
					if 'Size' not in json_data:
						print "No size infomration", f
						quit()
					#Format name
					author = json_data['Author'].replace('et al.','').replace(' ','_').replace('(','').replace(')','').replace(',','')
					while '__' in author:
						author = author.replace('__','_')
			#Move the file.
			print "Moving..."
			i = 0
			source_location= join(pickle_datasets_location,element)
			target_location = join(pickle_datasets_location,author)
			while os.path.exists(target_location):
				target_location+=random_elements[i]
				i+=1
			os.system('mv '+source_location+' '+target_location)
