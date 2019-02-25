#Updating a single dataset. Assuming we already have the .json files.

dataset=$1

#Move the data to the format we want.
python NameDatasets.py $dataset

#Remove any previous compressed files.
rm -f json.tar.gz nucleotides.tar.gz meta.tar.gz
#Now compress...
cd ../data
tar -zcvf data.tar.gz json nucleotides meta

#scp data.tar.gz [path on the server]
