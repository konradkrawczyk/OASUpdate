# OASUpdate

Scripts to update OAS.

# How-To.

Put the raw data in data/raw. e.g. data/raw/Sheng

Next run the script from the code directory specifying which raw dataset you want to process:

`./single_update.sh Sheng`

This should have created a data.tar.gz file in the data directory that contains the meta, json and nucleotide data on the new ds.

data.tar.gz should be uploaded to the folder on the server that has the same folder structure (json, nucleotides, meta).

Uncompress data.tar.gz in that directory and the folders with the new ds should appear next to these (Check!). 

Scripts to update the data server-side should still be there but for completeness, they should achieve the following:

Remove the old .tar.gz files for:

`meta.tar.gz`

`nucleotidies.tar.gz`

`json.tar.gz`

These files are the ones that are linked to the bulk upload.

There should be a script to recreate these json files from the json, nucleotides and meta folders - hence including the new ds.

