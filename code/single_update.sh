#Updating a single dataset. Assuming we already have the .json files.

dataset=$1

#Creaate the folder structure.
echo "Creating data structure for  $dataset";
python ProcessData.py parallel_datasets

#Split into genes in genes and create basic stats in stats
echo "Getting gene stats for $dataset";
python ProcessData.py read_single_dataset $dataset

#Create the gene counts in genes_stats.
echo "Preforming further analytics on $dataset";
python Analytics.py gene_counts $dataset

#Calculate cdr stats that are in data/cdrs_stats
echo "Calculating CDR stats for $dataset";
python CDR.py calculate_stats $dataset
