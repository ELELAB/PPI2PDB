# Mentha2PDB

#Requirements
python3 use module load python/3.7/modulefile
pypdb v2.1
pandas

#Usage
python mentha2pdb.py -i /data/databases/mentha-20210301/2021-03-01 -t target_uniprot_ID.txt - s 0.2 -o out.csv
-i mentha database to analyze
-t input file with list of the uniprot IDs of the target proteins
-s threshold of mentha score for filtering (i.e., remove all the entries with mentha score below the threshold)
-o name output file

#Example
run the script in example as
python mentha2pdb.py -i /data/databases/mentha-20210301/2021-03-01 -t target_uniprot_ID.txt - s 0.2 -o out.csv
