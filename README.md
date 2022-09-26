# Mentha2PDB

#Requirements

python3 (on local server: use module load python/3.7/modulefile) <br />
pypdb v2.2 (only version 2.2 works since it has been recently updated to support the RCSB Search API V2) <br />
pandas <br />

#Usage

python mentha2pdb.py -i /data/databases/mentha-20210301/2021-03-01 -t target_uniprot_ID.txt -s 0.2 -o out.csv <br />
-i mentha database to analyze <br />
-t input file with list of the uniprot IDs of the target proteins <br />
-s threshold of mentha score for filtering (i.e., remove all the entries with mentha score below the threshold) <br />
-o name output file <br />

#Example
run the script in example as <br />
python mentha2pdb.py -i /data/databases/mentha-20210301/2021-03-01 -t target_uniprot_ID.txt -s 0.2 -o out.csv <br />
