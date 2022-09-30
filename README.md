# Mentha2PDB

#Requirements

python3 (on local server: use module load python/3.7/modulefile) <br />
pypdb v2.2 (only version 2.2 works since it has been recently updated to support the RCSB Search API V2) <br />
pandas <br />

#Usage

python mentha2pdb.py -i /data/databases/mentha-20220530/2022-05-30 -t target_uniprot_ID.txt -s 0.2 -o out.csv -p <br />
-i mentha database to analyze <br />
-t input file with list of the uniprot IDs of the target proteins <br />
-s threshold of mentha score for filtering (i.e., remove all the entries with mentha score below the threshold) <br />
-o name output file <br />
-p include in the ouput the PMID of the relevant publications related to the interaction. The "unassigned<#code>" are broken publication annotations in the Mentha database generally coming from missing annotations in IntAct <br />
-x generate a single csv output file per each target uniprot ID named dataframe_<target_uniprot_ID>.csv (this option overrides option -o) <br />
-a  have in output input files for AlphaFold_multimer <br />

In case of incorrect or obsolete Uniprot ID or gene names annotations present in Mentha database, mentha2pdb write a log file reporting them, please check the log file carefully. 

#Example
run the script in example as <br />
python mentha2pdb.py -i /data/databases/mentha-20220530/2022-05-30 -t target_uniprot_ID.txt -s 0.2 -o out.csv -p -a <br />
