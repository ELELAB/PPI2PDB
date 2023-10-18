# Mentha2PDB

#References

Arnaudi M, Beltrame L, Degn K, Utichi M, Pettenella A, Scrima S, et al. MAVISp: Multi-layered Assessment of VarIants by Structure for proteins. bioRxiv. 2022. https://doi.org/10.1101/2022.10.22.513328.

contacts: matl-at-cancer.dk, elenap-at-cancer.dk, elpap-at-dtu.dk

#Requirements

python3 (on local server: use module load python/3.7/modulefile) <br />
pypdb v2.2 (only version 2.2 works since it has been recently updated to support the RCSB Search API V2) <br />
pandas <br />

#Usage

python mentha2pdb.py -i /data/databases/mentha-20220530/2022-05-30 -t target_uniprot_ID.txt -s 0.2 -o out.csv -p -extra /data/databases/AF_Huri_HuMAP/summary/HuRI.csv  /data/databases/AF_Huri_HuMAP/summary/humap.csv <br />
-i mentha database to analyze <br />
-t input file with list of the uniprot IDs of the target proteins <br />
-s threshold of mentha score for filtering (i.e., remove all the entries with mentha score below the threshold) <br />
-o name output file <br />
-p include in the ouput the PMID of the relevant publications related to the interaction. The "unassigned<#code>" are broken publication annotations in the Mentha database generally coming from missing annotations in IntAct <br />
-x generate a single csv output file per each target uniprot ID named dataframe_<target_uniprot_ID>.csv (this option overrides option -o) <br />
-a have in output input files for AlphaFold_multimer <br />
-c Config file containing manual annotations of PDBs or pair of partners not included in the mentha db to be annotated in the final output <br /> 
-extra AlphaFold2 dimeric complexes databases (HuRI.csv and humap.csv datasets) from Burke, D.F. et al.  Nat Struct Mol Biol 30, 216–225 (2023). https://doi.org/10.1038/s41594-022-00910-8  <br />

In case of incorrect or obsolete Uniprot ID or gene names annotations present in Mentha database, mentha2pdb write a log file reporting them, please check the log file carefully.
The argument -c can be used to give mentha2pdb in input a configuration .ini file with pairs of partners for which there are information in literature about their interactions but they are not automatically retrived since they are not yet in the menths database, there are issues in the annotation of the experimental structure (i.e. PDB with fusion constructs) or experimental structure not released yet. The entries from the configuration file should include the following format: 

[Q9GZQ8]
interactor1= Q9Y4G2,PLEKHM1,3X0W,25498145

[Q7Z3C6]
interactor1 = Q2TAZ0,ATG2A,,36347259
interactor2 = O75143,ATG13,,34369648

Where [Q9GZQ8] and [Q7Z3C6] are the uniprot IDs of the target proteins. The interactors to annotate should be listed one after the other, one per row. For each interactor is mandatory to specify its Uniprot ID, gene name and PMID associated with the publication in wich the interaction is experimentally validated. Optional a PDB could be specified and it will be annotated.

Furthermore the -extra and the -af arguments can be used to annotate dimeric complexes generated with AlphaFold2 from the HuRI and Hu.Map databases (HuRI.csv and humap.csv datasets) from Burke, D.F. et al. 2023,  Nat Struct Mol Biol 30, 216–225 (https://doi.org/10.1038/s41594-022-00910-8).The -extra argument add two columns "HuRI" and "humap" at the end of the output csv file that is annotated if a model of the complex of the target and interactor has been generated with confidence score (called pDockQ score) higher than a cutoff (by default the cutoff is set to pDockQ > 0.5 since it is reported to define high-confidence models). 
The -ec argument can be used to set a different pDockQ cutoff than the default one to filter the models of the complexes. 
The -af argument allows the script to generate a local copy of the corresponding files of the filtered models in the folder AF_Huri_HuMAP. 


#Example
run the bash script run.sh in example folder as bash run.sh it will perform <br />
python ../mentha2pdb.py -i /data/databases/mentha-20230417/2023-04-17 -t target_uniprot_ID.txt -s 0.2 -o out.csv -p -a <br />

run the bash script run.sh in example2 folder as bash run.sh it will perform <br />
python ../mentha2pdb.py -i /data/databases/mentha-20230417/2023-04-17 -t target_uniprot_ID.txt -s 0.2 -o out.csv -p -a -c config.ini -extra /data/databases/AF_Huri_HuMAP/summary/HuRI.csv  /data/databases/AF_Huri_HuMAP/summary/humap.csv -af /data/databases/AF_Huri_HuMAP -ec 0 <br />

run the bash script run.sh in example3 folder as bash run.sh it will perform <br />
python ../mentha2pdb.py -i /data/databases/mentha-20230417/2023-04-17 -t target_uniprot_ID.txt -s 0.2 -o Q9GZQ8.csv -p -a -c config.ini -extra /data/databases/AF_Huri_HuMAP/summary/HuRI.csv  /data/databases/AF_Huri_HuMAP/summary/humap.csv -af /data/databases/AF_Huri_HuMAP/ -ec 0.2 <br />
