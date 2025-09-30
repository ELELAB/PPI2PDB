#Requirements
#mentha2pdb.py
#target_uniprot_ID.txt  #including uniprot ID of protein of interest

uniprot=$1

echo $1 > target_uniprot_ID.txt

#load python3.10
module load python/3.10
#symbolic link to local version of Mentha to use
ln -s /data/databases/mentha-20250428/2025-04-28

python ../mentha2pdb.py -i 2025-04-28 -t target_uniprot_ID.txt -s 0.2 -o $1.csv -p -a  -extra /data/databases/AF_Huri_HuMAP/summary/huri_upac.csv  /data/databases/AF_Huri_HuMAP/summary/humap_upac.csv -af /data/databases/AF_Huri_HuMAP -ec 0.2

