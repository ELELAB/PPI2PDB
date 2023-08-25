module load python/3.7/modulefile
python ../mentha2pdb.py -i /data/databases/mentha-20230417/2023-04-17 -t target_uniprot_ID.txt -s 0.2 -o out.csv -p -a -c config.ini -extra /data/databases/AF_Huri_HuMAP/summary/HuRI.csv  /data/databases/AF_Huri_HuMAP/summary/humap.csv
