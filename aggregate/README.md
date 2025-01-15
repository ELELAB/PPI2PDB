# Aggregate Script

## **Description**
- The **Aggregate** script integrates protein-protein interaction (PPI) data from databases Mentha and STRING into a single aggregated CSV file of ranked interactors. 
- The output includes confidence scores of interaction along with available structures from the Protein Data Bank (PDB) or the Beltrao’s dataset.
- Structures not extracted from PPI databases but sourced via PDBminer are annotated as:
  - `**`: Filtered by **PDBminer_complexes**.
  - `*`: Sourced directly from **PDBminer** (not filtered by PDBminer_complexes).
- The interactors are ranked based on the Mentha confidence score whereas STRING confidence score is the secondary ranking criterion.

## **Requirements**
### **Software:**
- Python 3.x 
- Python packages:
  - `pandas`

### **Files:**
#### For the same target protein: ####
- A Mentha2PDB (`mentha2pdb.py`) output CSV
- A STRING2PDB (`string2pdb`) output CSV
- A PDBminer output csv
- A PDBminer_complexes (`find_PDBminer_complexes.py`) output csv


## **Arguments**
1. `-m <mentha_file>` (required): Path to the Mentha2PDB output file (CSV). 
   
2. `-s <string_file>` (required): Path to the STRING2PDB output file (CSV).

3. `-p <PDBminer_file>` (required): Path to the PDBminer output file (CSV).

4. `-pc <PDBminer_complexes_file>` (required): Path to the PDBminer_complexes output file (CSV).

5. `-o <output_filename>` (optional): Name for the output file. If not given, default name will be used, based on the Uniprot AC of the target.

---

## **How to Run**
1. **Activate Python Environment**:
   ```bash
   module load python/3.10/modulefile
   ```

2. **Execute the Script**:
   ```bash
   ./aggregate -m <path_to_mentha_file> -s <path_to_string_file> -p <path_to_PDBminer_file> -pc <path_to_PDBminer_complexes_file>
   ```

Replace with the respective file paths to the Mentha, STRING, PDBminer and PDBminer_complexes outputs.


## **Example of run**
Run the bash script in the folder `example/` with `bash run.sh`. It will perform:
   ```bash
    ./aggregate -m O15315.csv -s O15315_string_interactors.csv -p O15315_all.csv -pc O15315_filtered.csv
   ```

---

## **Output**
The output is an aggregated CSV file stored in the working directory with the default filename:
```plaintext
<target_uniprot_id>_aggregated.csv
```

### **Columns in the Output File**
| **Column Name**           | **Description**                                                      |
|----------------------------|---------------------------------------------------------------------|
| `target uniprot id`        | UniProt ID of the target protein                                    |
| `target uniprot gene`      | Gene name of the target protein                                     |
| `interactor uniprot id`    | UniProt ID of the interacting protein                               |
| `interactor uniprot gene`  | Gene name of the interacting protein                                |
| `mentha score`             | Confidence score of the interaction from Mentha                     |
| `string score`             | Confidence score of the interaction from STRING                     |
| `structure`                | Available PDB IDs of the PPI structures                             |
| `other resources`          | Report structures obtained from pdbminer/ pdbminer_complexes        |               

