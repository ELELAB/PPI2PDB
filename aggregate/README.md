# Aggregate Script

## **Description**
- The **Aggregate** script integrates protein-protein interaction (PPI) data from databases Mentha and STRING into a single aggregated CSV file of ranked interactors. 
- The output includes confidence scores of interaction along with available structures from the Protein Data Bank (PDB) or the Beltrao’s dataset.
- Structures sourced via PDBminer are annotated in the column **PDBminer_structure**. If they additionally pass the PDBminer_complexes filtering, they are reported in **PDBminer_complexes_structure** column.
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

4. `-c <PDBminer_complexes_file>` (required): Path to the PDBminer_complexes output file (CSV).

5. `-o <output_filename>` (optional): Name for the output file. If not given, default name will be used, based on the Uniprot AC of the target.

---

## **How to Run**
1. **Activate Python Environment**:
   ```bash
   module load python/3.10/modulefile
   ```

2. **Execute the Script**:
   ```bash
   ./aggregate -m <path_to_mentha_file> -s <path_to_string_file> -p <path_to_PDBminer_file> -c <path_to_PDBminer_complexes_file>
   ```

Replace with the respective file paths to the Mentha, STRING, PDBminer and PDBminer_complexes outputs.


## **Example of run**
Run the bash script in the folder `example/` with `bash run.sh`. It will perform:
   ```bash
    ./aggregate -m O15315.csv -s O15315_string_interactors.csv -c O15315_filtered.csv -p O15315_all.csv 
   ```

---

## **Output**
The output is an aggregated CSV file stored in the working directory with the default filename:
```plaintext
<target_uniprot_id>_aggregated.csv
```

### **Columns in the Output File**
| **Column Name**                | **Description**                                                     |
|--------------------------------|---------------------------------------------------------------------|
| `Target_Uniprot_AC`            | UniProt ID of the target protein                                    |
| `Target_protein`               | Gene name of the target protein                                     |
| `Interactor_Uniprot_AC`        | UniProt ID of the interacting protein                               |
| `Interactor_protein`           | Gene name of the interacting protein                                |
| `Mentha_score`                 | Confidence score of the interaction from Mentha                     |
| `String_score`                 | Confidence score of the interaction from STRING                     |
| `PPI_structure`                | Available PDB IDs of the PPI structures                             |
| `PDBminer_complexes_structure` | Structures obtained from pdbminer_complexes                         |               
| `PDBminer_structure`           | Structures obtained from pdbminer                                   |               

