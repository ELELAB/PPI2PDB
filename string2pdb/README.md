# DESCRIPTION
This script extracts filtered interactors for a specified target from the STRING database, using the STRING API. 
It provides the combined score of the confidence of each interaction (**String score**), along with 
the three components of the physical network interaction score: **Experimental score**, **Database score** and **Text mining score**. 
**Database score** > 0 threshold is being used to retrieve only reviewed interactors. 
The **String Score** threshold is adjustable, allowing the user to specify the desired confidence level for extracted interactions. 
The script also maps the interactors' STRING identifiers to their corresponding UniProt ACs using the STRING database’s human protein aliases file.

# REQUIREMENTS
Ensure you have the following installed:

- `requests`
- `pandas`
- Python 3 (on local server: use module load python/3.10/modulefile)

# ARGUMENTS
1. `<identifier>` (required, string): Uniprot Accession Code of the target protein.

2. `<aliases_file_path>` (required, string): Path to the STRING aliases file, that contains the information for mapping STRING ID to Uniprot AC. Default for running on the kbf-bioinfo01.cancer.dk the server: `/data/databases/STRING/9606.protein.aliases.v12.0.txt` 
    - **Note:** The default aliases file path applies to the kbf-bioinfo01.cancer.dk server. If running on a different server or environment, provide the correct path to a STRING aliases file. The aliases file can be downloaded from the **STRING database** and the path can be provided using the `<aliases_file_path>` argument.

3. `<threshold>` (optional, float): String_score threshold. Default for MAVISp: **String_score > 0.7** AND **Database_score** > 0 .
   To adjust the String_score threshold, specify with `-t <threshold>`.

   **Other recommended thresholds for customized searches along with the confidence of interaction:**

   | Confidence | Threshold |
   |------------|-----------|
   | Low        | 0.15      |
   | Medium     | 0.4       |
   | High       | 0.7       |
   | Highest    | 0.9       |

4. `<network>` (optional, string): 
   - Default for MAVISp: `'physical'` (Only interactors that are likely to bind to/form complexes with the target).
   - Choices: `'physical'` (subnetwork) or `'functional'` (entire network of the database).
   - To use the functional network, specify: `-n functional`.

# HOW TO RUN
1. Activate the Python environment:
   ```bash
   module load python/3.10/modulefile
   ```
2. Run the script:
   ```bash
   ./stringscore <identifier> <aliases_file_path> [-t <threshold>] [-n <network>]
   ```
## Example for MAVISp run:
### Run the bash script run.sh in `example/` folder as bash run.sh it will perform: <br />
   ```bash
   ./stringscore O15315
   ```
## Customized example of run:
### Run the bash script run.sh in `example2/` folder as bash run.sh it will perform: <br />
   ```bash
   ./stringscore Q99728 -t 0.9 -n functional
   ```
## Output:
The output is a csv file generated in the working directory and named after the input UniProt AC:
{UniProt_AC}_string_interactors.csv

The CSV file contains the following columns:
- **Target_protein**: Common name of the target protein
- **Target_Uniprot_AC**: Uniprot AC of the target protein
- **Target_id**: STRING identifier for the target protein
- **Interactor**: Common name of the interacting protein
- **Interactor_id**: STRING identifier for the interactor
- **Interactor_UniProt_AC**: Uniprot AC of the interactor
- **String_score**: Overall confidence score of the interaction
- **Experimental_score**: Score from experimental evidence
- **Database_score**: Score from database evidence
- **Textmining_score**: Score from text-mining evidence
