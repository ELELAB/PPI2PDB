# Description:
This script extracts filtered interactors for a specified target from the STRING database, using the STRING API. 
It provides the combined score of the confidence of each interaction (**String score**), along with 
the three components of the physical network interaction score: **Experimental score**, **Database score** and **Text mining score**. 
**Database score** > 0 threshold is being used to retrieve only reviewed interactors. 
The **String Score** threshold is adjustable, allowing the user to specify the desired confidence level for extracted interactions. 
The script also maps the interactors' gene names to their corresponding  primary UniProt ACs using the Uniprot API.

# Requirements:
Ensure you have the following installed:

- `requests`
- `pandas`
- Python 3 (on local server: use module load python/3.10/modulefile)

# Arguments:
1. `<identifier>` (required, string): Uniprot Accession Code of the target protein.

2. `<threshold>` (optional, float): String_score threshold. Default for MAVISp: **String_score > 0.7** AND **Database_score** > 0 .
   To adjust the String_score threshold, specify with `-t <threshold>`.

   **Other recommended thresholds for customized searches along with the confidence of interaction:**

   | Confidence | Threshold |
   |------------|-----------|
   | Low        | 0.15      |
   | Medium     | 0.4       |
   | High       | 0.7       |
   | Highest    | 0.9       |

3. `<network>` (optional, string): 
   - Default for MAVISp: `'physical'` (Only interactors that are likely to bind to/form complexes with the target).
   - Choices: `'physical'` (subnetwork) or `'functional'` (entire network of the database).
   - To use the functional network, specify: `-n functional`.

# How to run:
1. Activate the Python environment:
   ```bash
   module load python/3.10/modulefile
   ```
2. Run the script:
   ```bash
   ./stringscore <identifier> [-t <threshold>] [-n <network>]
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

### CSV File Columns:
| Column Name                | Description                                      |
|----------------------------|--------------------------------------------------|
| **Target_protein**         | Common gene name of the target protein           |
| **Target_Uniprot_AC**      | UniProt AC of the target protein                 |
| **Target_id**              | STRING identifier for the target protein         |
| **Interactor**             | Common gene name of the interacting protein      |
| **Interactor_id**          | STRING identifier for the interactor             |
| **Interactor_UniProt_AC**  | UniProt AC of the interactor                     |
| **String_score**           | Overall confidence score of the interaction      |
| **Experimental_score**     | Score from experimental evidence                 |
| **Database_score**         | Score from database evidence                     |
| **Textmining_score**       | Score from text-mining evidence                  |
