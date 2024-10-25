# DESCRIPTION
This script extracts filtered interactors for a specified target from the STRING database. 
It provides the combined score of the confidence of each interaction (**String_score**), 
along with the three components of the physical network score: **Experimental_score**, **Database_score** and **Text_mining score**. 
**String_score** is used as the threshold for filtering.

# REQUIREMENTS
Ensure you have the following installed:

- `requests`
- `pandas`
- Python 3 (on local server: use module load python/3.10/modulefile)

# ARGUMENTS
1. `<identifier>` (required, string): HUGO Gene Name of the target protein.

2. `<threshold>` (optional, float): String_score threshold. Default for MAVISp: **0.7**. 
   To adjust the threshold, specify with `-t <threshold>`.

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

# HOW TO RUN
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
   ./stringscore RAD51B
   ```
## Customized example of run:
### Run the bash script run.sh in `example2/` folder as bash run.sh it will perform: <br />
   ```bash
   ./stringscore BARD1 -t 0.9 -n functional
   ```
## Output:
The output is a csv file generated in the working directory and named after the input HUGO Gene Name:
{HUGO_NAME}_string_interactors.csv

The CSV file contains the following columns:
- **Target_protein**: Common name of the target protein
- **Target_id**: STRING identifier for the target protein
- **Interactor**: Common name of the interacting protein
- **Interactor_id**: STRING identifier for the interactor
- **String_score**: Overall confidence score of the interaction
- **Experimental_score**: Score from experimental evidence
- **Database_score**: Score from database evidence
- **Textmining_score**: Score from text-mining evidence
