# DESCRIPTION:
# Script that extracts filtered interactors of specified target from STRING database
# Provides the combined score of the confidence of each interaction (String_score)
# along with its 3 components (Experimental score, Database score, Text mining score)
# String_score is used as the threshold for filtering.

# REQUIREMENTS:
# python 3
# requests library

#PARAMETERS:

# 1. <identifier> (required): Uniprot accession number of target protein

# 2. <threshold> (optional): String_score threshold. Default for MAVISp: 0.7.

# Other recommended thresholds for customized searches:

# *Confidence*  *Threshold*    *Flag*    
# low               0.15       -t 150
# medium            0.4        -t 400
# high              0.7        -t 700 
# highest           0.9        -t 900

# 3. <network> (optional):
# Default for MAVISp: 'physical' (Only interactors that are likely to bind to/form complexes with target)
# Choices: 'physical' (subnetwork) or 'functional' (entire network of database) 
# TO use the functional, specify: -n functional

#HOW TO RUN:

# 1. Activate python env:
module load python

# 2. Run the stringscore.py:
python stringscore.py <identifier> [-t <threshold>] [-n <network>]

#Example for MAVISp run:
python stringscore.py O15315

#Other example:
python stringscore.py O15315 -t 900 -n functional

# OUTPUT:
# tsv file in working dir:
# {uniprot_ac}_string_interactors.tsv

# Example:
# O15315_string_interactors.tsv







