#!/usr/bin/env python3

# STRING2PDB 
# Copyright (C) 2024  Eleni Kiahaki and Matteo Tiberti, Cancer Structural Biology, Danish Cancer Institute
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import requests
import argparse
import pandas as pd
from io import StringIO

def main():
    parser = argparse.ArgumentParser(
        description="Retrieval of interaction data from the STRING database for a given gene name."
    )
    parser.add_argument(
        "identifier", 
        type=str,
        help="HUGO Gene name to retrieve interactors for."
    )
    parser.add_argument(
        "-t", 
        "--threshold",
        type=float,
        default=0.7,
        help="Minimum STRING confidence score for interaction filtering (default: 0.7)."
    )
    parser.add_argument(
        "-n",
        "--network",
        type=str,
        default="physical",
        choices=["functional", "physical"],
        help="STRING network type to be used: 'physical' for physical interactions, 'functional' for all interactions (default: 'physical')."
    )

    args = parser.parse_args()
    
    base_url = "https://string-db.org/api/tsv/get_string_ids"
    params = {
        'identifier': args.identifier,  
        'species': 9606,         
        'limit': 0,  
        'caller_identity': "MAVISp_web_app"     
    }
    try:
        string_response = requests.get(base_url, params=params)
        string_response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error: Unable to get data ({e})")
        return

    # Read tsv data into pandas dataframe:
    data = pd.read_csv(StringIO(string_response.text), sep='\t')

    # Check if identifier is in data:
    if args.identifier not in data['preferredName'].values:
        print(f"Error: No STRING identifier found for {args.identifier} in the results.")
        return

    # Handle case of multiple STRING IDs found:
    if data.shape[0] > 1:
        print(f"Warning: Multiple STRING IDs found for {args.identifier}. Selecting STRING ID based on preferred name.")
        matching_row = data[data['preferredName'] == args.identifier]
        string_id = matching_row.iloc[0]['stringId']
    else:
        string_id = data.iloc[0]['stringId']

    # Get interactors from STRING API
    interactors_url = "https://string-db.org/api/tsv/interaction_partners"
    interactors_params = {
        'identifier': string_id,  
        'species': 9606,
        'required_score': args.threshold,
        'limit': 0,
        'network_type': args.network,
        'caller_identity': "MAVISp_web_app"
    }

    try:
        interactors_response = requests.get(interactors_url, params=interactors_params)
        interactors_response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error: Unable to get interaction data ({e})")
        return

    # Read tsv data into pandas dataframe:
    interactors_df = pd.read_csv(StringIO(interactors_response.text), sep='\t')

    # Filter interactors:
    filtered_interactors = interactors_df[interactors_df['score'] >= args.threshold]

    if filtered_interactors.empty:
        print(f"No interactors found with score > {args.threshold}.")
        return

    # Extract columns:
    interactors = filtered_interactors[['preferredName_A', 'stringId_A', 'preferredName_B', 'stringId_B', 
                                        'score', 'escore', 'dscore', 'tscore']]

    # Output CSV file:
    output_file = f"{args.identifier}_string_interactors.csv"
    interactors.to_csv(output_file, index=False, header=['Target_protein', 'Target_id', 'Interactor', 'Interactor_id', 
                                                         'String_score', 'Experimental_score', 'Database_score', 
                                                         'Textmining_score'])

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
