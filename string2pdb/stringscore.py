#!/usr/bin/env python3

import requests
import argparse
import csv
import pandas as pd
from io import StringIO

def get_interactors(string_id, threshold, network):
    base_url = "https://string-db.org/api/tsv/interaction_partners"
    params = {
        'identifier': string_id,  
        'species': 9606,
        'required_score': threshold,
        'limit': 0,
        'network_type': network,
        'caller_identity': "MAVISp_web_app"
    }

    try:
        interactors_response = requests.get(base_url, params=params)
        interactors_response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Error: Unable to get interaction data ({e})")
        return None

    # Data parsing logic outside try-except to ensure it's executed if no exception occurs
    data = interactors_response.text.strip().splitlines()
    interactors = []

    for line in data[1:]:
        columns = line.split('\t')
        target_protein = columns[2]
        target_id = columns[0]
        interactor = columns[3]
        interactor_id = columns[1]
        score = float(columns[5])
        escore = float(columns[10])
        dscore = float(columns[11])
        tscore = float(columns[12])

        if score >= threshold:
            interactors.append((target_protein, target_id, interactor, interactor_id, score, escore, dscore, tscore))

    return interactors


def main():
    parser = argparse.ArgumentParser(
        description="Retrieval of interaction data from the STRING database for a given identifier."
    )
    parser.add_argument(
        "identifier", 
        type=str,
        help="STRING identifier or gene name to retrieve interactors for."
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

    # Check if identifier is in the dataframe
    if args.identifier in data['preferredName'].values:
        if data.shape[0] > 1:
            print(f"Warning: Multiple STRING IDs found for {args.identifier}. Using the first one: {data.iloc[0, 1]}")
        string_id = data.iloc[0]['stringId']
    else:
        print(f"Error: No STRING identifier found for {args.identifier} in the results.")
        return
    
    if string_id:
        interactors = get_interactors(string_id, args.threshold, args.network)
        if interactors:
            output_file = f"{args.identifier}_string_interactors.tsv"
            
            with open(output_file, 'w', newline='') as tsvfile:
                tsv_writer = csv.writer(tsvfile, delimiter='\t')
                
                tsv_writer.writerow(['Target_protein', 'Target_id', 'Interactor', 'Interactor_id', 'String_score', 'Experimental_score', 'Database_score', 'Textmining_score'])
                
                for target_protein, target_id, interactor, interactor_id, score, escore, dscore, tscore in interactors:
                    tsv_writer.writerow([target_protein, target_id, interactor, interactor_id, score, escore, dscore, tscore])
            
            print(f"Results saved to {output_file}")
        else:
            print("No interactors found.")
    else:
        print("No STRING ID found.")


if __name__ == "__main__":
    main()
