#!/usr/bin/env python3

import requests
import argparse
import csv


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "identifier", 
        type=str,
    )
    
    parser.add_argument(
        "-t", 
        "--threshold",
        type=int,
        default=700,
    )

    parser.add_argument(
        "-n",
        "--network",
        type=str,
        default= "physical",
        choices=["functional", "physical"]
    )

    return parser.parse_args()


def get_string_id(identifier):
    base_url = "https://string-db.org/api/tsv/get_string_ids"
    params = {
        'identifier': identifier,  
        'species': 9606,         
        'limit': 1,
        'caller_identity': "MAVISp_web_app"     
    }

    try:
        string_response = requests.get(base_url, params=params)
        string_response.raise_for_status()
        data = string_response.text.strip().splitlines()  
        if len(data) > 1:
            columns = data[1].split('\t')  
            return columns[1]  
        else:
            print(f"Error: No STRING identifier found for {identifier}")
            return None

    except requests.exceptions.RequestException as e:
        print(f"Error: Unable to get data ({e})")
        return None


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
            

            if score >= threshold / 1000:
                interactors.append((target_protein, target_id, interactor, interactor_id, score, escore, dscore, tscore))

        return interactors

    except requests.exceptions.RequestException as e:
        print(f"Error: Unable to get interaction data ({e})")
        return None

def main():
    args = parse_arguments()
    string_id = get_string_id(args.identifier)
    
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