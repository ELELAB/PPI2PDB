# -*- coding: utf-8 -*-
"""
Created on Thu May 12 09:14:52 2022

@author: Matteo
"""

import sys
import argparse
from decimal import Decimal
import pandas as pd
import pypdb
import re
import requests
import warnings
import csv

def main(argv):

    warnings.filterwarnings("ignore")
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--i', help='mentha database file')
    parser.add_argument('-t','--t', help='File with target uniprots')
    parser.add_argument('-s','--s', type=Decimal, help='Cutoff score')
    parser.add_argument('-o','--o', nargs='?',const='dataframe.csv',default='dataframe.csv',help='Output name') 
    parser.add_argument('-f', '--filter', action='store_true')
    
    args = parser.parse_args()
    
    filterSameProteinInteraction = False
    
    if args.filter:
        filterSameProteinInteraction = True
        
    args = parser.parse_args()
    
    #read data with pandas
    data = pd.read_csv(args.i, sep=';', converters={'Score':Decimal})
    
    #filtering for taxon.A = 9606 AND taxon.B = 9606 AND score >= cutoff (args.s)
    data = data[(data['Taxon A'] == 9606) & (data['Taxon B'] == 9606) & (data['Score'] >= args.s)]
    
    dataframeOut = pd.DataFrame(columns= ['target uniprot id', 'target uniprot gene',                                          #2 -> from csv
                                          'interactor uniprot id', 'interactor uniprot gene',                                  #2 -> from csv
                                          'mentha score',                                                                       #1 -> from csv
                                          'PDB id',                                                                            #1 -> from pypdb lib
                                          'fusion',                                                                            #1 -> from summary request
                                          'target chain id','target starting residue','target ending residue',                 #3 -> from mappings request
                                          'interactor chain id','interactor starting residue','interactor ending residue',     #3 -> from mappings request
                                          'other interactors',                                                                 #1 -> from mappings request
                                          'method',                                                                            #1 -> from summary request
                                          'resolution',                                                                        #1 -> from experiment request
                                          'dna chains','num ligands'])                                                         #2 -> from summary request
                                                                                                                               #-----------------------------
                                                                                                                               #18 columns total                                                                                                         
    
    #open uniprot target file and get lines
    with open(args.t,'r') as uniprotTargets:
        for uniprot in uniprotTargets:
            uniprot = uniprot.rstrip()
            
            #get data where protein A or protein b matches uniprot selected
            uniprotData = data[(data['Protein A'] == uniprot) | (data['Protein B'] == uniprot)]
            
            uniprotData = uniprotData.reset_index()  # make sure indexes pair with number of rows
            
            targetQueryResult = pypdb.Query(uniprot).search(num_attempts=10,sleep_time=0.9) 
            print('Target {}'.format(uniprot))
            for index, row in uniprotData.iterrows():
                targetProtein = ''
                interactorProtein = ''
                targetGene = ''
                interactorGene = ''
                score = Decimal('0')
                #list that will later added to dataframe out
                outRow = []
                
                if row['Protein A'] == uniprot:
                    targetProtein = row['Protein A']
                    interactorProtein = row['Protein B']
                    targetGene = row['Gene A']
                    interactorGene = row['Gene B']
                else:
                    targetProtein = row['Protein B']
                    interactorProtein = row['Protein A']
                    targetGene = row['Gene B']
                    interactorGene = row['Gene A']
                
                #filter protein interaction with self
                if  filterSameProteinInteraction and row['Protein A'] == row['Protein B']:
                    print('skipped protein interaction with self \n \t target {} interactor'.format(targetProtein,interactorProtein))
                    continue
                
                score = row['Score']
                
                #setup first 5 of outRow
                outRow.extend([targetProtein,targetGene,interactorProtein,interactorGene,score])
                
                #sending pypdb requests
                #targetQueryResult = pypdb.Query(targetProtein).search(num_attempts=10,sleep_time=0.9)   
                interactorQueryResult = pypdb.Query(interactorProtein).search(num_attempts=10,sleep_time=0.9)
                
                #check if something went wrong in pypdb -> set na and go next
                if interactorQueryResult is None or targetQueryResult is None:
                    #set output row to na (13 cause we had 5 set and 13 missing positions)
                    outRow.extend(['na']*13)
                    print('\t PYPDB nonetype returned {}                 '.format(interactorProtein),end='\r')
                    #append row to dataframe Out                  
                    #rowSerie = pd.Series(outRow, index = dataframeOut.columns)
                    #dataframeOut = dataframeOut.append(rowSerie, ignore_index=True)
                    dataframeOut.loc[len(dataframeOut)] = outRow
                else:
                    #get common pdbs to both proteins
                    commonPdbs = set(targetQueryResult).intersection(set(interactorQueryResult))
                
                    #if intersection is not empty
                    if commonPdbs != set():
                        print('\t protein interactor {} share pdbs -> {}'.format(interactorProtein,commonPdbs))
                        #intersection not empty -> run requests to get other columns
                        #do requests
                        # !! can be multiple pdbs !!
                        for pdb in commonPdbs:
                            fused = ''
                            dna = ''
                            ligands = ''
                            method = ''
                            targetChainIds = ''
                            targetStart = ''
                            targetEnd = '' 
                            interactorChainIds = '' 
                            interactorStart = '' 
                            interactorEnd = '' 
                            otherInteractors = ''
                            resolution = ''
                            
                            #do requests on pdb chosen
                            #request summary -> from summary we get fusion, method,dna chains, num ligands
                            fused, dna, ligands, method = get_summary(pdb)
                            #request mappings -> from mappings we get chain infos (id, start, stop) and other interactors
                            targetChainIds, targetStart, targetEnd, interactorChainIds, interactorStart, interactorEnd, otherInteractors = get_mappings_data(pdb,targetProtein,interactorProtein)
                            #request experiment -> from experiment we get resolution
                            resolution = get_experiment(pdb)
                            
                            #add data to output row
                            outRow.extend([pdb,fused,targetChainIds,targetStart,targetEnd,interactorChainIds,interactorStart,interactorEnd,otherInteractors,method,resolution,dna,ligands])
                            
                            #add out row to dataframe
                            #rowSerie = pd.Series(outRow, index = dataframeOut.columns)
                            #dataframeOut = dataframeOut.append(rowSerie, ignore_index=True)
                            dataframeOut.loc[len(dataframeOut)] = outRow
                            
                            #reset outrow to first five values -> first five are fixed until we don't change target - interactor pair
                            #other values change basing on the pdb selected
                            outRow = outRow[:5]
                    else:
                        print('\t interactor {} ->  NO COMMON PDBS'.format(interactorProtein),end='\r')
                        #intersection empty -> set na and go on
                        #empty  pdb list, no requests set na
                        outRow.extend(['na']*13)
                        #rowSerie = pd.Series(outRow, index = dataframeOut.columns)
                        #dataframeOut = dataframeOut.append(rowSerie, ignore_index=True)
                        dataframeOut.loc[len(dataframeOut)] = outRow
                        #print(dataframeOut)
                        

    
    print('\n')
    dataframeOut.to_csv(args.o, index=False, quoting=csv.QUOTE_NONE)
                        
            
def get_experiment(pdb):
    """
    This function retrieves PDB > experiment
        
    :param pdb_id: String,
    :return: resolution: String
    """ 
    
    url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/'
    data = make_request(url, "get", pdb)
    resolution = ''
    
    if data is None:
        resolution = 'none'
        #print("########")
        #print("EXPERIMENT CALL ERROR -> no data")
        #print("########")
    else:
        if 'resolution' in data[pdb.lower()][0].keys():
            resolution = Decimal(str(data[pdb.lower()][0]['resolution']))
        else:
            #print('resolution not in data')
            resolution = 'na'
            
    return resolution
    
def get_summary(pdb):
    """
    This function retrieves PDB > summary
        
    :param pdb_id: String,
    :return: fused: String ('yes'/'')
             dna: String 
             ligands: String
             method: String
    """ 
      
    url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/'
    data = make_request(url, "get", pdb)
    
    dna = ''
    ligands = ''
    method = ''
    title = ''
    fused = ''
    
    #data None -> means that we have no data from the request
    #no data -> no title
    #no title -> fused = ''
    if data is None:
        title = 'NO TITLE'
        dna = 'none'
        ligands = 'none'
        method = 'none'
        fused = 'none'
        #print("########")
        #print("PDB ", pdb, " has no title :(")
        #print("########")
    else:
        #get INFOs
        title = data[pdb.lower()][0]['title']       
        dna = data[pdb.lower()][0]['number_of_entities']['dna']
        ligands = data[pdb.lower()][0]['number_of_entities']['ligand']
        method = data[pdb.lower()][0]['experimental_method'][0]
        
    
        if 'fused' in title or 'fusion' in title:
            fused = 'yes'
        else:
            fused = 'na'
    
    return fused,dna,ligands,method

def get_mappings_data(pdb,targetProtein,interactorProtein):
    """
    This function will GET the mappings data from
    the PDBe API using the make_request() function
    
    :param pdb_id: String
    :return: targetChainIds: String
             targetStart: String
             targetEnd: String
             interactorChainIds: String
             interactorStart: String
             interactorEnd: String
             otherInteractors: String 
    """
    
    base_url = "https://www.ebi.ac.uk/pdbe/"
    api_base = base_url + "api/"
    uniprot_mapping_url = api_base + 'mappings/uniprot/'
    # Check if the provided PDB id is valid
    # There is no point in making an API call
    # with bad PDB ids
    if not re.match("[0-9][A-Za-z][A-Za-z0-9]{2}", pdb):
        #print("Invalid PDB id")
        return 'none','none','none','none','none','none','none'
    
    # GET the mappings data
    mappings_data = make_request(uniprot_mapping_url, "get", pdb)
    
    targetChainIds = []
    interactorChainIds = []
    targetStart = []
    targetEnd = []
    interactorStart = []
    interactorEnd = []
    otherInteractors = []
    
    # Check if there is data
    if not mappings_data:
        #print("NA")
        mappings_data = ['NA']
        return 'none','none','none','none','none','none','none'
    else:
        #extract chains data
        #dict that contains uniprots
        uniprot = mappings_data[pdb.lower()]['UniProt']
        for uID in uniprot.keys():
            if uID != targetProtein and uID != interactorProtein:
                otherInteractors.append(uID)
            else:
                if uID == targetProtein:
                    uniprotTarget = uniprot[targetProtein]['mappings']
                    for mapping in uniprotTarget:
                        targetChainIds.append(mapping['chain_id'])
                        targetStart.append(mapping['unp_start'])
                        targetEnd.append(mapping['unp_end'])
                if uID == interactorProtein:
                    uniprotInteractor = uniprot[interactorProtein]['mappings']
                    for mapping in uniprotInteractor:
                        interactorChainIds.append(mapping['chain_id'])
                        interactorStart.append(mapping['unp_start'])
                        interactorEnd.append(mapping['unp_end'])
    
    #convert lists to strings
    targetChainIds = ';'.join([str(x) for x in targetChainIds])
    targetStart = ';'.join([str(x) for x in targetStart])
    targetEnd = ';'.join([str(x) for x in targetEnd])
    
    interactorChainIds = ';'.join([str(x) for x in interactorChainIds])
    interactorStart = ';'.join([str(x) for x in interactorStart])
    interactorEnd = ';'.join([str(x) for x in interactorEnd])
    
    otherInteractors = ';'.join([str(x) for x in otherInteractors])
    if otherInteractors == '':
        otherInteractors = 'na'
        
    return targetChainIds,targetStart,targetEnd,interactorChainIds,interactorStart,interactorEnd,otherInteractors     
  
def make_request(url, mode, pdb_id):
    """
    This function can make GET and POST requests to
    the PDBe API
    
    :param url: String,
    :param mode: String,
    :param pdb_id: String
    :return: JSON or None
    """
    if mode == "get":
        response = requests.get(url=url+pdb_id)
    elif mode == "post":
        response = requests.post(url, data=pdb_id)

    if response.status_code == 200:
        return response.json()
    else:
        print("NA from ", url," for pdb ", pdb_id)
    
    return None 
  
if __name__ == "__main__":
    main(sys.argv[1:])
