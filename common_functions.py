#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import requests
import tarfile
import io
import json
import re
import obonet
import networkx as nx
import numpy as np
import time
import ftplib
import os
import zipfile
from git import Repo
import shutil

def download_ftp_directory(ftp_url, download_dir, exclude_files=None):
    """
    This function downloads all files from the specified FTP directory and saves them to the specified destination folder,
    excluding subdirectories and the specified files.

    Parameters:
    ftp_url (str): The FTP URL of the directory containing the files to be downloaded.
    download_dir (str): The folder where the downloaded files will be saved.
    exclude_files (list): A list of filenames to be excluded from downloading.
    """
    if exclude_files is None:
        exclude_files = []

    if not os.path.exists(download_dir):
        os.makedirs(download_dir)
    
    # Parse the FTP URL
    ftp_host = ftp_url.split('/')[2]
    ftp_dir = '/'.join(ftp_url.split('/')[3:])
    
    # Connect to the FTP server
    ftp = ftplib.FTP(ftp_host)
    ftp.login()
    ftp.cwd(ftp_dir)
    
    # List files in the directory
    items = ftp.mlsd()
    
    # Exclude subdirectories and the specified files
    files = []
    for item, facts in items:
        if item in exclude_files:
            continue
        if facts['type'] == 'file':  # Check if item is a file
            files.append(item)
    
    # Download each file
    for file in files:
        local_file_path = os.path.join(download_dir, file)
        with open(local_file_path, 'wb') as local_file:
            ftp.retrbinary('RETR ' + file, local_file.write)
        print(f"Downloaded: {file}")

    # Close the FTP connection
    ftp.quit()

def clone_repo(bitbucket_url, clone_dir):
    """
    This function clones a Bitbucket repository to the specified directory using GitPython.

    Parameters:
    bitbucket_url (str): The URL of the Bitbucket repository to be cloned.
    clone_dir (str): The directory where the repository will be cloned.
    """
    if os.path.exists(clone_dir):
        # Remove the existing directory and its contents
        shutil.rmtree(clone_dir)
    
    os.makedirs(clone_dir)
    
    try:
        # Clone the repository
        Repo.clone_from(bitbucket_url, clone_dir)
        print(f"Cloned repository to: {clone_dir}")
    except Exception as e:
        print(f"Failed to clone repository: {e}")

def download_and_extract_file(url, download_dir):
    """
    This function downloads a tar.bz2 file from a given URL and extracts its contents to the specified destination folder
    without saving the tar.bz2 file to disk.

    Parameters:
    url (str): The URL of the tar.bz2 file to be downloaded.
    download_dir (str): The folder where the extracted files will be saved.

    Returns:
    str: The path to the extracted contents.
    """
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)
    
    response = requests.get(url)
    if response.status_code == 200:
        # Create a TarFile object from the downloaded content
        file_like_object = io.BytesIO(response.content)
        with tarfile.open(fileobj=file_like_object, mode='r:bz2') as tar:
            tar.extractall(path=download_dir)
        print(f"Extracted files to: {download_dir}")
        return download_dir
    else:
        print(f"Failed to download {url}")
        return None

def card_data_to_content(url, file):
    """
    This function downloads a tar.bz2 file from the provided URL, extracts a specific file from the archive, 
    reads its content, and returns it as a string.

    Parameters:
    url (str): The URL of the tar.bz2 file to be downloaded.
    file (str): The relative path of the specific file to be extracted from the tar.bz2 archive.

    Returns:
    str: The content of the extracted file as a UTF-8 encoded string, or None if the file is not found 
         or there is an error during the process.
    """
    response = requests.get(url)
    if response.status_code == 200:
        # Create a TarFile object from the downloaded content
        file_like_object = io.BytesIO(response.content)
        with tarfile.open(fileobj=file_like_object, mode='r:bz2') as tar:
            # List files in the tar.bz2 archive
            files_in_tar = tar.getnames()
            # print("Files in the tar.bz2:", files_in_tar)

            file = './'+file

            if file in files_in_tar:
                extracted_file = tar.extractfile(file)
                if extracted_file is not None:
                    # Read the content of the extracted file
                    content = extracted_file.read().decode('utf-8')
                    print(f"Data extracted and readed from: {file}")
                    return content
            else:
                print(f"The file {file} was not found in the tar.bz2.")
    else:
        print(f"Error downloading the file: {response.status_code}")

def parse_obo_file(obo_path):
    graph = obonet.read_obo(obo_path)
    return graph

def extract_subgraph(graph, node_id):
    """
    Extracts a subgraph starting from a specific node and all its descendants.
    
    :param graph: The original graph (networkx.DiGraph)
    :param node_id: The identifier of the root node for the subgraph
    :return: The corresponding subgraph (networkx.DiGraph)
    """
    # Ensure the node is in the graph
    if node_id not in graph:
        raise ValueError(f"The node {node_id} is not in the graph.")
    
    # Get all descendants of the node
    descendants = nx.descendants(graph, node_id)
    descendants.add(node_id)  # Include the root node
    
    # Create the subgraph
    subgraph = graph.subgraph(descendants).copy()
    
    return subgraph

def card_graph(obo_file, json_file, categories_file, map_file, acc, colors):
    """
    This function generates a graph visualization for a specific antimicrobial resistance (AMR) gene
    using data from OBO, JSON files. It filters the relevant data, constructs a subgraph,
    and adds nodes and edges with appropriate labels and attributes.

    Parameters:
    obo_file (str): Path to the OBO file containing ontology information.
    json_file (str): Path to the JSON file containing CARD data.
    map_file (str): Path to the TSV file mapping UniProtKB accessions to ARO accessions.
    acc (str): The UniProtKB accession number to search for.
    color (str): The color to be used for nodes in the graph.

    Returns:
    graph (Graph): The constructed graph object with nodes and edges.
    aro (str): The ARO accession number associated with the given UniProtKB accession.
    """
    
    ### get aro from map file
    map_df = pd.read_csv(map_file, sep='\t')
    map_aro_df = map_df.loc[map_df['UniProtKB_acc'].str.contains(acc, na=False)]
    if map_aro_df.empty: # check accession 
        return None, None
    else:
        aro = map_aro_df['ARO Accession'].iloc[0]
        print(aro)
    
    obo_graph = parse_obo_file(obo_path=obo_file)
    graph = extract_subgraph(graph=obo_graph, node_id=aro)
    
    ## categories files
    cat_df = pd.read_csv(categories_file, sep='\t')

    for node, data in graph.nodes(data=True):
        # get short name
        syn = data.get('synonym', None)
        label = None
        if syn:
            for i in syn:
                if 'CARD_Short_Name' in i:
                    label = re.match(r'"(.+)"', i).group(1)

        if label:
            data['label'] = label
        else:
            data['label'] = data.get('name')
        
        data['title'] = f'{node}: {data.get('name')}; {data.get('def')}'
        data['group'] = 'card'
        data['color'] = colors['card']

        if node in cat_df['ARO Accession'].values:
            cat = cat_df.loc[cat_df['ARO Accession'] == node, 'ARO Category'].iloc[0]
            data['color'] = colors[cat]
            data['title'] = f'{node}: {cat}: {data.get('name')}; {data.get('def')}'


    for source, target, attr in graph.edges:
        graph.edges[source, target, attr]['label'] = attr

    ## json file parse to get SNPs
    with open(json_file, 'r') as file:
        json_data = json.load(file)

    for key, value in json_data.items():
        if isinstance(value, dict):
            if 'ARO_accession' in value.keys(): ### take into account protein overexpresion model as well
                if str(value.get('ARO_accession', False)) in aro and value['model_type'] == 'protein variant model':
                    graph.add_node('SNPs', **{
                        'name': 'SNPs',
                        'label': 'SNPs' ,
                        'title': value['model_param']['snp']['param_description'],
                        'group': 'card',
                        'color': colors['card']
                    })
                    
                    graph.add_edge(aro, 'SNPs', label=value['model_type'], title=value['model_description'])

                    for snp in value['model_param']['snp']['param_value'].values():
                        graph.add_node(snp, **{
                        'name': snp,
                        'label': snp ,
                        'title': snp,
                        'group': 'card',
                        'color': colors['card']
                        })
                        graph.add_edge('SNPs', snp, label=value['model_param']['snp']['param_type'])

    ## remove some nodes for a clearer graph
    graph.remove_node('ARO:1000001') # process or component of antibiotic biology or chemistry
    graph.remove_node('ARO:1000002') # mechanism of antibiotic resistance
    graph.remove_node('ARO:1000003') # antibiotic molecule
    graph.remove_node('ARO:3000000') # determinant of antibiotic resistance

    return graph, aro

def amrfinderplus_graph(database_dir, map_file, acc, color):
    """
    This function generates a graph visualization for a specific antimicrobial resistance (AMR) gene
    using data from database files and a map file. It constructs a subgraph based on the AMR gene's
    hierarchy and adds nodes and edges with appropriate labels and attributes.

    Parameters:
    database_dir (str): Path to the directory containing database files.
    map_file (str): Path to the TSV file mapping UniProtKB accessions to RefSeq protein accessions.
    acc (str): The UniProtKB accession number to search for.
    color (str): The color to be used for nodes in the graph.

    Returns:
    graph (Graph): The constructed graph object with nodes and edges.
    node_name (str): The name of the root node of the graph.
    """
    
    ### create dataframe from files
    cat = pd.read_csv(database_dir + '/ReferenceGeneCatalog.txt', sep='\t')
    hier = pd.read_csv(database_dir + '/ReferenceGeneHierarchy.txt', sep='\t')
    map_df = pd.read_csv(map_file, sep='\t')

    ### get RefSeq Protein accession from map file
    filt_map = map_df.loc[map_df['UniProtKB_acc'].str.contains(acc, na=False)]
    if filt_map.empty:
        return None, None
    else:
        ext_acc = filt_map['refseq_protein_accession'].iloc[0]
        print(ext_acc)
    
    ### get data of refseq protein accession
    tsv_acc = cat.loc[cat['refseq_protein_accession'] == ext_acc]
    tsv_acc = tsv_acc.replace({np.nan: None})

    ### define graph
    graph = nx.MultiDiGraph()

    nodes = ['gene_family', 'allele', 'synonyms', 'subclass', 'class', 'subtype', 'type', 'scope', 'SNPs']
    prefix = 'afp_'
    product_name = tsv_acc.iloc[0]['product_name']
    pubmed_reference = 'PMID:' + str(tsv_acc.iloc[0]['pubmed_reference'])
    subtype = tsv_acc.iloc[0]['subtype']
    allele = tsv_acc.iloc[0]['allele']
    group = 'amrfinderplus'

    for node in nodes:
        
        node_id = prefix + node
        title = node + ': ' + product_name + ' [' + pubmed_reference + ']'
        if node != 'SNPs':
            label = tsv_acc.iloc[0][node]
        else:
            label = ''

        data_dict = {'name': node, 'label': label, 'title': node, 'group': group, 'color': color}

        if node == 'SNPs':
            if subtype == 'POINT':
                snps = tsv_acc['allele'].str.split('_', expand=True).iloc[:,-1].to_list()
                data_dict['label'] = '; '.join(snps)
                graph.add_node(node_id, **data_dict)
        
        elif node == 'gene_family':
            if allele == None:
                data_dict['title'] = title
                graph.add_node(node_id, **data_dict)
            else:
                graph.add_node(node_id, **data_dict)
        
        elif node == 'synonyms':
            synonym = tsv_acc.iloc[0][node]
            data_dict['label'] = synonym
            if synonym != None:
                graph.add_node(node_id, **data_dict)
        
        elif node == 'allele':
            if subtype != 'POINT':
                if allele != None:
                    data_dict['label'] = allele
                    graph.add_node(node_id, **data_dict)
        
        else:
            graph.add_node(node_id, **data_dict)

    ### get hierarchical data
    hier_array = []
    if allele != None and subtype != 'POINT':
        node_name = tsv_acc.iloc[0]['allele']
    else:
        node_name = tsv_acc.iloc[0]['gene_family']

    hier_num = 0
    while True:
        node_id = 'hier' + str(hier_num)
        filt_hier = hier.loc[hier['node_id'] == node_name]
        
        if filt_hier.empty:
            print('No hierarchical data')
            break
        node_name = hier.loc[hier['node_id'] == node_name].iloc[0]['parent_node_id']
        node_data = hier.loc[hier['node_id'] == node_name].iloc[0]['name']
        
        if node_data != node_data:
            node_data = ''
        hier_array.append((node_id, node_name, node_data))
        hier_num += 1
        
        if node_name == 'ALL':
            break

    ### define hierarchical nodes
    for item in hier_array:
        graph.add_node(item[0], **{
                'name': item[1],
                'label': item[1],
                'title': item[2],
                'group': group,
                'color': color
            })
    ### define hierarchical edges:
    for i in range(0, len(hier_array)-1):
        graph.add_edge(hier_array[i][0], hier_array[i+1][0], label='is_a', title='is_a')

    ### define edges
    graph.add_edge(prefix + 'subtype', prefix + 'type', label='part_of', title='part_of')
    graph.add_edge(prefix + 'type', prefix + 'scope', label='part_of', title='part_of')
    graph.add_edge(prefix + 'subclass', prefix + 'class', label='is_a', title='is_a')
    
    if len(hier_array) > 0:
        graph.add_edge(prefix + 'gene_family', 'hier0', label='is_a', title='is_a')

    if allele == None:
        graph.add_edge(prefix + 'gene_family', prefix + 'subtype', label='part_of', title='part_of')
        graph.add_edge(prefix + 'gene_family', prefix + 'subclass', label='confers_resistace_to_subclass', title='confers_resistace_to_subclass')
        if synonym != None:
            graph.add_edge(prefix + 'gene_family', prefix + 'synonyms', label='synonym', title='synonym')

        return graph, prefix + 'gene_family'

    elif subtype == 'POINT':
        graph.add_edge(prefix + 'gene_family', prefix + 'SNPs', label='SNPs', title='AMRFinderPlus SNPs')
        graph.add_edge(prefix + 'gene_family', prefix + 'subtype', label='part_of', title='part_of')
        graph.add_edge(prefix + 'gene_family', prefix + 'subclass', label='confers_resistace_to_subclass')
        if synonym != None:
            graph.add_edge(prefix + 'gene_family', prefix + 'synonyms', label='synonym', title='synonym')

        return graph, prefix + 'gene_family'

    elif allele != None:
        graph.add_edge(prefix + 'allele', prefix + 'gene_family', label='is_a', title='is_a')
        graph.add_edge(prefix + 'allele', prefix + 'subtype', label='part_of', title='part_of')
        graph.add_edge(prefix + 'allele', prefix + 'subclass', label='confers_resistace_to_subclass', title='confers_resistace_to_subclass')
        if synonym != None:
            graph.add_edge(prefix + 'allele', prefix + 'synonyms', label='synonym', title='synonym')

        return graph, prefix + 'allele'

def resfinder_graph(phenotype, map_file, acc, color):
    """
    This function generates a graph visualization for a specific antimicrobial resistance (AMR) gene from ResFinder database
    using data from phenotypes.txt and map files. It constructs a subgraph based on the AMR gene's relationships
    and adds nodes and edges with appropriate labels and attributes.

    Parameters:
    phenotype (str): Path to the directory containing phenotype data.
    map_file (str): Path to the TSV file mapping UniProtKB accessions to ResFinder accessions.
    acc (str): The UniProtKB accession number to search for.
    color (str): The color to be used for nodes in the graph.

    Returns:
    graph (Graph): The constructed graph object with nodes and edges.
    gene_accession (str): The gene accession number associated with the given UniProtKB accession.
    """

    phen = pd.read_csv(phenotype + '/phenotypes.txt', sep='\t')
    prot_ids_col = 'resfinder_accession'
    phen[prot_ids_col] = phen['Gene_accession no.'].str.split('_').apply(lambda x: x[-1])
    for i in range(0, phen.shape[0]):
            prefix = phen.loc[i,'Gene_accession no.'].split('_')[-2]
            if prefix == 'NG' or prefix == 'NC':
                phen.loc[i,prot_ids_col] = f'{prefix}_{phen.loc[i,'Gene_accession no.'].split('_')[-1]}'
    
    map_df = pd.read_csv(map_file, sep='\t')

    ### get resfinder accession from map file
    filt_df = map_df.loc[map_df['UniProtKB_acc'].str.contains(acc, na=False)]
    if filt_df.empty:
        return None, None
    else:
        resfinder_acc = filt_df[prot_ids_col].iloc[0]
        print(resfinder_acc)

    ### get data of protein accession
    tsv_acc = phen.loc[phen[prot_ids_col] == resfinder_acc]
    tsv_acc = tsv_acc.replace({np.nan: None})

    gene_accession = tsv_acc.iloc[0]['Gene_accession no.']
    res_class = tsv_acc.iloc[0]['Class']
    phenotype = tsv_acc.iloc[0]['Phenotype']
    pmids = 'PMIDS: ' + tsv_acc.iloc[0]['PMID']
    mechanism = tsv_acc.iloc[0]['Mechanism of resistance']
    notes = tsv_acc.iloc[0]['Notes']
    required_gene = tsv_acc.iloc[0]['Required_gene']

    graph = nx.MultiDiGraph()
    graph.add_node(gene_accession, **{
            'name': gene_accession,
            'label': gene_accession,
            'title': pmids,
            'group': 'resfinder',
            'color': color
        })
    graph.add_node(res_class, **{
            'name': res_class,
            'label': res_class,
            'title': 'Class',
            'group': 'resfinder',
            'color': color
        })
    graph.add_node(phenotype, **{
            'name': phenotype,
            'label': phenotype,
            'title': 'phenotype',
            'group': 'resfinder',
            'color': color
        })
    graph.add_node(mechanism, **{
            'name': mechanism,
            'label': mechanism,
            'title': 'mechanism',
            'group': 'resfinder',
            'color': color
        })

    if required_gene != None:

        graph.add_node(required_gene, **{
            'name': required_gene,
            'label': required_gene,
            'title': 'required_gene',
            'group': 'resfinder',
            'color': color
            })
        graph.add_edge(gene_accession, required_gene, label='required_of')
        graph.add_edge(required_gene, phenotype, label='resistace_to')

    if notes != None:
        graph.add_node(notes, **{
            'name': notes,
            'label': notes,
            'title': 'notes',
            'group': 'resfinder',
            'color': color
        })
        graph.add_edge(gene_accession, notes, label='is_a')

    graph.add_edge(gene_accession, res_class, label='is_a')
    graph.add_edge(gene_accession, mechanism, label='participates_in')
    graph.add_edge(mechanism, phenotype, label='confers_resistace_to')

    return graph, gene_accession

def simplify_graph(graph):
    """
    ¡¡¡ On develop !!!
    This function simplifies a given graph by merging nodes with identical predecessors and successors. 
    It creates a new node to represent the merged nodes and updates the edges accordingly.

    Parameters:
    graph (Graph): The original graph to be simplified.

    Returns:
    simplified_graph (Graph): The simplified graph with merged nodes.
    """
    # Create a copy of the graph to make modifications
    simplified_graph = graph.copy()
    
    # Find nodes with the same parent and successor
    nodes_to_merge = {}
    for node in graph.nodes():
        
        predecessors = tuple(sorted(graph.predecessors(node)))
        successors = tuple(sorted(graph.successors(node)))

        if len(predecessors) == 0 or len(successors) == 0:
            continue

        key = (predecessors, successors)
        # print(key)
        if key not in nodes_to_merge:
            nodes_to_merge[key] = []
        nodes_to_merge[key].append(node)

    # Merge nodes
    for (predecessors, successors), nodes in nodes_to_merge.items():
        if len(nodes) > 1:
            # Create a new node that represents all merged nodes
            new_node = '; '.join(nodes)
            simplified_graph.add_node(new_node)
            
            # Connect predecessors and successors to the new node
            for pred in predecessors:
                simplified_graph.add_edge(pred, new_node, label='is_a', title='is_a')
            for succ in successors:
                simplified_graph.add_edge(new_node, succ, label='is_a', title='is_a')
            
            # Remove old nodes
            for node in nodes:
                simplified_graph.remove_node(node)
    
    return simplified_graph

def insert_newline_every_n_words(text, n):
    words = text.split()
    return '\n'.join([' '.join(words[i:i + n]) for i in range(0, len(words), n)])

def insert_newline_every_n_chars(text, n):
    result = []
    current_line = []

    for word in text.split():
        current_line.append(word)
        if len(' '.join(current_line)) >= n:
            result.append(' '.join(current_line))
            current_line = []

    # Append the last line if it's not empty
    if current_line:
        result.append(' '.join(current_line))

    return '\n'.join(result)

def run_idmapping(accesion, from_db, to_db):
    """
    This function initiates an ID mapping job using the UniProt API, monitors the job status, 
    and returns the job ID once the job is finished.

    Parameters:
    accesion (str): The accession ID to be mapped.
    from_db (str): The source database of the accession ID.
    to_db (str): The target database for the ID mapping.

    Returns:
    jobid (str): The job ID of the ID mapping job.
    """

    # Define the URL of the API of IDmapping
    url = 'https://rest.uniprot.org/idmapping/run'
    # Define parameter
    data = {
        'ids': accesion,
        'from': from_db,
        'to': to_db
    }
    
    # do POST
    response = requests.post(url, data=data)
    if response.status_code == 200:
        jobid = response.json()['jobId']
        print(f'{jobid} POST')
        
    else:
        print(f'error running ID mapping {response.text}')

    status_url = f'https://rest.uniprot.org/idmapping/status/{jobid}'

    while True:
        status_response = requests.get(status_url)
        if status_response.status_code == 200:
            status = status_response.json()
            if 'results' in status.keys():
                print(f'{jobid} FINISHED')
                break
            elif status['jobStatus'] == 'RUNNING':
                print(f"{jobid} RUNNING")
                time.sleep(5)
            elif status['jobStatus'] == 'FAILED':
                print(f"{jobid} FAILED")
                break
            else:
                print(status)
        else:
            print("Error, status code: ", status_response.status_code)
            print("Error: ", status_response.text)
            break
        
    return response.json()['jobId']

def get_idmapping_results(jobid):
    """
    This function retrieves the results of an ID mapping job from the UniProt API 
    and returns them as a DataFrame.

    Parameters:
    jobid (str): The job ID of the completed ID mapping job.

    Returns:
    df (DataFrame): A DataFrame containing the mapping results, with 'from' IDs grouped 
                    and concatenated 'to' IDs.
    """

    # Define URL
    print(jobid, 'GETTING')
    results_url = f'https://rest.uniprot.org/idmapping/stream/{jobid}'
    
    # get results
    results_response = requests.get(results_url)
    if results_response.status_code == 200:
        results = results_response.json()
        # print(results['results'])
        df = pd.DataFrame(results['results'])
        df = df.groupby('from').agg({'to': lambda x: ';'.join(x)}).reset_index()
        return df
    else:
        print("Error: ", results_response.status_code)
        print("Error: ", results_response.text)
        return None

def get_idmapping_suggestedIds(jobid):
    """
    This function retrieves the suggested IDs from an ID mapping job using the UniProt API
    and returns them as a DataFrame.

    Parameters:
    jobid (str): The job ID of the completed ID mapping job.

    Returns:
    df (DataFrame): A DataFrame containing the suggested IDs, with 'from' IDs grouped 
                    and concatenated 'to' IDs. Returns None if no suggested IDs are found.
    """

    # Define URL
    print(jobid, 'GETTING NO STREAM')
    results_url = f'https://rest.uniprot.org/idmapping/results/{jobid}'
    
    # get results
    results_response = requests.get(results_url)
    if results_response.status_code == 200:
        results = results_response.json()
        if 'suggestedIds' not in results.keys():
            print(jobid, 'NO SUGGESTED IDS')
            return None
        else:
            print(results['suggestedIds'])
            df = pd.DataFrame(results['suggestedIds'])
            df = df.groupby('from').agg({'to': lambda x: ';'.join(x)}).reset_index()
            return df
    else:
        print("Error: ", results_response.status_code)
        print("Error: ", results_response.text)
        return None

def mapping_card_protein_accession(card_data_dir):
    """
    This function maps CARD protein accessions to UniProtKB accessions using the UniProt ID mapping API. 
    It reads the CARD data from a TSV file, performs the mapping, and returns a DataFrame with the results.

    Parameters:
    card_data_dir (str): Path to the directory containing CARD data.

    Returns:
    df (DataFrame): A DataFrame containing ARO accessions, original protein accessions, and mapped UniProtKB accessions.
    """

    tsv = card_data_dir + '/aro_index.tsv'
    prot_ids_col = 'Protein Accession'

    df = pd.read_csv(tsv, sep='\t')
    df['UniProtKB_acc'] = np.nan
    df['UniProtKB_acc'] = df['UniProtKB_acc'].astype('object')

    from_to = [
    ('EMBL-GenBank-DDBJ_CDS', 'UniProtKB'),
    ('RefSeq_Protein', 'UniProtKB')
    ]
    
    sg_array = []

    for from_db, to_db in from_to:
        pids = df.loc[df['UniProtKB_acc'].isnull()][prot_ids_col].to_list()

        jobid = run_idmapping(
            accesion=pids, 
            from_db=from_db, 
            to_db=to_db
            )
        map_df = get_idmapping_results(jobid=jobid)
    
        df = df.merge(
            map_df, 
            how='left', 
            left_on=prot_ids_col, 
            right_on='from')
        
        df.loc[df['UniProtKB_acc'].isnull(), 'UniProtKB_acc'] = df['to']
        df.drop(columns=['from', 'to'], inplace=True)

        ### get suggested ids
        sg_array.append(get_idmapping_suggestedIds(jobid=jobid))

    ### get UniProtKB from UPIs
    sgmap_df = pd.concat(sg_array)
    sgmap_df.reset_index(inplace=True, drop=True)
    
    upis = sgmap_df['to'].to_list()
    jobid = run_idmapping(
            accesion=upis, 
            from_db='UniParc', 
            to_db='UniProtKB'
            )
    map_upis = get_idmapping_results(jobid=jobid)
    sgmap_df = sgmap_df.merge(
                            map_upis, 
                            how='left', 
                            left_on='to', 
                            right_on='from', 
                            suffixes=('', '_upi'))
    sgmap_df.loc[sgmap_df['to_upi'].notnull(), 'to'] = sgmap_df['to_upi']
    
    df = df.merge(
                sgmap_df, 
                how='left', 
                left_on=prot_ids_col,
                right_on='from')
    df.loc[df['UniProtKB_acc'].isnull(), 'UniProtKB_acc'] = df['to']
    df.drop(columns=['from', 'to'], inplace=True)

    return df[['ARO Accession', prot_ids_col, 'UniProtKB_acc']]

def mapping_armfinderplus_protein_accession(amrfinderplus_dir):
    """
    This function maps AMRFinderPlus protein accessions to UniProtKB accessions using the UniProt ID mapping API. 
    It reads the AMRFinderPlus data from a TSV file, performs the mapping, and returns a DataFrame with the results.

    Parameters:
    amrfinderplus_dir (str): Path to the directory containing AMRFinderPlus data.

    Returns:
    df (DataFrame): A DataFrame containing original protein accessions, mapped UniProtKB accessions.
    """
    
    tsv = amrfinderplus_dir + '/ReferenceGeneCatalog.txt'
    
    df = pd.read_csv(tsv, sep='\t')
    df['UniProtKB_acc'] = np.nan
    df['UniProtKB_acc'] = df['UniProtKB_acc'].astype('object')

    from_to = [
    ('RefSeq_Protein', 'UniProtKB')
    ]

    sg_array = []

    for from_db, to_db in from_to:
        pids = df.loc[df['UniProtKB_acc'].isnull()][prot_ids_col].to_list()

        jobid = run_idmapping(
            accesion=pids, 
            from_db=from_db, 
            to_db=to_db
            )
        map_df = get_idmapping_results(jobid=jobid)
    
        df = df.merge(
            map_df, 
            how='left', 
            left_on=prot_ids_col, 
            right_on='from')
        
        df.loc[df['UniProtKB_acc'].isnull(), 'UniProtKB_acc'] = df['to']
        df.drop(columns=['from', 'to'], inplace=True)

        ### get suggested ids
        sg_array.append(get_idmapping_suggestedIds(jobid=jobid))

    ### get UniProtKB from UPIs
    sgmap_df = pd.concat(sg_array)
    sgmap_df.reset_index(inplace=True, drop=True)

    upis = sgmap_df['to'].to_list()
    jobid = run_idmapping(
            accesion=upis, 
            from_db='UniParc', 
            to_db='UniProtKB'
            )
    map_upis = get_idmapping_results(jobid=jobid)

    sgmap_df = sgmap_df.merge(
                            map_upis, 
                            how='left', 
                            left_on='to', 
                            right_on='from', 
                            suffixes=('', '_upi'))
    sgmap_df.loc[sgmap_df['to_upi'].notnull(), 'to'] = sgmap_df['to_upi']
    
    df = df.merge(
                sgmap_df, 
                how='left', 
                left_on=prot_ids_col,
                right_on='from')
    df.loc[df['UniProtKB_acc'].isnull(), 'UniProtKB_acc'] = df['to']
    df.drop(columns=['from', 'to'], inplace=True)

    return df[[prot_ids_col, 'UniProtKB_acc']]

def mapping_resfinder_protein_accession(resfinder_dir):
    """
    This function maps ResFinder protein accessions to UniProtKB accessions using the UniProt ID mapping API. 
    It reads the ResFinder data from a TSV file, performs the mapping, and returns a DataFrame with the results.

    Parameters:
    resfinder_dir (str): Path to the directory containing ResFinder data.

    Returns:
    df (DataFrame): A DataFrame containing original protein accessions and mapped UniProtKB accessions.
    """
    
    print('resfinder')
    tsv = resfinder_dir + '/phenotypes.txt'
    prot_ids_col = 'resfinder_accession'
    df = pd.read_csv(tsv, sep='\t')
    df[prot_ids_col] = df['Gene_accession no.'].str.split('_').apply(lambda x: x[-1])
    for i in range(0, df.shape[0]):
            prefix = df.loc[i,'Gene_accession no.'].split('_')[-2]
            if prefix == 'NG' or prefix == 'NC':
                df.loc[i,prot_ids_col] = f'{prefix}_{df.loc[i,'Gene_accession no.'].split('_')[-1]}'

    df['UniProtKB_acc'] = np.nan
    df['UniProtKB_acc'] = df['UniProtKB_acc'].astype('object')

    from_to = [
    ('EMBL-GenBank-DDBJ', 'UniProtKB'),
    ('RefSeq_Nucleotide', 'UniProtKB')
    ]

    sg_array = []

    for from_db, to_db in from_to:
        pids = df.loc[df['UniProtKB_acc'].isnull()][prot_ids_col].to_list()
        # print(pids)

        jobid = run_idmapping(
            accesion=pids, 
            from_db=from_db, 
            to_db=to_db
            )
        map_df = get_idmapping_results(jobid=jobid)
        
        df = df.merge(
            map_df, 
            how='left', 
            left_on=prot_ids_col, 
            right_on='from')
        
        df.loc[df['UniProtKB_acc'].isnull(), 'UniProtKB_acc'] = df['to']
        df.drop(columns=['from', 'to'], inplace=True)

        ### get suggested ids
        sg = get_idmapping_suggestedIds(jobid=jobid)
        if sg:
            sg_array.append()


    ### get UniProtKB from UPIs
    if len(sg_array) > 0:
        sgmap_df = pd.concat(sg_array)
        sgmap_df.reset_index(inplace=True, drop=True)

        upis = sgmap_df['to'].to_list()
        jobid = run_idmapping(
                accesion=upis, 
                from_db='UniParc', 
                to_db='UniProtKB'
                )
        map_upis = get_idmapping_results(jobid=jobid)

        sgmap_df = sgmap_df.merge(
                                map_upis, 
                                how='left', 
                                left_on='to', 
                                right_on='from', 
                                suffixes=('', '_upi'))
        sgmap_df.loc[sgmap_df['to_upi'].notnull(), 'to'] = sgmap_df['to_upi']
        
        df = df.merge(
                    sgmap_df, 
                    how='left', 
                    left_on=prot_ids_col,
                    right_on='from')
        df.loc[df['UniProtKB_acc'].isnull(), 'UniProtKB_acc'] = df['to']
        df.drop(columns=['from', 'to'], inplace=True)

    return df[[prot_ids_col, 'UniProtKB_acc']]



