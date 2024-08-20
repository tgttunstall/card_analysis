#!/usr/bin/env python
# coding: utf-8 

import pandas as pd
import networkx as nx
import argparse
import webbrowser
from pyvis.network import Network
from common_functions import card_graph
from common_functions import amrfinderplus_graph
from common_functions import resfinder_graph
from common_functions import simplify_graph
from common_functions import insert_newline_every_n_chars
from common_functions import create_pyvis_html
from common_functions import graph_to_cytoscape_json
from common_functions import graph_to_cytoscape_desktop_json
from common_functions import save_graph_as_png
from config import BASE_DIR, DATA_DIR, CARD_DIR, MAP_DIR
from pathlib import Path
import json
import py4cytoscape as p4c

parser = argparse.ArgumentParser(description='Create knowledge graph from UniProt KB accession with CARD ontology')
parser.add_argument('--accession', type=str, help='UniProt KB accession', required=True)
parser.add_argument('--outdir', type=str, help='Path to html file with the knowledge graph', required=True)
parser.add_argument('--visualization', type=str, choices=['pyvis', 'cytoscape', 'cytoscape_desktop', 'png'], default='pyvis', help='Visualization tool to use (pyvis or cytoscape)', required=True)
args = parser.parse_args()

acc = args.accession
outdir = Path(args.outdir)

### Define group color mapping
group_colors = {
    'card': '#0000FF', #blue
    'uniprot': 'red',
    'amrfinderplus': 'green',
    'resfinder': 'orange'
    # Add more groups and their colors as needed
}
card_colors = {
    'card': '#0000FF', # blue
    'AMR Gene Family': '#8A2BE2', # violet blue
    'Drug Class': '#00BFFF', # Depskyblue
    'Resistance Mechanism': '#B0E0E6' # pale blue
}

### create graph
print(f'>>> {acc} <<<')
G = nx.MultiDiGraph()

### Create CARD KG
card_G, aro = card_graph(
    obo_file=CARD_DIR / 'ontology/aro.obo',
    json_file=CARD_DIR / 'data/card.json',
    categories_file=CARD_DIR / 'data/aro_categories.tsv',
    map_file=MAP_DIR / 'card_map.tsv',
    acc=acc, 
    colors=card_colors
    )
### Create AMRFinderPlus graph
amrfinderplus_G, gene_fam = amrfinderplus_graph(
    database_dir=DATA_DIR / 'AMRFinderPlus', 
    map_file=MAP_DIR / 'amrfinderplus_map.tsv',
    acc=acc,
    color=group_colors.get('amrfinderplus', 'black')  
    )

### Create AMRFinderPlus graph
resfinder_G, gene_accession = resfinder_graph(
    phenotype=DATA_DIR / 'resfinder_db', 
    map_file=MAP_DIR / 'resfinder_map.tsv',
    acc=acc,
    color=group_colors.get('resfinder', 'black')  
    )

### Add UniProt accession node and edge manually
G.add_node(acc, **{
    'name': acc,
    'def': '',
    'label': acc,
    'title': acc,
    'group': 'uniprot',
    'color': group_colors.get('uniprot', 'black')
})

### graphs union
if card_G is not None:
    G = nx.union(G, card_G)
    G.add_edge(acc, aro, label='is')
if amrfinderplus_G is not None:
    G = nx.union(G, amrfinderplus_G)
    G.add_edge(acc, gene_fam, label='is')
if resfinder_G is not None:
    G = nx.union(G, resfinder_G)
    G.add_edge(acc, gene_accession, label='is')

# Customize the nodes and edges
for node_id in G.nodes:
    in_degree = G.in_degree(node_id)  # Get the number of incoming edges

    # Get existing attributes and modify or add new ones
    G.nodes[node_id]['label'] = insert_newline_every_n_chars(G.nodes[node_id].get('label', ''), 35)
    G.nodes[node_id]['title'] = insert_newline_every_n_chars(G.nodes[node_id].get('title', ''), 80)
    G.nodes[node_id]['color'] = G.nodes[node_id].get('color', 'black')
    G.nodes[node_id]['font_size'] = 45 + in_degree * 2  # Increase font size based on in-degree
    G.nodes[node_id]['size'] = 25 + in_degree * 6  # Adjust size based on in-degree
    G.nodes[node_id]['font_color'] = 'white'  # Change font color to white

for u, v, key, data in G.edges(keys=True, data=True):
    data['width'] = 2
    data['font_size'] = 18  # Decrease font size
    data['font_face'] = 'arial'
    data['font_color'] = 'gray'  # Change font color to gray

    if data.get('label') == 'is_a':
        data['color'] = 'olive'
        data['font_color'] = 'olivedrab'
    elif 'confers_resistance_to' in data.get('label', ''):
        data['color'] = 'firebrick'
        data['font_color'] = 'indianred'

### merge nodes with same predecesor and sucesor nodes !!! on develop
# G = simplify_graph(graph=G)


### Generate the pyvis network graph
net = create_pyvis_html(graph=G)
# net.show_buttons()
# Redefine the open method of webbrowser to do nothing
webbrowser.open = lambda url, new=0, autoraise=True: None
net.save_graph(str(outdir / f"{acc}.html"))

if args.visualization == 'cytoscape':
    # Generate the Cytoscape.js compatible JSON
    cytoscape_json = graph_to_cytoscape_json(G)
    
    # Save the JSON to a file
    json_file = outdir / f"{acc}_cytoscape_js.json"
    with open(json_file, 'w') as f:
        json.dump(cytoscape_json, f)
    
    print(f'Graph saved in Cytoscape.js JSON format at {json_file}')

elif args.visualization == 'cytoscape_desktop':
    # Generate the Cytoscape desktop compatible JSON
    cytoscape_desktop_json = graph_to_cytoscape_desktop_json(G)
    
    # Save the JSON to a file
    json_file = outdir / f"{acc}_cytoscape_desktop.json"
    with open(json_file, 'w') as f:
        json.dump(cytoscape_desktop_json, f)
    
    print(f'Graph saved in Cytoscape desktop JSON format at {json_file}')


png_file = outdir / f"{acc}.png"
save_graph_as_png(graph=G, file=png_file, layout='spring')
print(f'Graph saved in PNG format at {png_file}')