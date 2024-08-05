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
from config import BASE_DIR, DATA_DIR, CARD_DIR, MAP_DIR

parser = argparse.ArgumentParser(description='Create knowledge graph from UniProt KB accession with CARD ontology')
parser.add_argument('--accession', type=str, help='UniProt KB accession', required=True)
parser.add_argument('--outdir', type=str, help='Path to html file with the knowledge graph', required=True)
args = parser.parse_args()

acc = args.accession
outdir = args.outdir

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
	obo_file=CARD_DIR + '/ontology/aro.obo',
  json_file=CARD_DIR + '/data/card.json',
  categories_file=CARD_DIR + '/data/aro_categories.tsv',
	map_file=MAP_DIR + '/card_map.tsv',
	acc=acc, 
	colors=card_colors
	)
### Create AMRFinderPlus graph
amrfinderplus_G, gene_fam = amrfinderplus_graph(
  database_dir=DATA_DIR + '/AMRFinderPlus', 
  map_file=BASE_DIR + '/map_tsv/amrfinderplus_map.tsv',
  acc=acc,
  color=group_colors.get('amrfinderplus', 'black')  
  )

### Create AMRFinderPlus graph
resfinder_G, gene_accession = resfinder_graph(
  phenotype=DATA_DIR + '/resfinder_db', 
  map_file=BASE_DIR + '/map_tsv/resfinder_map.tsv',
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

### merge nodes with same predecesor and sucesor nodes !!! on develop
# G = simplify_graph(graph=G)

# Create a PyVis network
net = Network(notebook=False, height='800px', width='100%', bgcolor='#222222', font_color='white')

# Load the NetworkX graph into the PyVis network
net.from_nx(G)

# Customize the nodes and edges
for node in net.nodes:
  # print(node)
  # node['title'] = insert_newline_every_n_chars(node['title'], 80)  # Display label on hover
  node['label'] = insert_newline_every_n_chars(node.get('label', ''), 80)  # Show wrap label on the node
  node['title'] = insert_newline_every_n_chars(node.get('title', ''), 80)  # Show wrap title on the node
  node['color'] = G.nodes[node['id']].get('color', 'black') # get color from networkx graph node
  node['font'] = {'size': 45, 'color': 'white', 'vadjust': 0, 'multi': 'html'}  # Increase font size and change color to white
  node['size'] = 25  # Increase node size

for edge in net.edges:
  edge['title'] = insert_newline_every_n_chars(edge.get('title', ''), 80)  # Display edge label on hover
  edge['color'] = 'gray'
  edge['width'] = 2
  edge['font'] = {'size': 18, 'face': 'arial', 'color': 'gray'}  # Decrease font size and change color to white


# Apply physics layout for better node separation
net.set_options("""
var options = {
  "nodes": {
    "color": {
      "highlight": {
        "border": "white",
        "background": "black"
      },
      "hover": {
        "border": "white",
        "background": "black"
      }
    }
  },
  "edges": {
    "color": {
      "color": "gray"
    },
    "font": {
      "color": "gray",
      "size": 16,
      "face": "arial",
      "background": "none",
      "strokeWidth": 0,
      "strokeColor": "none",
      "multi": true
    },
    "smooth": true,
    "arrows": {
      "to": {
        "enabled": true,
        "scaleFactor": 1
      }
    }
  },
  "physics": {
    "barnesHut": {
      "gravitationalConstant": -20000,
      "centralGravity": 0.3,
      "springLength": 200,
      "springConstant": 0.01,
      "damping": 0.09
    },
    "minVelocity": 0.75
  }
}
""")

# Generate the network graph
# net.show_buttons()
# Redefine the open method of webbrowser to do nothing
webbrowser.open = lambda url, new=0, autoraise=True: None
net.save_graph(f"{outdir}/{acc}.html")


