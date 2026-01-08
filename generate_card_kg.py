#!/usr/bin/env python

import argparse
import networkx as nx
import json
from pathlib import Path
#os.chdir('/home/tanu/git/card_analysis')
from common_functions import card_graph, insert_newline_every_n_chars, create_pyvis_html, graph_to_cytoscape_json, graph_to_cytoscape_desktop_json, save_graph_as_png, capture_canvas_as_png

parser = argparse.ArgumentParser(description='Create CARD-only knowledge graph from UniProtKB accession.')
parser.add_argument('--accession', type=str, required=True, help='UniProtKB accession')
parser.add_argument('--outdir', type=str, required=True, help='Output directory')
parser.add_argument('--map_file', type=str, required=True, help='TSV mapping UniProtKB to ARO')
parser.add_argument('--obo_file', type=str, required=True, help='ARO ontology file (obo)')
parser.add_argument('--json_file', type=str, required=True, help='CARD JSON file')
#parser.add_argument('--aro_index', type=str, required=True, help='ARO index file')
parser.add_argument('--categories_file', type=str, required=True, help='ARO categories TSV')
parser.add_argument('--visualization', nargs='*', default=['pyvis'], choices=['pyvis', 'cytoscape', 'cytoscape_desktop', 'png', 'all'], help='Visualization formats')
args = parser.parse_args()

# Default to all visualizations if "all" is specified
if 'all' in args.visualization:
    visualizations = ['pyvis', 'cytoscape', 'cytoscape_desktop', 'png']
else:
    visualizations = args.visualization

acc = args.accession
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

# Colors
card_colors = {
    'card': 'blue',
    'Antibiotic': 'rebeccapurple',
    'Drug Class': 'mediumorchid',
    'AMR Gene Family': 'steelblue',
    'Resistance Mechanism': 'deepskyblue'
}

# Build graph
G = nx.MultiDiGraph()
print(f">>> Building graph for {acc}")
card_G, aro = card_graph(
    obo_file=args.obo_file,
    json_file=args.json_file,
    categories_file=args.categories_file,
    map_file=args.map_file,
    #aro_index=args.aro_index,
    acc=acc,
    colors=card_colors
)

if card_G is None:
    print(f"[WARN] No CARD graph found for {acc}")
    exit()

# Add UniProtKB node and connect to ARO root
G.add_node(acc, **{
    'name': acc,
    'label': acc,
    'title': acc,
    'group': 'uniprot',
    'color': 'red'
})
G = nx.union(G, card_G)
G.add_edge(acc, aro, label='is')

# Style nodes
for node in G.nodes:
    in_deg = G.in_degree(node)
    G.nodes[node]['label'] = insert_newline_every_n_chars(G.nodes[node].get('label', ''), 35)
    G.nodes[node]['title'] = insert_newline_every_n_chars(G.nodes[node].get('title', ''), 80)
    G.nodes[node]['font_size'] = 45 + in_deg * 2
    G.nodes[node]['size'] = 25 + in_deg * 6
    G.nodes[node]['font_color'] = 'white'

# Style edges
for u, v, key, data in G.edges(keys=True, data=True):
    data['width'] = 2
    data['font_size'] = 18
    data['font_face'] = 'arial'
    data['font_color'] = 'gray'
    if data.get('label') == 'is_a':
        data['color'] = 'olive'
        data['font_color'] = 'olivedrab'
    elif 'confers_resistance_to' in data.get('label', ''):
        data['color'] = 'firebrick'
        data['font_color'] = 'indianred'

# Output files using plain string paths
if 'pyvis' in visualizations:
    net = create_pyvis_html(graph=G)
    #html_file = str(outdir) + "/" + acc + ".html"
    #png_path = str(outdir) + "/" + acc + ".png"
    #net.save_graph(html_file)
    #net.write_html(html_file)
    
    #capture_canvas_as_png(html_file, png_path, browser="chrome", wait_time=5)
    
    # instead of net.write_html(html_file)
    html_path = str(outdir) + "/" + acc + ".html"
    png_path = str(outdir) + "/" + acc + ".png"

    net_html = net.generate_html(name=acc + ".html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(net_html)

    capture_canvas_as_png(html_path, png_path, browser="chrome", wait_time=5)
    print(f"[✓] Pyvis HTML and PNG saved for {acc}")

if 'png' in visualizations:
    png_file = str(outdir) + "/" + acc + "_static.png"
    save_graph_as_png(G, file=png_file)
    print(f"[✓] PNG (static) saved for {acc}")

if 'cytoscape' in visualizations:
    json_data = graph_to_cytoscape_json(G)
    json_file = str(outdir) + "/" + acc + "_cytoscape_js.json"
    with open(json_file, 'w') as f:
        json.dump(json_data, f)
    print(f"[✓] Cytoscape.js JSON saved for {acc}")

if 'cytoscape_desktop' in visualizations:
    json_data = graph_to_cytoscape_desktop_json(G)
    json_file = str(outdir) + "/" + acc + "_cytoscape_desktop.json"
    with open(json_file, 'w') as f:
        json.dump(json_data, f)
    print(f"[✓] Cytoscape desktop JSON saved for {acc}")

