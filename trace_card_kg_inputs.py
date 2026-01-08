#!/usr/bin/env python3

import argparse
import pandas as pd
import json
# from common_functions import parse_obo_file, extract_subgraph, card_colors
from common_functions import parse_obo_file, extract_subgraph

card_colors = {
    'card': 'blue',
    'Antibiotic': 'rebeccapurple',
    'Drug Class': 'mediumorchid',
    'AMR Gene Family': 'steelblue',
    'Resistance Mechanism': 'deepskyblue'
}


def trace_inputs(accession, map_file, obo_file, json_file, categories_file, outdir):
    print(f"\n🔍 Tracing node-centric CARD KG data for UniProtKB accession: {accession}")

    # Step 1: Get root ARO
    map_df = pd.read_csv(map_file, sep="\t")
    acc_row = map_df.loc[map_df.iloc[:, 0] == accession]
    if acc_row.empty:
        print(f"❌ Accession {accession} not found in mapping file.")
        return
    aro_root = acc_row["ARO"].iloc[0]
    print(f"✅ Root ARO: {aro_root}")

    # Step 2: Load data
    obo_graph = parse_obo_file(obo_file)
    graph = extract_subgraph(obo_graph, aro_root)
    cat_df = pd.read_csv(categories_file, sep="\t")
    cat_map = dict(zip(cat_df["ARO Accession"], cat_df["ARO Category"]))
    with open(json_file, "r") as f:
        card_data = json.load(f)

    # Step 3: Identify antibiotic nodes
    antibiotic_nodes = set()
    for src, tgt, label in graph.edges:
        if label == "confers_resistance_to_antibiotic":
            antibiotic_nodes.add(tgt)

    # Step 4: Collect node and edge summary per node
    rows = []
    for node, data in graph.nodes(data=True):
        name = data.get("name", node)
        #category = ""
        #source = "default"
        #color = card_colors["card"]

        #if node in antibiotic_nodes:
        #    category = "Antibiotic"
        #    source = "obo edge: confers_resistance_to_antibiotic"
        #    color = card_colors["Antibiotic"]
        #elif node in cat_map:
        #    category = cat_map[node]
        #    source = "aro_categories.tsv"
        #    color = card_colors.get(category, color)
        #elif any(isinstance(entry, dict) and str(entry.get("ARO_accession")) == node for entry in card_data.values()):
         #   category = "From card.json"
         #   source = "card.json"

        # Track all known sources
        sources = []
        category = ""
        color = card_colors["card"]

        # Check for ontology presence
        if node in obo_graph.nodes:
            sources.append("aro.obo")

        # Check for antibiotic-related nodes (by edge)
        if node in antibiotic_nodes:
            category = "Antibiotic"
            color = card_colors["Antibiotic"]
            if "aro.obo" not in sources:
                sources.append("aro.obo")  # came from edge in ontology

        # Check for category file match
        if node in cat_map:
            category = cat_map[node]
            color = card_colors.get(category, color)
            sources.append("aro_categories.tsv")

        # Check for card.json match
        in_card_json = any(
            isinstance(entry, dict) and str(entry.get("ARO_accession")) == node
            for entry in card_data.values()
        )
        if in_card_json:
            sources.append("card.json")

        # Deduplicate and join sources
        source = ";".join(sorted(set(sources)))         

        # Edge summaries
        edges = list(graph.edges(node, keys=True))
        targets = [t for (_, t, _) in edges]
        labels = [k for (_, _, k) in edges]
        target_names = [graph.nodes[t].get("name", t) for t in targets]

        rows.append({
            "UniProtKB": accession,
            "ARO": node,
            "Name": name,
            "Category": category,
            "Source": source,
            "Color": color,
            "Edge Targets (ARO)": ";".join(targets) if targets else "",
            "Edge Labels": ";".join(labels) if labels else "",
            "Target Names": ";".join(target_names) if target_names else ""
        })

    # Step 5: Output
    df = pd.DataFrame(rows)
    #out_file = f"traced_output_{accession}.csv"
    out_file = outdir + f"/traced_output_{accession}.csv"

    df.to_csv(out_file, index=False)

    print(f"\n📄 Node summary (with edge targets) written to: {out_file}")
    print(f"🔢 Total nodes: {len(df)}")
    print(f"📊 Categories: {df['Category'].value_counts().to_dict()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trace node-level data for a CARD KG.")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--map_file", required=True)
    parser.add_argument("--obo_file", required=True)
    parser.add_argument("--json_file", required=True)
    parser.add_argument("--categories_file", required=True)
    parser.add_argument("--outdir", default=".", help="Directory to write output CSV file (default: current directory)")

    args = parser.parse_args()

    trace_inputs(
        accession=args.accession,
        map_file=args.map_file,
        obo_file=args.obo_file,
        json_file=args.json_file,
        categories_file=args.categories_file,
        outdir=args.outdir
    )

