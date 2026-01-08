#!/usr/bin/env bash
CARD_DIR="/home/tanu/amr/databases/card"
OUT_DIR="/home/tanu/amr"
#MAP_FILE="/home/tanu/amr/card_map.tsv"
#MAP_FILE="/home/tanu/amr/map_tsv/card_map.tsv"
MAP_FILE="/home/tanu/amr/map_tsv/CARD-UniProt-Mapping.tsv"
#############################################
#ACC="P28585"
#ACC="Q182T3"
ACC="A6T5M6"
##############################################
python generate_card_kg.py \
  --accession "${ACC}" \
  --outdir "${OUT_DIR}" \
  --map_file "${MAP_FILE}" \
  --obo_file "${CARD_DIR}/ontology/aro.obo" \
  --json_file "${CARD_DIR}/data/card.json" \
  --categories_file "${CARD_DIR}/data/aro_categories.tsv" \
  --visualization pyvis png

#  --aro_index "${CARD_DIR}/data/aro_index.tsv" \


python trace_card_kg_inputs.py \
  --accession  "${ACC}" \
  --map_file "${MAP_FILE}" \
  --obo_file "${CARD_DIR}/ontology/aro.obo" \
  --json_file "${CARD_DIR}/data/card.json" \
  --categories_file "${CARD_DIR}/data/aro_categories.tsv" \
  --outdir "${OUT_DIR}"

