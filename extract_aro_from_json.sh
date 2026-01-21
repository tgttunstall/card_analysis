#!/usr/bin/env bash

BASEDIR="${HOME}"
DATADIR="${BASEDIR}/amr/databases/card/data"
CARD_DATA="${DATADIR}/card.json"
OUTDIR="${DATADIR}/ind_json_output"

echo 
echo "Reading card data file: ${CARD_DATA}"
echo 
echo "Output dir: ${OUTDIR}"

myaro="3007637"
OUTFILE="${OUTDIR}/${myaro}.json"
mkdir -pv ${OUTDIR}

echo "Output file written to: ${OUTFILE}"

jq --arg aro "${myaro}" '
  .[]
  | objects
  | select(.ARO_accession? == $aro)
' ${CARD_DATA} > ${OUTFILE}

# direct cmd test
#jq '                
#  .[]
#  | objects
#  | select(.ARO_accession? == "3007637")
#' card.json 
