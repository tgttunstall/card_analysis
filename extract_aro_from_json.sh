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
OUTFILE="${OUTDIR}/aro_${myaro}.json"
mkdir -pv ${OUTDIR}


#jq --arg aro "${myaro}" '
#  .[]
#  | objects
#  | select(.ARO_accession? == $aro)
#' ${CARD_DATA} > ${OUTFILE}

echo 
echo "Sanity check to see if ARO is unique"
echo 
jq --arg aro "${myaro}" '
  [.[] | objects | select(.ARO_accession? == $aro)] | length
' "${CARD_DATA}"
############################
# fail if no match
jq -e --arg aro "${myaro}" '
  .[] 
  | objects 
  | select(.ARO_accession? == $aro) # not a bash variable
' "${CARD_DATA}" > "${OUTFILE}" || {
  echo "ERROR: ARO ${myaro} not found" >&2
  exit 1
}
#FIXME: output file still gets written even if jq fails
echo "Output file written to: ${OUTFILE}"
#############################
#Extract specific fields: to ensure it valid json, wrap in array

jq -e --arg aro "${myaro}" '
 [
  .[]
  | objects
  | select(.ARO_accession? == $aro)
  | .ARO_category? // empty
  | .[]
  | select (
  	.category_aro_class_name == "AMR Gene Family"
	or .category_aro_class_name == "Drug Class"
	or .category_aro_class_name == "Resistance Mechanism"
    )
 ]
' "${CARD_DATA}" >  "${OUTDIR}/aro_${myaro}_sfields.json" || {
  echo "ERROR: ARO ${myaro} not found" >&2
  exit 1
}

echo "Output file for specific fields written to: ${OUTDIR}/aro_${myaro}_sfields.json"

#############################
# FILE 2: aro_categories.tsv
cut -f1 ${DATADIR}/aro_categories.tsv | sort | uniq -c > ${OUTDIR}/aro_categories_type.txt


############################
# direct cmd test
#jq '                
#  .[]
#  | objects
#  | select(.ARO_accession? == "3007637")
#' card.json 
