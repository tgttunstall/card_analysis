#!/usr/bin/env bash
#!/usr/bin/env bash

BASEDIR="$HOME"
DATADIR="$BASEDIR/amr/databases/card/data"
CARD_DATA="$DATADIR/card.json"
OUTDIR="$DATADIR/ind_json_output"

mkdir -pv "$OUTDIR"

myaro="3007637"
OUTFILE="$OUTDIR/aro_${myaro}_filtered.json"

echo "Reading card data file: $CARD_DATA"
echo "Output file: $OUTFILE"

tmpfile="$(mktemp)"

# Filter the JSON
jq --arg aro "$myaro" '
  [
    .[]
    | objects
    | select(.ARO_accession? == $aro)
    | {
        model_id,
        model_name,
        model_type,
        ARO_category: (
          .ARO_category? // {}
          | with_entries(
              select(
                .value.category_aro_class_name == "AMR Gene Family"
                or .value.category_aro_class_name == "Drug Class"
                or .value.category_aro_class_name == "Resistance Mechanism"
              )
            )
        )
      }
  ]
' "$CARD_DATA" > "$tmpfile"

# Check if any matches were found
if [[ $(jq 'length' "$tmpfile") -eq 0 ]]; then
    rm -f "$tmpfile"
    echo "ERROR: ARO $myaro not found or no matching categories" >&2
    exit 1
fi

# Optional: info if multiple models exist for the ARO
num_models=$(jq --arg aro "$myaro" '[.[] | objects | select(.ARO_accession? == $aro)] | length' "$CARD_DATA")
if [[ "$num_models" -gt 1 ]]; then
    echo "INFO: Multiple models ($num_models) found for ARO $myaro"
fi

mv "$tmpfile" "$OUTFILE"
echo "Filtered JSON written to $OUTFILE"



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
