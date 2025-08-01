#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from config import MAP_DIR, CARD_DIR, AMRFINDERPLUS_DIR, RESFINDER_DIR, POINTFINDER_DIR
from common_functions import id_mapping_protein_accessions
from common_functions import complete_id_mapping_with_api_results
 
def main():
    ########
    ### CARD
    ### create dataframe from file
    card_data_file = CARD_DIR + '/data/aro_index.tsv'
    df = pd.read_csv(card_data_file, sep='\t')
    df['UniProtKB_acc'] = np.nan
    df['UniProtKB_acc'] = df['UniProtKB_acc'].astype('object')
    ### from to databases
    card_from_to = [
    ('EMBL-GenBank-DDBJ_CDS', 'UniProtKB'),
    ('RefSeq_Protein', 'UniProtKB')
    ]
    ### run id mapping
    df = id_mapping_protein_accessions(
        df=df, 
        prot_ids_col='Protein Accession',
        from_to=card_from_to
        )
    ### check accession not found with id mapping
    df = complete_id_mapping_with_api_results(df)
    df.to_csv(MAP_DIR + '/card_map.tsv', sep='\t', index=False)

    #################
    ### AMRFinderPlus
    ### create dataframe from file
    amrfinderplis_data_file = AMRFINDERPLUS_DIR + '/ReferenceGeneCatalog.txt'
    df = pd.read_csv(amrfinderplis_data_file, sep='\t')
    df['UniProtKB_acc'] = np.nan
    df['UniProtKB_acc'] = df['UniProtKB_acc'].astype('object')
    ### from to databases
    amrfinderplus_from_to = [
    ('RefSeq_Protein', 'UniProtKB')
    ]
    ### run id mapping
    df = id_mapping_protein_accessions(
        df=df, 
        prot_ids_col='refseq_protein_accession',
        from_to=amrfinderplus_from_to
        )
    ### check accession not found with id mapping
    df = complete_id_mapping_with_api_results(df)
    df.to_csv(MAP_DIR + '/amrfinderplus_map.tsv', sep='\t', index=False)

    #############
    ### ResFinder
    ### create dataframe from file
    resfinder_data_file = RESFINDER_DIR + '/phenotypes.txt'
    df = pd.read_csv(resfinder_data_file, sep='\t')
    # creating database id columns
    df['resfinder_accession'] = df['Gene_accession no.'].str.split('_').apply(lambda x: x[-1])
    for i in range(0, df.shape[0]):
            prefix = df.loc[i,'Gene_accession no.'].split('_')[-2]
            if prefix == 'NG' or prefix == 'NC':
                df.loc[i,'resfinder_accession'] = f'{prefix}_{df.loc[i,'Gene_accession no.'].split('_')[-1]}'
    df['UniProtKB_acc'] = np.nan
    df['UniProtKB_acc'] = df['UniProtKB_acc'].astype('object')
    ### from to databases
    resfinder_from_to = [
    ('EMBL-GenBank-DDBJ', 'UniProtKB'),
    ('RefSeq_Nucleotide', 'UniProtKB')
    ]
    ### run id mapping
    df = id_mapping_protein_accessions(
        df=df, 
        prot_ids_col='resfinder_accession',
        from_to=resfinder_from_to
        )
    ### check accession not found with id mapping
    df = complete_id_mapping_with_api_results(df)
    df.to_csv(MAP_DIR + '/resfinder_map.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()











