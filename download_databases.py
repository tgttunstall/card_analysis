#!/usr/bin/env python
# coding: utf-8

from config import CARD_DIR, AMRFINDERPLUS_DIR, RESFINDER_DIR, POINTFINDER_DIR
from common_functions import download_and_extract_file
from common_functions import download_ftp_directory
from common_functions import clone_repo

def main():

  print('downloading CARD data')
  download_and_extract_file(
      url='https://card.mcmaster.ca/latest/data', 
      download_dir=CARD_DIR + '/data/')

  print('downloading CARD ontology')
  download_and_extract_file(
      url='https://card.mcmaster.ca/latest/ontology', 
      download_dir=CARD_DIR + '/ontology/')

  print('downloading ARMFinderPlus')
  download_ftp_directory(
      ftp_url='ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/', 
      download_dir=AMRFINDERPLUS_DIR,
      exclude_files=['AMR.LIB']
      )

  print('cloning ResFinder')  
  clone_repo(
      bitbucket_url='https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/',
      clone_dir=RESFINDER_DIR
      )

  print('downloading PointFinder')
  clone_repo(
      bitbucket_url='https://bitbucket.org/genomicepidemiology/pointfinder_db/src/master/',
      clone_dir=POINTFINDER_DIR
      )

if __name__ == "__main__":
  main()

