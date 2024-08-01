import os

current_dir = os.path.dirname(os.path.abspath(__file__))
 
BASE_DIR = os.path.dirname(current_dir)
DATA_DIR = os.path.join(BASE_DIR, 'databases')
CARD_DIR = os.path.join(DATA_DIR, 'card')
AMRFINDERPLUS_DIR = os.path.join(DATA_DIR, 'amrfinderplus')
RESFINDER_DIR = os.path.join(DATA_DIR, 'resfinder')
POINTFINDER_DIR = os.path.join(DATA_DIR, 'pointfinder')

CODE_DIR = os.path.join(BASE_DIR, 'code')