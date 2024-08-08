from pathlib import Path

# Current directory
current_dir = Path(__file__).resolve().parent

# Base directory
BASE_DIR = current_dir.parent

# Specific directories
MAP_DIR = BASE_DIR / 'map_tsv'
DATA_DIR = BASE_DIR / 'databases'
CARD_DIR = DATA_DIR / 'card'
AMRFINDERPLUS_DIR = DATA_DIR / 'amrfinderplus'
RESFINDER_DIR = DATA_DIR / 'resfinder'
POINTFINDER_DIR = DATA_DIR / 'pointfinder'
CODE_DIR = BASE_DIR / 'code'
