import os

current_dir = os.path.dirname(os.path.abspath(__file__))
 
BASE_DIR = os.path.dirname(current_dir)
DATA_DIR = os.path.join(BASE_DIR, 'databases')
CODE_DIR = os.path.join(BASE_DIR, 'code')