import sys
from os import getcwd
from pathlib import Path
sys.path.append(str(Path(getcwd()).parent.joinpath('src')))
from src.helper_functions import add_metadata_to_barcodes, write_features_and_mtx, create_h5ad

#Existing SJs to mtx
output_dir = "../results/SJ_from_existing_BAMs"
with open("../results/directory_list.txt", 'r') as f:
    lines = f.readlines()
directory_list = [line.strip('\n') for line in lines]
add_metadata_to_barcodes(
    [Path(directory).name for directory in directory_list],
    output_dir)
write_features_and_mtx(directory_list, output_dir, "_SJ.out.tab")
create_h5ad(output_dir, "features_ann.tsv")