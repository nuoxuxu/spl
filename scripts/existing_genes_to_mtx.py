from pathlib import Path
from src.helper_functions import add_metadata_to_barcodes, write_features_and_mtx
#Existing Genes to mtx
output_dir = "../results/Gene_from_existing_BAMs"
with open("../results/directory_list.txt", 'r') as f:
    lines = f.readlines()
directory_list = [line.strip('\n') for line in lines]
add_metadata_to_barcodes(
    [Path(directory).name for directory in directory_list],
    output_dir)
write_features_and_mtx(directory_list, output_dir, "_ReadsPerGene.out.tab")