from pathlib import Path
from src.helper_functions import add_metadata_to_barcodes, create_h5ad

output_dir = Path("/external/rprshnas01/netdata_kcni/stlab/Nuo/output/SJ/raw")
with open("../results/corrupted_gz_removed.txt", "r") as f:
    lines = f.readlines()
directory_list = [line.strip('\n') for line in lines]

# transform the directory list to a filemanifest.tsv
file_manifest = "\n".join(
    [(lambda directory: "\t".join([str(file) for file in Path(directory).iterdir()][::-1] +
    [Path(directory).name]))(directory_path) for directory_path in directory_list]
    )
with open("../results/file_manifest_STARsolo.tsv", 'w') as f:
    f.write(file_manifest)

#SJ
add_metadata_to_barcodes(directory_list, output_dir)

create_h5ad(output_dir, "barcodes_ann.tsv")

#Gene
add_metadata_to_barcodes(directory_list, "/external/rprshnas01/netdata_kcni/stlab/Nuo/output/Gene/raw")