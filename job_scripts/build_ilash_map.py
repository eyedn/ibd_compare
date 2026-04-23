###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_ilash_map.py
###############################################################################


from pathlib import Path
from sys import argv

from python_utils.prepare_ilash_maps_utils import build_ilash_map

if len(argv) != 4:
    raise SystemExit(
        "Usage: build_ilash_map.py <reference_map> <ped_map> <output_map>"
    )

reference_map = Path(argv[1])
ped_map = Path(argv[2])
output_map = Path(argv[3])

build_ilash_map(reference_map, ped_map, output_map)
