# Data source: https://www.10xgenomics.com/datasets/preview-data-xenium-prime-gene-expression

import spatialdata as sd
from spatialdata import read_zarr
from spatialdata_io import xenium
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq

xenium_path = "./Raw/Xenium_Prime_Human_Lymph_Node_Reactive_FFPE_outs"
sdata = xenium(xenium_path)
zarr_path = "./Derived/Xenium_5K_Human_Lymph_Node"
sdata.write(zarr_path+'/20250609_data.zarr')
sdata = read_zarr(zarr_path+'/20250609_data.zarr')

sdata

SpatialData object, with associated Zarr store: /Users/xxx/Downloads/Xenium_5K/Data/Derived/Xenium_5K_Human_Lymph_Node/20250609_data.zarr
├── Images
│     ├── 'he_image': DataTree[cyx] (3, 30987, 28915), (3, 15493, 14457), (3, 7746, 7228), (3, 3873, 3614), (3, 1936, 1807)
│     └── 'morphology_focus': DataTree[cyx] (5, 34119, 39776), (5, 17059, 19888), (5, 8529, 9944), (5, 4264, 4972), (5, 2132, 2486)
├── Labels
│     ├── 'cell_labels': DataTree[yx] (34119, 39776), (17059, 19888), (8529, 9944), (4264, 4972), (2132, 2486)
│     └── 'nucleus_labels': DataTree[yx] (34119, 39776), (17059, 19888), (8529, 9944), (4264, 4972), (2132, 2486)
├── Points
│     └── 'transcripts': DataFrame with shape: (<Delayed>, 13) (3D points)
├── Shapes
│     ├── 'cell_boundaries': GeoDataFrame shape: (708983, 1) (2D shapes)
│     ├── 'cell_circles': GeoDataFrame shape: (708983, 2) (2D shapes)
│     └── 'nucleus_boundaries': GeoDataFrame shape: (702300, 1) (2D shapes)
└── Tables
      └── 'table': AnnData (708983, 4624)
with coordinate systems:
    ▸ 'global', with elements:
        he_image (Images), morphology_focus (Images), cell_labels (Labels), nucleus_labels (Labels), transcripts (Points), cell_boundaries (Shapes), cell_circles (Shapes), nucleus_boundaries (Shapes)
