
import os
import pandas as pd
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from scipy.io import mmread
from scipy.stats import ttest_ind
import json
base_path = "X:/projects/bsc83/Data/spatial/Myometrium_NatCom2023/extracted1/"
# --- 1. Load Metadata ---
metadata = pd.read_csv(os.path.join(base_path,"metadata.csv"))
# columns: folder_name, Status (e.g., "peri", "post"), [other columns]
genes = None  # Will be set after first load
all_samples = []

# --- 2. Load Data for All Samples ---
for _, row in metadata.iterrows():
    folder = row['folder_name']
    folder_path = os.path.join(base_path,folder)
    spatial_path = os.path.join(folder_path, "spatial")
    matrix_path = os.path.join(folder_path, "filtered_feature_bc_matrix")
    
    # Load matrix
    matrix = mmread(os.path.join(matrix_path, "matrix.mtx.gz")).T.tocsr()
    barcodes = pd.read_csv(os.path.join(matrix_path, "barcodes.tsv.gz"), header=None)[0].tolist()
    features = pd.read_csv(os.path.join(matrix_path, "features.tsv.gz"), sep="\t", header=None)
    if genes is None:
        genes = features[1].tolist()
    
    # Load positions
    positions = pd.read_csv(
        os.path.join(spatial_path, "tissue_positions.csv"),
        sep=",", header=0,
       # names=["barcode", "tissue", "row", "col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    )
    positions = positions[positions['in_tissue'] == 1]
    
    # Load scale factors, image
    with open(os.path.join(spatial_path, "scalefactors_json.json")) as f:
        scale_factors = json.load(f)
    lowres_scale = scale_factors['tissue_lowres_scalef']
    lowres_img = Image.open(os.path.join(spatial_path, "tissue_lowres_image.png"))
    
    # Scale to lowres image
    positions['x_lowres'] = positions['pxl_col_in_fullres'] * lowres_scale
    positions['y_lowres'] = positions['pxl_row_in_fullres'] * lowres_scale
    
    # Match barcodes
    barcode_to_index = {b: i for i, b in enumerate(barcodes)}
    positions = positions[positions['barcode'].isin(barcode_to_index)]
    indices = [barcode_to_index[bc] for bc in positions['barcode']]
    
    # Per-sample mean expression
    mean_expr = np.array(matrix[indices].mean(axis=0)).flatten()
    
    all_samples.append({
        'folder': folder,
        'Status': row['Status'],
        'positions': positions,
        'matrix': matrix,
        'indices': indices,
        'img': lowres_img,
        'mean_expr': mean_expr
    })

def compare_spatial_with_shared_scale(peri_sample, post_sample, gene, genes):
    gene_idx = genes.index(gene)
    peri_expr = peri_sample['matrix'][peri_sample['indices'], gene_idx].toarray().flatten()
    post_expr = post_sample['matrix'][post_sample['indices'], gene_idx].toarray().flatten()
    all_expr = np.concatenate([peri_expr, post_expr])
    vmin, vmax = all_expr.min(), all_expr.max()

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    axs[0].imshow(peri_sample['img'])
    im0 = axs[0].scatter(peri_sample['positions']['x_lowres'], peri_sample['positions']['y_lowres'], c=peri_expr, cmap='hot', s=25, vmin=vmin, vmax=vmax)
    axs[0].set_title(f'PERI: {gene}({peri_sample["folder"]})')
    axs[0].axis('off')
    fig.colorbar(im0, ax=axs[0])

    axs[1].imshow(post_sample['img'])
    im1 = axs[1].scatter(post_sample['positions']['x_lowres'], post_sample['positions']['y_lowres'], c=post_expr, cmap='hot', s=25, vmin=vmin, vmax=vmax)
    axs[1].set_title(f'POST: {gene}({post_sample["folder"]})')
    axs[1].axis('off')
    fig.colorbar(im1, ax=axs[1])
     # Save as PDF
    plt.savefig("X:/projects/bsc83/Projects/GTEx_v8/Ole/paper/Spatial/Spatial_MOFA_factor4.pdf", format="pdf", bbox_inches="tight")

    # You can also save as PNG with DPI control
    plt.savefig("X:/projects/bsc83/Projects/GTEx_v8/Ole/paper/Spatial/Spatial_MOFA_factor4.png", format="png", dpi=300, bbox_inches="tight")
    plt.suptitle(f'{gene} expression (shared color scale)')
    plt.show()

top_gene = "VEGFA"
peri_samples = [sample for sample in all_samples if sample['Status'] == 'peri']
post_samples = [sample for sample in all_samples if sample['Status'] == 'post']
compare_spatial_with_shared_scale(peri_samples[0], post_samples[0], top_gene, genes)