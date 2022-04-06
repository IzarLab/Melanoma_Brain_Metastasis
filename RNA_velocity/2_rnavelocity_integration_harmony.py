#!/usr/bin/env python3

# Merging loom files, integration using harmony and RNAvelocity pipeline
# Authors: Jana Biermann, PhD; Somnath Tagore, PhD

# Import packages
import scvelo as scv
import matplotlib
import matplotlib.pyplot as pl
import pandas as pd
import sys
import os
import scanpy as sc
import scanpy.external as sce

print("\nArguments passed:",sys.argv[1])

# Provide MDM_MBM or MDM_MPM
label = sys.argv[1]

current_path = "data/RNA_velocity/"
if not os.path.exists(current_path):
    os.makedirs(current_path)


def read_loom(sample):
    # Read-in loom files, add info from Seurat object, and subset
    print('Reading in',sample)

    # read loom file
    adata = scv.read("".join(["data/", sample, "/", sample, ".loom"]), sparse=True, cache=True)
    adata.var_names_make_unique()

    # Add cell type info to adata
    anno = pd.read_csv("".join(["data/", sample, "/", sample, "_ct_myeloid.csv"]))
    tmp = anno.iloc[:, 0].str.split('-', expand=True)
    tmp['CellID'] = "".join([sample, ":"]) + tmp.iloc[:, 0]
    anno.index = tmp['CellID']
    adata = adata[tmp['CellID'], :]
    adata.obs['cell_type_fine'] = anno['cell_type_fine']
    adata.obs['M1_ada1'] = anno['M1_ada1']
    adata.obs['M2_ada1'] = anno['M2_ada1']
    adata.obs['patient'] = sample

    # Subset to MDMs
    adata = adata[adata.obs['cell_type_fine'].isin(["MDM FTL+", "MDM"])]
    return adata

if label in ['MDM_MBM']:
    MBM05_sn = read_loom('MBM05_sn',sys.argv[1])
    MBM06_sn = read_loom('MBM06_sn',sys.argv[1])
    MBM08_sn = read_loom('MBM08_sn',sys.argv[1])
    MBM09_sn = read_loom('MBM09_sn',sys.argv[1])
    MBM10_sn = read_loom('MBM10_sn',sys.argv[1])
    MBM11_sn = read_loom('MBM11_sn',sys.argv[1])
    MBM14_sn = read_loom('MBM14_sn',sys.argv[1])
    MBM15_sn = read_loom('MBM15_sn',sys.argv[1])
    MBM17_sn = read_loom('MBM17_sn',sys.argv[1])
    MBM20_sn = read_loom('MBM20_sn',sys.argv[1])
    MBM21_sn = read_loom('MBM21_sn',sys.argv[1])
    # Merge into an already existing AnnData object
    adata = MBM05_sn.concatenate(MBM06_sn, MBM08_sn, MBM09_sn, MBM10_sn, MBM11_sn, MBM14_sn,
                                 MBM15_sn, MBM17_sn, MBM20_sn, MBM21_sn)

if label in ['MDM_MPM']:
    MPM01_sn = read_loom('MPM01_sn',sys.argv[1])
    MPM02_sn = read_loom('MPM02_sn',sys.argv[1])
    MPM04_sn = read_loom('MPM04_sn',sys.argv[1])
    MPM05_sn = read_loom('MPM05_sn',sys.argv[1])
    MPM06_sn = read_loom('MPM06_sn',sys.argv[1])
    MPM07_sn = read_loom('MPM07_sn',sys.argv[1])
    MPM08_sn = read_loom('MPM08_sn',sys.argv[1])
    MPM09_sn = read_loom('MPM09_sn',sys.argv[1])
    MPM10_sn = read_loom('MPM10_sn',sys.argv[1])
    MPM11_sn = read_loom('MPM11_sn',sys.argv[1])
    # Merge into an already existing AnnData object
    adata = MPM01_sn.concatenate(MPM02_sn, MPM04_sn, MPM05_sn,
                                 MPM06_sn, MPM07_sn, MPM08_sn, MPM09_sn, MPM10_sn, MPM11_sn)

# Save merged object
adata.write("".join([current_path,"Merged_",sys.argv[1],".h5ad"]))

# Read-in merged object
#adata = sc.read_h5ad("".join([current_path,"Merged_",sys.argv[1],".h5ad"]))

# Run harmony
sc.tl.pca(adata)
sce.pp.harmony_integrate(adata, 'patient')
print('X_pca_harmony' in adata.obsm)

# Replace PCA with harmony
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']

# Neighbors and louvain clustering
sc.pp.neighbors(adata)
scv.tl.louvain(adata, resolution=0.7)

# UMAP
scv.tl.umap(adata)

# Run RNAvelocity
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata, n_jobs=16)

# Recover dynamics
scv.tl.recover_dynamics(adata, n_jobs=16)

# Latent time
scv.tl.latent_time(adata)

# Cell fate
scv.tl.terminal_states(adata)

# Speed and coherence
scv.tl.velocity_confidence(adata)

# Pseudotime
scv.tl.velocity_pseudotime(adata)

# PAGA graph abstraction for trajectory inference.
# It provides a graph-like map of the data topology with weighted edges corresponding to the
# connectivity between two clusters. Here, PAGA is extended by velocity-inferred directionality.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='cell_type_fine')

# Save h5ad object
adata.write("".join([current_path,"Harmony_",sys.argv[1],"_velocity.h5ad"]))