#!/usr/bin/env python
# coding: utf-8

import cellrank as cr
import scvelo as scv
import scanpy as sc
import numpy as np
import sys


adata_path = sys.argv[1]
path = sys.argv[2]

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2

scv.settings.figdir = path
sc.settings.figdir = path
################
adata = sc.read_h5ad(adata_path)


# # prepare layers
# hvg annotation
sc.pp.highly_variable_genes(adata)
print(f"This detected {np.sum(adata.var['highly_variable'])} highly variable genes. ")


# use scVelo's `moments` function for imputation - note that hack we're using here:
# we're copying our `.X` matrix into the layers because that's where `scv.tl.moments`
# expects to find counts for imputation
adata.layers["spliced"] = adata.X
adata.layers["unspliced"] = adata.X
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


# # Initialize the CytoTRACE kernel

from cellrank.tl.kernels import CytoTRACEKernel

ctk = CytoTRACEKernel(adata)


sc.pl.violin(adata, keys=["ct_pseudotime"], groupby='Ident', rotation=90,
    save='_cytotrace_pseudotime_violin.png')


# # Compute & visualize a transition matrix

ctk.compute_transition_matrix(threshold_scheme="soft", nu=0.5)
ctk.compute_projection(basis="umap")



scv.pl.velocity_embedding_stream(
    adata, color="ct_pseudotime", vkey="T_fwd", basis="umap", legend_loc="right",
     color_map="gnuplot2",
    save='cytotrace_embedding_stream_pseudotime.png'
)


scv.pl.velocity_embedding_stream(
    adata, color="Ident", vkey="T_fwd", basis="umap", legend_loc="right",
    save='cytotrace_embedding_stream_Ident.png')
