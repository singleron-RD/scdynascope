#!/usr/bin/env python

import argparse
import scanpy as sc
import numpy as np

import utils
from __init__ import ASSAY

class Labeled_summary:
    """
    Features
    - Labeled summary of per gene and per cell. 
    """

    def __init__(self, args):
        self.args = args

        # input files
        self.outdir = args.outdir
        self.sample = args.sample
        self.dir_labeled = args.labeled_matrix
        self.dir_total = args.filtered_matrix
        self.min_cells = args.min_cells
        self.min_genes = args.min_genes
        
        # output json 
        self.report_dict = {}

    
    def run(self):
        adata = self.read_mtx()
        self.tor_cal(adata)
        utils.write_multiqc(self.report_dict, self.sample, ASSAY, "tor")

    def read_mtx(self):
        bdata = sc.read_10x_mtx(self.dir_labeled, var_names='gene_ids')
        adata = sc.read_10x_mtx(self.dir_total, var_names='gene_ids')
        adata.layers['total']=adata.X.copy()
        adata.layers['labeled']=bdata.X.copy()
        return self.dropna_gene(adata)

    def dropna_gene(self, adata):
        gene_sum = adata.X.sum(axis=0)
        nonzero_genes = gene_sum.nonzero()[1]
        adata = adata[:, nonzero_genes]
        return adata

    def tor_cal(self, adata):
        adata = self.dropna_gene(adata)
        cell_ntr = adata.layers['labeled'].sum(axis=1) / adata.layers['total'].sum(axis=1)
        cell_ntr = np.array(cell_ntr).flatten().tolist()
        gene_num = adata.layers['labeled'].getnnz(axis=1)
        gene_ntr = adata.layers['labeled'].sum(axis=0) / adata.layers['total'].sum(axis=0)
        gene_ntr = np.array(gene_ntr).flatten().tolist()
        cell_num = adata.layers['labeled'].getnnz(axis=0)
        #adata.obs['TOR'] = cell_ntr
        #adata.var['TOR'] = gene_ntr.T
        self.report_dict['Cells'] = [elem2 for elem1, elem2 in zip(gene_num, cell_ntr) if elem1 > self.min_genes]
        self.report_dict['Genes'] = [elem2 for elem1, elem2 in zip(cell_num, gene_ntr) if elem1 > self.min_cells]


def get_opts_labeled_summary():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--outdir', help='Output diretory.', required=True)
    parser.add_argument('--sample', help='Sample name.', required=True)
    parser.add_argument("--labeled_matrix", help='labeled matrix dir', required=True)
    parser.add_argument("--filtered_matrix", help='filtered matrix dir', required=True)
    parser.add_argument("--min_cells", help='min cell number to call a gene\'s tor', 
                        type=int, default=10, required=False)
    parser.add_argument("--min_genes", help='min gene number to call a cell\'s tor', 
                        type=int, default=10, required=False)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args=get_opts_labeled_summary()
    Labeled_summary(args).run()