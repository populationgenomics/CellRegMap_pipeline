import click
# import os
# import sys
# import scanpy as sc
import pandas as pd
import xarray as xr
# from numpy import ones
# from pandas_plink import read_plink1_bin
from numpy.linalg import cholesky
# import time
# from limix.qc import quantile_gaussianize

# from cellregmap import run_interaction 

@click.Click()
@click.parameter('--chrom', flag=True)
@click.parameter('--gene-index', flag=True)
@click.parameter('--sample-mapping-file', required=True)
@click.parameter('--genotype-file', required=True)
@click.parameter('--phenotype-file', required=True)
@click.parameter('--context-file', required=True)
@click.parameter('--kinship-file', required=True)
@click.parameter('--feature-variant-file', required=True)
@click.parameter('--covariate-file', required=False)

def main(chrom, gene_index, sample_mapping_file, genotype_file, phenotype_file, context_file, kinship_file, feature_variant_file, covariate_file):
    
    ######################################
    ###### sample mapping file (SMF) #####
    ######################################

    ## this file will map cells to donors 
    sample_mapping = pd.read_csv(sample_mapping_file, dtype={"individual_long": str, "genotype_individual_id": str, "phenotype_sample_id": str}, index_col=0)

    ## extract unique individuals
    donors0 = sample_mapping["genotype_individual_id"].unique()
    donors0.sort()
    print("Number of unique donors: {}".format(len(donors0)))

    ######################################################
    ###### check if gene output file already exists ######
    ######################################################

    fvf = pd.read_csv(feature_variant_file, index_col = 0)
    genes = fvf[fvf['chrom']==int(chrom)]['feature'].unique()
    gene_name = genes[gene_index]
    
    outfilename = f"{gene_name}.tsv"

    if os.path.exists(outfilename):
        print("File already exists, exiting")
        sys.exit()

    ######################################
    ############ kinship file ############
    ######################################

    ## read in GRM (genotype relationship matrix; kinship matrix)
    K = pd.read_csv(kinship_file, index_col=0)
    K.index = K.index.astype('str')
    assert all(K.columns == K.index) #symmetric matrix, donors x donors

    K = xr.DataArray(K.values, dims=["sample_0", "sample_1"], coords={"sample_0": K.columns, "sample_1": K.index})
    K = K.sortby("sample_0").sortby("sample_1")
    donors = sorted(set(list(K.sample_0.values)).intersection(donors0))
    print("Number of donors after kinship intersection: {}".format(len(donors)))

    ## subset to relevant donors
    K = K.sel(sample_0=donors, sample_1=donors)
    assert all(K.sample_0 == donors)
    assert all(K.sample_1 == donors)

    ## and decompose such as K = hK @ hK.T (using Cholesky decomposition)
    hK = cholesky(K.values)
    hK = xr.DataArray(hK, dims=["sample", "col"], coords={"sample": K.sample_0.values})
    assert all(hK.sample.values == K.sample_0.values)

    del K
    print("Sample mapping number of rows BEFORE intersection: {}".format(sample_mapping.shape[0]))
    ## subsample sample mapping file to donors in the kinship matrix
    sample_mapping = sample_mapping[sample_mapping["genotype_individual_id"].isin(donors)]
    print("Sample mapping number of rows AFTER intersection: {}".format(sample_mapping.shape[0]))

    ## use sel from xarray to expand hK (using the sample mapping file)
    hK_expanded = hK.sel(sample=sample_mapping["genotype_individual_id"].values)
    assert all(hK_expanded.sample.values == sample_mapping["genotype_individual_id"].values)

"""

Task A:
    cwd: Gene A:
    /Users/anna/cromwell-workfig-dir/workflow-id/workflonwame/taskname/shard-1/execution/

    cwd: Gene X
    /Users/anna/cromwell-workfig-dir/workflow-id/workflonwame/taskname/shard-25/execution/


Summarise Task:
    cwd: TaskName/execution
    inputs:
        .../TaskName/inputs/12jasd123/gene1.txt
        .../TaskName/inputs/12u6i98ca/gene2.txt
"""

# os.path.splitext(os.path.basename("/big/long/myfilename.txt"))[0] # myfilename

