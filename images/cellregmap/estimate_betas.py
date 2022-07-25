import click
import os
import sys
import numpy as np
import scanpy as sc
import pandas as pd
import xarray as xr
from numpy import ones
from pandas_plink import read_plink1_bin
from numpy.linalg import cholesky
from limix.qc import quantile_gaussianize

from cellregmap import estimate_betas

@click.command()
@click.option('--chrom', required=True)
@click.option('--gene-name', required=True)
@click.option('--sample-mapping-file', required=True)
@click.option('--genotype-file', required=True)
@click.option('--phenotype-file', required=True)
@click.option('--context-file', required=True)
@click.option('--kinship-file', required=True)
@click.option('--beta-feature-variant-file', required=True)
@click.option(
    "--output-folder", required=False, default=""
)  # by default current directory, where you are running your script from
@click.option('--n-contexts', required=False, type=int)
@click.option('--maf-file', required=False)

def main(
    chrom, 
    gene_name, 
    sample_mapping_file, 
    genotype_file, 
    phenotype_file, 
    context_file, 
    kinship_file, 
    beta_feature_variant_file, 
    output_folder, 
    maf_file,
    n_contexts=10):
    
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

    outfilename = f"{output_folder}{gene_name}"
    outfilename_betaGxC = outfilename+"_betaGxC.csv"

    if os.path.exists(outfilename_betaGxC):
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

    ######################################
    ############ genotype file ###########
    ######################################

    ## read in genotype file (plink format)
    G = read_plink1_bin(genotype_file)

    ## select relavant SNPs based on feature variant filter file 
    fvf = pd.read_csv(beta_feature_variant_file, index_col = 0)
    leads = fvf[fvf['gene']==gene_name]['snp_id'].unique()
    G_sel = G[:,G['snp'].isin(leads)]

    # expand out genotypes from cells to donors (and select relevant donors in the same step)
    G_expanded = G_sel.sel(sample=sample_mapping["individual_long"].values)
    #assert all(hK_expanded.sample.values == G_expanded.sample.values)

    del G

    ######################################
    ############ context file ############
    ######################################

    # cells by contexts
    C = pd.read_pickle(context_file)
    C = xr.DataArray(C.values, dims=["cell", "pc"], coords={"cell": C.index.values, "pc": C.columns.values})
    C = C.sel(cell=sample_mapping["phenotype_sample_id"].values)
    assert all(C.cell.values == sample_mapping["phenotype_sample_id"].values)

    ######################################
    ########### phenotype file ###########
    ######################################

    # open anndata 
    adata = sc.read(phenotype_file)
    # sparse to dense
    mat = adata.raw.X.todense()
    # make pandas dataframe
    mat_df = pd.DataFrame(data=mat.T, index=adata.raw.var.index, columns=adata.obs.index)
    # turn into xr array
    phenotype = xr.DataArray(mat_df.values, dims=["trait", "cell"], coords={"trait": mat_df.index.values, "cell": mat_df.columns.values})
    phenotype = phenotype.sel(cell=sample_mapping["phenotype_sample_id"].values)

    del mat
    del mat_df

    ######################################
    ########### Prepare model ############
    ######################################

    n_cells = phenotype.shape[1]
    W = ones((n_cells, 1)) # just intercept as covariates

    # select gene
    y = phenotype.sel(trait=gene_name)

    y = quantile_gaussianize(y)
    y = y.values.reshape(y.shape[0],1)

    cells = phenotype["cell"].values
    del phenotype

    GG = G_expanded.values
    snps = G_expanded["snp"].values

    del G_sel

    # get MAF
    df_maf = pd.read_csv(maf_file, sep="\t")

    mafs = np.array([])
    for snp in snps:
        mafs = np.append(mafs, df_maf[df_maf["SNP"] == snp]["MAF"].values)

    ##################################
    ########### Run model ############
    ##################################

    print("Running for gene {}".format(gene_name))

    betas = estimate_betas(y=y, W=W, E=C.values[:,0:10], G=GG, maf=mafs, hK=hK_expanded)
    beta_G = betas[0]
    beta_GxC = betas[1][0]

    beta_G_df = pd.DataFrame({"chrom":G_expanded.chrom.values,
                "betaG":beta_G,
                "variant":G_expanded.snp.values})

    beta_G_df.to_csv(outfilename+"_betaG.csv")

    beta_GxC_df = pd.DataFrame(data = beta_GxC, columns=snps, index=cells)
    beta_GxC_df.to_csv(outfilename_betaGxC)


if __name__ == '__main__':
    main()  # pylint disable=no-value-for-argument
