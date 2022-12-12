# pylint: disable=missing-function-docstring,no-value-for-parameter,too-many-locals

import os
import sys
import logging

import click
import numpy as np
import pandas as pd
import scanpy as sc
import xarray as xr

from limix.qc import quantile_gaussianize
from numpy import ones
from numpy.linalg import cholesky
from pandas_plink import read_plink1_bin

from cellregmap import estimate_betas

# use logging to print statements, display at info level
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


@click.command()
@click.option('--gene-name', required=True)
@click.option('--sample-mapping-file', required=True)
@click.option('--genotype-file', required=True)
@click.option('--phenotype-file', required=True)
@click.option('--context-file', required=True)
@click.option('--kinship-file', required=True)
@click.option('--beta-feature-variant-file', required=True)
@click.option(
    '--output-folder', required=False, default=''
)  # by default current directory, where you are running your script from
@click.option('--n-contexts', required=False, type=int)
@click.option('--maf-file', required=False)
def main(
    gene_name: str,
    sample_mapping_file: str,
    genotype_file: str,
    phenotype_file: str,
    context_file: str,
    kinship_file: str,
    beta_feature_variant_file: str,
    output_folder: str,
    maf_file: str = None,
    n_contexts: int = 10,
):

    # region SAMPLE_MAPPING_FILE

    ## this file will map cells to donors
    sample_mapping = pd.read_csv(
        sample_mapping_file,
        dtype={
            'individual_long': str,
            'genotype_individual_id': str,
            'phenotype_sample_id': str,
        },
        index_col=0,
    )
    # extract unique individuals
    donors0 = sample_mapping['genotype_individual_id'].unique()
    donors0.sort()
    logging.info(f'Number of unique donors: {len(donors0)}')

    # endregion SAMPLE_MAPPING_FILE

    # check if gene output file already exists

    outfilename = f'{output_folder}{gene_name}'
    outfilename_betaGxC = outfilename + '_betaGxC.csv'

    if os.path.exists(outfilename_betaGxC):
        logging.info('File already exists, exiting')
        sys.exit()

    # region KINSHIP_FILE

    ## read in GRM (genotype relationship matrix; kinship matrix)
    K = pd.read_csv(kinship_file, index_col=0)
    K.index = K.index.astype('str')
    assert all(K.columns == K.index)  # symmetric matrix, donors x donors

    K = xr.DataArray(
        K.values,
        dims=['sample_0', 'sample_1'],
        coords={'sample_0': K.columns, 'sample_1': K.index},
    )
    K = K.sortby('sample_0').sortby('sample_1')
    donors = sorted(set(list(K.sample_0.values)).intersection(donors0))
    logging.info(f'Number of donors after kinship intersection: {len(donors)}')

    # subset to relevant donors
    K = K.sel(sample_0=donors, sample_1=donors)
    assert all(K.sample_0 == donors)
    assert all(K.sample_1 == donors)

    # and decompose such as K = hK @ hK.T (using Cholesky decomposition)
    hK = cholesky(K.values)
    hK = xr.DataArray(hK, dims=['sample', 'col'], coords={'sample': K.sample_0.values})
    assert all(hK.sample.values == K.sample_0.values)

    del K  # delete to free up memory
    logging.info(
        f'Sample mapping number of rows BEFORE intersection: {sample_mapping.shape[0]}'
    )
    # subsample sample mapping file to donors in the kinship matrix
    sample_mapping = sample_mapping[
        sample_mapping['genotype_individual_id'].isin(donors)
    ]
    logging.info(
        f'Sample mapping number of rows AFTER intersection: {sample_mapping.shape[0]}'
    )

    # use sel from xarray to expand hK (using the sample mapping file)
    hK_expanded = hK.sel(sample=sample_mapping['genotype_individual_id'].values)
    assert all(
        hK_expanded.sample.values == sample_mapping['genotype_individual_id'].values
    )

    # region GENOTYPE_FILE

    # read in genotype file (plink format)
    G = read_plink1_bin(genotype_file)

    # select relevant SNPs based on feature variant filter file
    # these are variants that were identified to have GxC effects
    # we want to estimate single-cell effect sizes (betas) for
    fvf = pd.read_csv(beta_feature_variant_file, index_col=0)
    leads = fvf[fvf['gene'] == gene_name]['snp_id'].unique()
    G_sel = G[:, G['snp'].isin(leads)]

    # expand out genotypes from cells to donors (and select relevant donors in the same step)
    G_expanded = G_sel.sel(sample=sample_mapping['individual_long'].values)

    # delete large files to free up memory
    del G
    del G_sel

    # endregion GENOTYPE_FILE

    # region CONTEXT_FILE

    # cells by contexts
    C = pd.read_pickle(context_file)
    C = xr.DataArray(
        C.values,
        dims=['cell', 'pc'],
        coords={'cell': C.index.values, 'pc': C.columns.values},
    )
    C = C.sel(cell=sample_mapping['phenotype_sample_id'].values)
    assert all(C.cell.values == sample_mapping['phenotype_sample_id'].values)

    # endregion CONTEXT_FILE

    # region PHENOTYPE_FILE

    # open anndata
    adata = sc.read(phenotype_file)
    # sparse to dense
    mat = adata.raw.X.todense()
    # make pandas dataframe
    mat_df = pd.DataFrame(
        data=mat.T, index=adata.raw.var.index, columns=adata.obs.index
    )
    # turn into xr array
    phenotype = xr.DataArray(
        mat_df.values,
        dims=['trait', 'cell'],
        coords={'trait': mat_df.index.values, 'cell': mat_df.columns.values},
    )
    phenotype = phenotype.sel(cell=sample_mapping['phenotype_sample_id'].values)

    # delete large files to free up memory
    del mat
    del mat_df

    # endregion PHENOTYPE_FILE

    # region PREPARE_MODEL

    n_cells = phenotype.shape[1]
    W = ones((n_cells, 1))  # just intercept as covariates

    # select gene
    y = phenotype.sel(trait=gene_name)

    y = quantile_gaussianize(y)
    y = y.values.reshape(y.shape[0], 1)

    cells = phenotype['cell'].values
    del phenotype  # delete to free up memory

    GG = G_expanded.values
    snps = G_expanded['snp'].values

    # get minor allele frequencies (MAFs)
    # ideally these are precomputed, if not they can be calculated

    mafs = None
    if maf_file:
        df_maf = pd.read_csv(maf_file, sep='\t')

        mafs = np.array([])
        for snp in snps:
            mafs = np.append(mafs, df_maf[df_maf['SNP'] == snp]['MAF'].values)

    # endregion PREPARE_MODEL

    # region RUN_MODEL

    logging.info(f'Running for gene {gene_name}')

    betas = estimate_betas(
        y=y, W=W, E=C.values[:, 0:n_contexts], G=GG, maf=mafs, hK=hK_expanded
    )
    beta_G = betas[0]
    beta_GxC = betas[1][0]

    beta_G_df = pd.DataFrame(
        {
            'chrom': G_expanded.chrom.values,
            'betaG': beta_G,
            'variant': G_expanded.snp.values,
        }
    )

    beta_G_df.to_csv(outfilename + '_betaG.csv')

    beta_GxC_df = pd.DataFrame(data=beta_GxC, columns=snps, index=cells)
    beta_GxC_df.to_csv(outfilename_betaGxC)

    # endregion RUN_MODEL


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
