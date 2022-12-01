import os
import sys
import logging
import click
import scanpy as sc
import pandas as pd
import xarray as xr
from numpy import ones
from pandas_plink import read_plink1_bin
from numpy.linalg import cholesky
from limix.qc import quantile_gaussianize

from cellregmap import run_interaction, __version__

# use logging to print statements, display at info level
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


@click.command()
# @click.version_option(version=__version__)
@click.option('--chrom', required=True, help='More info here')
@click.option('--gene-name', required=True)
@click.option('--sample-mapping-file', required=True)
@click.option('--genotype-file', required=True)
@click.option('--phenotype-file', required=True)
@click.option('--context-file', required=True)
@click.option('--kinship-file', required=True)
@click.option('--feature-variant-file', required=True)
@click.option(
    '--output-folder', required=False, default=''
)  # by default current directory, where you are running your script from
@click.option('--n-contexts', required=False, type=int)
def main(
    chrom: str,
    gene_name: str,
    sample_mapping_file: str,
    genotype_file: str,
    phenotype_file: str,
    context_file: str,
    kinship_file: str,
    feature_variant_file: str,
    output_folder: str,
    n_contexts: int = 10,
):

    ######################################
    ###### sample mapping file (SMF) #####
    ######################################

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

    ## extract unique individuals
    donors0 = sample_mapping['genotype_individual_id'].unique()
    donors0.sort()
    logging.info(f'Number of unique donors: {len(donors0)}')

    ######################################################
    ###### check if gene output file already exists ######
    ######################################################

    outfilename = os.path.join(output_folder, f'{gene_name}.csv')

    if os.path.exists(outfilename):
        logging.info('File already exists, exiting')
        sys.exit()

    ######################################
    ############ kinship file ############
    ######################################

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

    ## subset to relevant donors
    K = K.sel(sample_0=donors, sample_1=donors)
    assert all(K.sample_0 == donors)
    assert all(K.sample_1 == donors)

    ## and decompose such as K = hK @ hK.T (using Cholesky decomposition)
    hK = cholesky(K.values)
    hK = xr.DataArray(hK, dims=['sample', 'col'], coords={'sample': K.sample_0.values})
    assert all(hK.sample.values == K.sample_0.values)

    del K  # delete K to free up memory
    logging.info(
       f'Sample mapping number of rows BEFORE intersection: {sample_mapping.shape[0]}'
    )
    ## subsample sample mapping file to donors in the kinship matrix
    sample_mapping = sample_mapping[
        sample_mapping['genotype_individual_id'].isin(donors)
    ]
    logging.info(
        f'Sample mapping number of rows AFTER intersection: {sample_mapping.shape[0]}'
    )

    ## use sel from xarray to expand hK (using the sample mapping file)
    hK_expanded = hK.sel(sample=sample_mapping['genotype_individual_id'].values)
    assert all(
        hK_expanded.sample.values == sample_mapping['genotype_individual_id'].values
    )

    ######################################
    ############ genotype file ###########
    ######################################

    ## read in genotype file (plink format)
    G = read_plink1_bin(genotype_file)

    ## select relavant SNPs based on feature variant filter file
    fvf = pd.read_csv(feature_variant_file, index_col=0)
    leads = fvf[fvf['feature'] == gene_name]['snp_id'].unique()
    G_sel = G[:, G['snp'].isin(leads)]

    # expand out genotypes from cells to donors (and select relevant donors in the same step)
    G_expanded = G_sel.sel(sample=sample_mapping['individual_long'].values)
    # assert all(hK_expanded.sample.values == G_expanded.sample.values)

    # delete large files to free up memory
    del G
    del G_sel

    ######################################
    ############ context file ############
    ######################################

    # cells by contexts
    C = pd.read_pickle(context_file)
    C = xr.DataArray(
        C.values,
        dims=['cell', 'pc'],
        coords={'cell': C.index.values, 'pc': C.columns.values},
    )
    C = C.sel(cell=sample_mapping['phenotype_sample_id'].values)
    assert all(C.cell.values == sample_mapping['phenotype_sample_id'].values)

    ######################################
    ########### phenotype file ###########
    ######################################

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

    ######################################
    ########### Prepare model ############
    ######################################

    n_cells = phenotype.shape[1]
    W = ones((n_cells, 1))  # just intercept as covariates

    # select gene
    y = phenotype.sel(trait=gene_name)

    y = quantile_gaussianize(y)
    y = y.values.reshape(y.shape[0], 1)

    del phenotype  # delete to free up memory

    GG = G_expanded.values

    ##################################
    ########### Run model ############
    ##################################

    logging.info(f'Running for gene {gene_name}')

    pvals = run_interaction(
        y=y, W=W, E=C.values[:, 0:n_contexts], G=GG, hK=hK_expanded
    )[0]

    pv = pd.DataFrame(
        {
            'chrom': G_expanded.chrom.values,
            'pv': pvals,
            'variant': G_expanded.snp.values,
        }
    )
    pv.to_csv(outfilename)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
