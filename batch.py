#!/usr/bin/env python3
# pylint: disable=import-error,import-outside-toplevel,missing-function-docstring,no-value-for-parameter,too-many-arguments,too-many-locals,wrong-import-position

"""
Hail Batch workflow for the rare-variant association analysis, including:

- get relevant variants around a gene and export genotypes as plink files
- generate input files for association tests
- run association tests
"""

# import python modules
import os
import re

import click
import logging

from cloudpathlib import AnyPath
from cpg_utils import to_path
from cpg_utils.hail_batch import (
    copy_common_env,
    dataset_path,
    get_config,
    init_batch,
    output_path,
    remote_tmpdir,
)

import numpy as np
import pandas as pd
import xarray as xr

from limix.qc import quantile_gaussianize
from scipy.stats import shapiro

import hail as hl
import hailtop.batch as hb

from cellregmap import (  # figure out how to import this from github
    run_gene_set_association,
    run_burden_association,
    omnibus_set_association,
)

# from ._utils import qv_estimate as qvalue

# use logging to print statements, display at info level
logging.basicConfig(
    format='%(levelname)s:%(message)s', level=logging.INFO
)  # consider removing some print statements?

DEFAULT_JOINT_CALL_MT = dataset_path('mt/v7.mt')
DEFAULT_ANNOTATION_HT = dataset_path(
    'tob_wgs_vep/104/vep104.3_GRCh38.ht'
)  # atm VEP only


CELLREGMAP_IMAGE = get_config()['workflow'][
    'driver_image'
]  # australia-southeast1-docker.pkg.dev/cpg-common/images/cellregmap:dev

# region SUBSET_VARIANTS


def filter_variants(
    mt_path: str,  # 'mt/v7.mt'
    samples: list[str],
    output_mt_path: str,  # 'tob_wgs_rv/densified_rv_only.mt'
):
    """Subset hail matrix table

    Input:
    joint call hail matrix table
    set of samples for which we have scRNA-seq data

    Output:
    subset hail matrix table, containing only variants that:
    1) are not ref-only, 2) biallelic, 3) meet QC filters, 4) are rare (MAF<5%)
    """
    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    # subset to relevant samples (samples we have scRNA-seq data for)
    mt = mt.filter_cols(hl.set(samples).contains(mt.s))

    # densify
    mt = hl.experimental.densify(mt)

    # filter out low quality variants and consider biallelic variants only (no multi-allelic, no ref-only)
    mt = mt.filter_rows(  # check these filters!
        (hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)  # QC
        & (hl.len(mt.alleles) == 2)  # remove hom-ref
        & (mt.n_unsplit_alleles == 2)  # biallelic
        & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNVs
    )

    # filter rare variants only (MAF < 5%)
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] < 0.05) & (mt.variant_qc.AF[1] > 0)
        | (mt.variant_qc.AF[1] > 0.95) & (mt.variant_qc.AF[1] < 1)
    )
    mt.write(output_mt_path, overwrite=True)
    logging.info(
        f'Number of rare (freq<5%) and QCed biallelic variants: {mt.count()[0]}'
    )


# endregion SUBSET_VARIANTS


# region GET_GENE_SPECIFIC_VARIANTS


def get_promoter_variants(
    mt_path: str,  # ouput path from function above
    ht_path: str,
    gene_details: dict[str, str],  # output of make_gene_loc_dict
    window_size: int,
    plink_file: str,  # 'tob_wgs_rv/pseudobulk_rv_association/plink_files/GENE'
):
    """Subset hail matrix table

    Input:
    mt_path: path to already subsetted hail matrix table
    ht_path: path to VEP HT
    gene_details: dict of info for current gene
    window_size: int, size of flanking region around genes
    plink_file: str, file prefix for writing plink data

    Output:
    For retained variants, that are: 1) in promoter regions and
    2) within 50kb up or down-stream of the gene body (or in the gene body itself)
    (on top of all filters done above)

    returns nothing
    """

    # read hail matrix table object (pre-filtered)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    gene_name = gene_details['gene_name']

    # get relevant chromosome
    chrom = gene_details['chr']

    # subset to window
    # get gene body position (start and end) and build interval
    left_boundary = max(1, int(gene_details['start']) - window_size)
    right_boundary = min(
        int(gene_details['end']) + window_size,
        hl.get_reference('GRCh38').lengths[chrom],
    )
    # get gene-specific genomic interval
    gene_interval = f'{chrom}:{left_boundary}-{right_boundary}'
    logging.info(f'Interval considered: {gene_interval}')  # 'chr22:23219960-23348287'

    # include variants up to {window size} up- and downstream
    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval(gene_interval, reference_genome='GRCh38')]
    )
    mt_path = output_path(f'{gene_name}_in_window.mt', 'tmp')
    mt = mt.checkpoint(
        mt_path, overwrite=True
    )  # add checkpoint to avoid repeat evaluation
    logging.info(f'Number of variants within interval: {mt.count()[0]}')

    # annotate using VEP
    vep_ht = hl.read_table(ht_path)
    mt = mt.annotate_rows(vep=vep_ht[mt.row_key].vep)

    # filter variants found to be in promoter regions
    mt = mt.filter_rows(
        mt.vep.regulatory_feature_consequences['biotype'].contains('promoter')
    )
    promoter_path = output_path(f'{gene_name}promoter_variants.mt', 'tmp')
    mt = mt.checkpoint(promoter_path, overwrite=True)  # checkpoint
    logging.info(
        f'Number of rare (freq<5%) QC-passing, biallelic SNPs in promoter regions: {mt.count()[0]}'
    )

    # export this as a Hail table for downstream analysis
    ht_path = output_path(f'summary_hts/{gene_name}_rare_promoter_summary.ht')
    ht = mt.rows()
    ht.write(ht_path, overwrite=True)

    # export MT object to PLINK (promoter variants)
    # pylint: disable=import-outside-toplevel
    from hail.methods import export_plink

    export_plink(mt, plink_file, ind_id=mt.s)


# endregion GET_GENE_SPECIFIC_VARIANTS

# region PREPARE_INPUT_FILES


def prepare_input_files(
    gene_name: str,
    cell_type: str,
    genotype_file_bed: str,
    genotype_file_bim: str,
    genotype_file_fam: str,
    phenotype_file: str,
    kinship_file: str,
    sample_mapping_file: str,
):
    """Prepare association test input files

    Input:
    genotype: plink files
    phenotype: gene expression

    Output:
    genotype matrix
    phenotype vector
    """
    from pandas_plink import read_plink1_bin

    expression_filename = AnyPath(output_path(f'{gene_name}_{cell_type}.csv'))
    genotype_filename = AnyPath(output_path(f'{gene_name}_rare_regulatory.csv'))
    kinship_filename = AnyPath(output_path(f'{gene_name}_kinship_common_samples.csv'))

    # read in phenotype file (tsv)
    phenotype = pd.read_csv(phenotype_file, sep='\t', index_col=0)

    phenotype = xr.DataArray(
        phenotype.values,
        dims=['sample', 'gene'],
        coords={'sample': phenotype.index.values, 'gene': phenotype.columns.values},
    )

    # read in genotype file (plink format)
    to_path(genotype_file_bed).copy('temp.bed')  # bed
    to_path(genotype_file_bim).copy('temp.bim')  # bim
    to_path(genotype_file_fam).copy('temp.fam')  # fam
    geno = read_plink1_bin('temp.bed')

    if kinship_file is not None:
        # read in GRM (genotype relationship matrix; kinship matrix)
        kinship = pd.read_csv(kinship_file, index_col=0)
        kinship.index = kinship.index.astype('str')
        assert all(
            kinship.columns == kinship.index
        )  # symmetric matrix, donors x donors
        kinship = xr.DataArray(
            kinship.values,
            dims=['sample_0', 'sample_1'],
            coords={'sample_0': kinship.columns, 'sample_1': kinship.index},
        )
        kinship = kinship.sortby('sample_0').sortby('sample_1')

    # this file will map different IDs (and OneK1K ID to CPG ID)
    sample_mapping = pd.read_csv(dataset_path(sample_mapping_file), sep='\t')

    # ensure samples are the same and in the same order across input files
    # samples with expression data
    donors_exprs = set(phenotype.sample.values).intersection(
        set(sample_mapping['OneK1K_ID'].unique())
    )

    logging.info(f'Number of unique donors with expression data: {len(donors_exprs)}')

    # samples with genotype data
    donors_geno = set(geno.sample.values).intersection(
        set(sample_mapping['InternalID'].unique())
    )
    logging.info(f'Number of unique donors with genotype data: {len(donors_geno)}')

    # samples with both (can this be done in one step?)
    sample_mapping1 = sample_mapping.loc[sample_mapping['OneK1K_ID'].isin(donors_exprs)]
    sample_mapping_both = sample_mapping1.loc[
        sample_mapping1['InternalID'].isin(donors_geno)
    ]
    donors_e = sample_mapping_both['OneK1K_ID'].unique()
    donors_g = sample_mapping_both['InternalID'].unique()
    assert len(donors_e) == len(donors_g)

    if kinship_file is not None:
        # samples in kinship
        donors_e_short = [re.sub('.*_', '', donor) for donor in donors_e]
        donors_k = sorted(
            set(list(kinship.sample_0.values)).intersection(donors_e_short)
        )

    logging.info(f'Number of unique common donors: {len(donors_g)}')

    # subset files

    # phenotype
    phenotype = phenotype.sel(sample=donors_e)
    # select gene
    y = phenotype.sel(gene=gene_name)
    y = quantile_gaussianize(y)
    del phenotype  # delete to free up memory
    # make data frame to save as csv
    y_df = pd.DataFrame(
        data=y.values.reshape(y.shape[0], 1), index=y.sample.values, columns=[gene_name]
    )

    # genotype
    geno = geno.sel(sample=donors_g)
    # make data frame to save as csv
    data = geno.values
    geno_df = pd.DataFrame(data, columns=geno.snp.values, index=geno.sample.values)
    geno_df = geno_df.dropna(axis=1)
    # delete large files to free up memory
    del geno

    if kinship_file is not None:
        # kinship
        kinship = kinship.sel(sample_0=donors_k, sample_1=donors_k)
        assert all(kinship.sample_0 == donors_k)
        assert all(kinship.sample_1 == donors_k)
        # make data frame to save as csv
        kinship_df = pd.DataFrame(
            kinship.values, columns=kinship.sample_0, index=kinship.sample_1
        )
        del kinship  # delete kinship to free up memory

    # save files
    with expression_filename.open('w') as ef:
        y_df.to_csv(ef, index=False)

    with genotype_filename.open('w') as gf:
        geno_df.to_csv(gf, index=False)

    if kinship_file is not None:
        with kinship_filename.open('w') as kf:
            kinship_df.to_csv(kf, index=False)
    else:
        kinship_df = None

    return y_df, geno_df, kinship_df


# endregion PREPARE_INPUT_FILES


# region GET_CRM_PVALUES


def get_crm_pvs(pheno, covs, genotypes, contexts=None):
    """
    CellRegMap-RV tests
    * score test (variance)
    * burden test (max, sum, comphet)
    * omnibus (Cauchy) test

    Input:
    input files

    Output:
    list of p-values from the three tests
    """
    pv_norm = shapiro(pheno).pvalue
    pv0 = run_gene_set_association(y=pheno, G=genotypes, W=covs, E=contexts)[0]
    pv1 = run_burden_association(
        y=pheno, G=genotypes, W=covs, E=contexts, mask='mask.max'
    )[0]
    pv2 = run_burden_association(
        y=pheno, G=genotypes, W=covs, E=contexts, mask='mask.sum'
    )[0]
    pv3 = run_burden_association(
        y=pheno, G=genotypes, W=covs, E=contexts, mask='mask.comphet'
    )[0]
    pv4 = omnibus_set_association(np.array([pv0, pv1]))
    pv5 = omnibus_set_association(np.array([pv0, pv2]))
    pv6 = omnibus_set_association(np.array([pv0, pv3]))

    return np.array([pv_norm, pv0, pv1, pv2, pv3, pv4, pv5, pv6])


# endregion GET_CRM_PVALUES

# region RUN_ASSOCIATION


def run_gene_association(
    gene_name: str,  # 'VPREB3'
    prepared_inputs: hb.resource.PythonResult,
    # genotype_mat_path: str,  # 'VPREB3_50K_window/SNVs.csv'
    # phenotype_vec_path: str,  # 'Bnaive/VPREB3_pseudocounts.csv'
    pv_path: str,  # 'Bnaive/VPREB3_pvals.csv'
):
    """Run gene-set association test

    Input:
    input files  (genotype, phenotype)
    * already in matrix / vector form
    * only matching samples, correct irder

    Output:
    table with p-values
    """
    from numpy import eye, ones

    # read the 3 dataframes generated by the previous job
    p_df, g_df, _ = prepared_inputs

    # because the current matrix is counting the copies of the reference allele
    # while we are interested in the alternative allele, flip the genotypes
    genotypes = 2 - g_df

    # get phenptypes
    pheno = p_df.values

    # covs
    covs = ones((genotypes.shape[0], 1))  # intercept of ones as covariates

    # contexts
    contexts = eye(genotypes.shape[0])

    # TODO: kinship

    cols = np.array(
        [
            'P_shapiro',
            'P_CRM_VC',
            'P_CRM_burden_max',
            'P_CRM_burden_sum',
            'P_CRM_burden_comphet',
            'P_CRM_omnibus_max',
            'P_CRM_omnibus_sum',
            'P_CRM_omnibus_comphet',
        ]
    )

    # create p-values data frame
    pvalues = get_crm_pvs(pheno, covs, genotypes, contexts)
    pv_df = pd.DataFrame(
        data=pvalues.reshape(pvalues.shape[0], 1).T,
        columns=cols,
        index=[gene_name],
    )

    pv_filename = AnyPath(output_path(pv_path))
    # print(pv_filename)
    with pv_filename.open('w') as pf:
        pv_df.to_csv(pf)

    return str(pv_filename)


# endregion RUN_ASSOCIATION


# region AGGREGATE_RESULTS


def summarise_association_results(
    *pv_dfs: list[str],  
    pv_all_filename: str,
):
    """Summarise results

    Input:
    p-values from all association tests

    Ouput:
    one csv table per cell type,
    combining results across all genes in a single file
    """
    pv_all_df = pd.concat([pd.read_csv(AnyPath(pv_df), index_col=0) for pv_df in pv_dfs])

    # run qvalues for all tests
    pv_all_df['Q_CRM_VC'] = qvalue(pv_all_df['P_CRM_VC'])
    pv_all_df['Q_CRM_burden_max'] = qvalue(pv_all_df['P_CRM_burden_max'])
    pv_all_df['Q_CRM_burden_sum'] = qvalue(pv_all_df['P_CRM_burden_sum'])
    pv_all_df['Q_CRM_burden_comphet'] = qvalue(pv_all_df['P_CRM_burden_comphet'])
    pv_all_df['Q_CRM_omnibus_max'] = qvalue(pv_all_df['P_CRM_omnibus_max'])
    pv_all_df['Q_CRM_omnibus_sum'] = qvalue(pv_all_df['P_CRM_omnibus_sum'])
    pv_all_df['Q_CRM_omnibus_comphet'] = qvalue(pv_all_df['P_CRM_omnibus_comphet'])

    print(pv_all_df.head)

    pv_all_filename = AnyPath(pv_all_filename)
    print(pv_all_filename)
    with pv_all_filename.open('w') as pf:
        pv_all_df.to_csv(pf)


# endregion AGGREGATE_RESULTS

# region MISCELLANEOUS


def make_gene_loc_dict(file) -> dict[str, dict]:
    """
    Turn gene information into a dictionary
    to avoid opening this file for every gene
    """
    from csv import DictReader

    gene_dict = {}

    with open(to_path(file)) as handle:
        reader = DictReader(handle, delimiter='\t')

        for row in reader:
            gene_dict[row['gene_name']] = row

    return gene_dict


# can probably be merged with below
def extract_genes(gene_list, expression_tsv_path) -> list[str]:
    """
    Takes a list of all genes and subsets to only those
    present in the expression file of interest
    """
    expression_df = pd.read_csv(AnyPath(expression_tsv_path), sep='\t')
    expression_df = filter_lowly_expressed_genes(expression_df)
    gene_ids = set(list(expression_df.columns.values)[1:])
    genes = set(gene_list).intersection(gene_ids)

    print("Extracting genes")
    print(gene_list)
    print(gene_ids)
    print(list(sorted(genes)))

    return list(sorted(genes))


# copied from https://github.com/populationgenomics/tob-wgs/blob/main/scripts/eqtl_hail_batch/launch_eqtl_spearman.py
# generalised to specify min pct samples as input
def filter_lowly_expressed_genes(expression_df, min_pct=10):
    """Remove genes with low expression in all samples
    Input:
    expression_df: a data frame with samples as rows and genes as columns,
    containing normalised expression values (i.e., the average number of molecules
    for each gene detected in each person).
    Returns:
    A filtered version of the input data frame, after removing columns (genes)
    with 0 values in more than 10% of the rows (samples).
    """
    # Remove genes with 0 expression in all samples
    expression_df = expression_df.loc[:, (expression_df != 0).any(axis=0)]
    genes_not_equal_zero = expression_df.iloc[:, 1:].values != 0
    n_expr_over_zero = pd.DataFrame(genes_not_equal_zero.sum(axis=0))
    percent_expr_over_zero = (n_expr_over_zero / len(expression_df.index)) * 100
    percent_expr_over_zero.index = expression_df.columns[1:]

    # Filter genes with less than 10 percent individuals with non-zero expression
    atleastNpercent = percent_expr_over_zero[(percent_expr_over_zero > min_pct)[0]]
    sample_ids = expression_df['sampleid']
    expression_df = expression_df[atleastNpercent.index]
    expression_df.insert(loc=0, column='sampleid', value=sample_ids)

    return expression_df


def remove_sc_outliers(df, outliers=None):
    """
    Remove outlier samples, as identified by single-cell analysis
    """
    if outliers is None:
        outliers = ['966_967', '88_88']
    else:
        outliers = outliers.extend(['966_967', '88_88'])
    df = df[-df['OneK1K_ID'].isin(outliers)]

    return df


def qvalue(pv, m=None, verbose=False, lowmem=False, pi0=None):
    """
    Estimates q-values from p-values
    Args
    =====
    m: number of tests. If not specified m = pv.size
    verbose: print verbose messages? (default False)
    lowmem: use memory-efficient in-place algorithm
    pi0: if None, it's estimated as suggested in Storey and Tibshirani, 2003.
         For most GWAS this is not necessary, since pi0 is extremely likely to be
         1

    Taken from https://github.com/nfusi/qvalue
    but updated for python 3 (e.g., xrange -> range)

    Copyright (c) 2012, Nicolo Fusi, University of Sheffield
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the organization nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    """
    import scipy as sp
    from scipy import interpolate

    assert pv.min() >= 0 and pv.max() <= 1, 'p-values should be between 0 and 1'

    original_shape = pv.shape
    pv = pv.ravel()  # flattens the array in place, more efficient than flatten()

    if m is None:
        m = float(len(pv))
    else:
        # the user has supplied an m
        m *= 1.0

    # if the number of hypotheses is small, just set pi0 to 1
    if len(pv) < 100 and pi0 is None:
        pi0 = 1.0
    elif pi0 is not None:
        pi0 = pi0
    else:
        # evaluate pi0 for different lambdas
        pi0 = []
        lam = sp.arange(0, 0.90, 0.01)
        counts = sp.array([(pv > i).sum() for i in sp.arange(0, 0.9, 0.01)])
        for l in range(len(lam)):
            pi0.append(counts[l] / (m * (1 - lam[l])))

        pi0 = sp.array(pi0)

        # fit natural cubic spline
        tck = interpolate.splrep(lam, pi0, k=3)
        pi0 = interpolate.splev(lam[-1], tck)
        if verbose:
            logging.info(f'qvalues pi0={pi0}, estimated proportion of null features')

        if pi0 > 1:
            if verbose:
                logging.info(
                    f'got pi0 > 1 ({pi0}) while estimating qvalues, setting it to 1'
                )
            pi0 = 1.0

    assert pi0 >= 0 and pi0 <= 1, f'pi0 is not between 0 and 1: {pi0}'

    if lowmem:
        # low memory version, only uses 1 pv and 1 qv matrices
        qv = sp.zeros((len(pv),))
        last_pv = pv.argmax()
        qv[last_pv] = (pi0 * pv[last_pv] * m) / float(m)
        pv[last_pv] = -sp.inf
        prev_qv = last_pv
        for i in range(int(len(pv)) - 2, -1, -1):
            cur_max = pv.argmax()
            qv_i = pi0 * m * pv[cur_max] / float(i + 1)
            pv[cur_max] = -sp.inf
            qv_i1 = prev_qv
            qv[cur_max] = min(qv_i, qv_i1)
            prev_qv = qv[cur_max]

    else:
        p_ordered = sp.argsort(pv)
        pv = pv[p_ordered]
        qv = pi0 * m / len(pv) * pv
        qv[-1] = min(qv[-1], 1.0)

        for i in range(len(pv) - 2, -1, -1):
            qv[i] = min(pi0 * m * pv[i] / (i + 1.0), qv[i + 1])

        # reorder qvalues
        qv_temp = qv.copy()
        qv = sp.zeros_like(qv)
        qv[p_ordered] = qv_temp

    # reshape qvalues
    qv = qv.reshape(original_shape)

    return qv


# endregion MISCELLANEOUS

config = get_config()


@click.command()
@click.option('--celltypes')
@click.option('--expression-files-prefix', default='scrna-seq/grch38_association_files')
@click.option(
    '--sample-mapping-file-tsv',
    default='scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv',
)
@click.option('--mt-path', default=DEFAULT_JOINT_CALL_MT)
@click.option('--anno-ht-path', default=DEFAULT_ANNOTATION_HT)
@click.option(
    '--chromosomes',
    help='List of chromosome numbers to run rare variant association analysis on. '
    'Space separated, as one argument (Default: all)',
)
@click.option('--genes', default=None)
@click.option(
    '--max-gene-concurrency',
    type=int,
    default=50,
    help=(
        'To avoid resource starvation, set this concurrency to limit horizontal scale. '
        'Higher numbers have a better walltime, but risk jobs that are stuck (which are expensive)'
    ),
)
def crm_pipeline(
    celltypes: str,
    expression_files_prefix: str,
    sample_mapping_file_tsv: str,
    mt_path: str,
    anno_ht_path: str,
    chromosomes: str = 'all',
    genes: str | None = None,
    window_size: int = 50000,
    max_gene_concurrency=50,  # redundant default?
):

    sb = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch('CellRegMap pipeline', backend=sb)

    # extract samples for which we have single-cell (sc) data
    sample_mapping_file = pd.read_csv(dataset_path(sample_mapping_file_tsv), sep='\t')
    sample_mapping_file = remove_sc_outliers(sample_mapping_file)
    sc_samples = sample_mapping_file['InternalID'].unique()

    # filter to QC-passing, rare, biallelic variants
    output_mt_path = output_path('densified_rv_only.mt')
    if not to_path(output_mt_path).exists():

        filter_job = batch.new_python_job(name='MT filter job')
        copy_common_env(filter_job)
        filter_job.image(CELLREGMAP_IMAGE)
        filter_job.call(
            filter_variants,
            mt_path=mt_path,
            samples=list(sc_samples),
            output_mt_path=output_mt_path,
        )

    else:
        logging.info('File already exists no need to filter')
        filter_job = None

    # grab all relevant genes across all chromosomes
    # simpler if gene details are condensed to one file
    gene_dict: dict[str, dict] = {}
    if chromosomes == 'all':
        # autosomes only for now
        chromosomes_list = list(np.arange(22) + 1)
    else:
        chromosomes_list = chromosomes.split(' ')
    for chromosome in chromosomes_list:
        print(chromosome)
        geneloc_tsv_path = dataset_path(
            os.path.join(
                expression_files_prefix,
                'gene_location_files',
                f'GRCh38_geneloc_chr{chromosome}.tsv',
            )
        )

        # concatenating across chromosomes to have a single dict
        gene_dict.update(make_gene_loc_dict(geneloc_tsv_path))

    print(gene_dict.keys)

    # isolate to the genes we are interested in
    if genes is not None:
        genes_of_interest = genes.split(' ')
    else:
        genes_of_interest = list(
            gene_dict.keys()
        )

    print(genes_of_interest)

    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.job.Job] = []

    # move out of main?
    def manage_concurrency_for_job(job: hb.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_gene_concurrency:
            job.depends_on(_dependent_jobs[-max_gene_concurrency])
        _dependent_jobs.append(job)

    # for each gene, extract relevant variants (in window + with some annotation)
    # submit a job for each gene (export genotypes to plink)
    genotype_jobs = []
    for gene in genes_of_interest:
        print(gene)
        # final path for this gene - generate first (check syntax)
        plink_file = output_path(f'plink_files/{gene}')
        print(plink_file)
        gene_dict[gene]['plink'] = plink_file

        # if the plink output exists, do not re-generate it
        if to_path(f'{plink_file}.bim').exists():
            continue

        plink_job = batch.new_python_job(f'Create plink files for: {gene}')
        # gene_dict[gene]['plink_job'] = plink_job
        manage_concurrency_for_job(plink_job)
        copy_common_env(plink_job)
        if filter_job:
            plink_job.depends_on(filter_job)

        plink_job.image(CELLREGMAP_IMAGE)
        plink_job.call(
            get_promoter_variants,
            mt_path=output_mt_path,
            ht_path=anno_ht_path,
            gene_details=gene_dict[gene],
            window_size=window_size,
            plink_file=plink_file,
        )
        genotype_jobs.append(plink_job)

    # the next phase will be done for each cell type
    celltype_list = celltypes.split(' ')
    print(celltype_list)
    for celltype in celltype_list:
        expression_tsv_path = dataset_path(
            os.path.join(
                expression_files_prefix,
                'expression_files',
                f'{celltype}_expression.tsv',
            )
        )

        genes_list = extract_genes(genes_of_interest, expression_tsv_path)
        print(f"Genes to run: {genes_list}")
        if len(genes_list) == 0:
            print("No genes to run!")
            continue
        gene_run_jobs = []
        pv_files = []
        for gene in genes_list:

            pv_file = dataset_path(f'{celltype}/{gene}_pvals.csv')
            print(pv_file)
            if to_path(pv_file).exists():
                print("We already ran associations for this gene!")
                continue

            print(f"Preparing inputs for: {gene}")
            # TODO: add checks to not re-run genes if files already exist
            print(gene_dict[gene]['plink'])
            if gene_dict[gene]['plink'] is None:
                print("No plink files for this gene, exit!")
                continue
            plink_output_prefix = gene_dict[gene]['plink']
            # prepare input files
            prepare_input_job = batch.new_python_job(f'Prepare inputs for: {gene}')
            manage_concurrency_for_job(prepare_input_job)
            copy_common_env(prepare_input_job)
            prepare_input_job.depends_on(*genotype_jobs)
            # plink_dep = gene_dict[gene].get('plink_job')
            # if plink_dep:
            #     prepare_input_job.depends_on(plink_dep)
            prepare_input_job.image(CELLREGMAP_IMAGE)
            # the python_job.call only returns one object
            # the object is a file containing y_df, geno_df, kinship_df
            # all pickled into a file
            input_results = prepare_input_job.call(
                prepare_input_files,
                gene_name=gene,
                cell_type=celltype,
                genotype_file_bed=plink_output_prefix + '.bed',
                genotype_file_bim=plink_output_prefix + '.bim',
                genotype_file_fam=plink_output_prefix + '.fam',
                phenotype_file=expression_tsv_path,
                kinship_file=None,
                sample_mapping_file=sample_mapping_file_tsv,
            )
            print(f"Running association for: {gene}")
            # run association
            run_job = batch.new_python_job(f'Run association for: {gene}')
            manage_concurrency_for_job(run_job)
            copy_common_env(run_job)
            run_job.depends_on(prepare_input_job)
            run_job.image(CELLREGMAP_IMAGE)
            gene_run_jobs.append(run_job)
            run_job.call(
                run_gene_association,
                gene_name=gene,
                prepared_inputs=input_results,
                pv_path=pv_file
                # genotype_mat_path=geno_path,
                # phenotype_vec_path=pheno_path,
            )
            # save pv filename as gene attribute
        #     gene_dict[gene]['pv_file'] = pv_file
            pv_files.append(pv_file)

        print(pv_files)
        # combine all p-values across all chromosomes, genes (per cell type)
        summarise_job = batch.new_python_job(f'Summarise all results for {celltype}')
        copy_common_env(summarise_job)
        summarise_job.depends_on(*gene_run_jobs)
        summarise_job.image(CELLREGMAP_IMAGE)
        pv_all_filename_csv = str(output_path(f'{celltype}_all_pvalues.csv'))
        # print(pv_all_filename_csv)
        # print([gene_dict[gene]['pv_file'] for gene in genes_list])
        summarise_job.call(
            summarise_association_results,
            # *[gene_dict[gene]['pv_file'] for gene in genes_list],
            *pv_files,
            pv_all_filename=str(pv_all_filename_csv),
        )

    # set jobs running
    batch.run(wait=False)


if __name__ == '__main__':
    crm_pipeline()
