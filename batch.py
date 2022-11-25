#!/usr/bin/env python3
# pylint: disable=wrong-import-position,import-error

"""
Hail Batch workflow for the rare-variant association analysis, including:

- get relevant variants around a gene and export genotypes as plink files
- generate input files for association tests 
- run association tests
"""

# import python modules # is there a specific order, rationale for grouping imports??
import os
import re

import click
import logging

from cloudpathlib import AnyPath
from cpg_utils import to_path
from cpg_utils.hail_batch import (
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
from numpy import eye, ones
from pandas_plink import read_plink1_bin
from scipy.stats import shapiro

# qvalue imports
import scipy as sp
from scipy import interpolate

import hail as hl
import hailtop.batch as hb
from hail.methods import export_plink

from cellregmap import (  # figure out how to import this from github
    run_gene_set_association,
    run_burden_association,
    omnibus_set_association,
)

# use logging to print statements, display at info level
logging.basicConfig(
    format="%(levelname)s:%(message)s", level=logging.INFO
)  # consider removing some print statements?

DEFAULT_JOINT_CALL_MT = dataset_path("mt/v7.mt")
DEFAULT_ANNOTATION_HT = dataset_path(
    "tob_wgs_vep/104/vep104.3_GRCh38.ht"
)  # atm VEP only

CELLREGMAP_IMAGE = (
    "australia-southeast1-docker.pkg.dev/cpg-common/images/cellregmap:0.0.3"
)

# region SUBSET_VARIANTS


def filter_variants(
    mt_path: str,  # "mt/v7.mt"
    samples: list[str],
    output_mt_path: str,  # "tob_wgs_rv/densified_rv_only.mt"
):
    """Subset hail matrix table

    Input:
    joint call hail matrix table
    set of samples for which we have scRNA-seq data

    Output:
    subsetted hail matrix table, containing only variants that:
    1) are not ref-only, 2) biallelic, 3) meet QC filters, 4) are rare (MAF<5%)
    """
    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    # subset to relevant samples (samples we have scRNA-seq data for)
    mt = mt.filter_cols(hl.set(samples).contains(mt.s))

    # densify (can this be done at the end?)
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
    mt = mt.checkpoint(output_mt_path, overwrite=True)
    logging.info(
        f"Number of rare (freq<5%) and QC'd biallelic variants: {mt.count()[0]}"
    )


# endregion SUBSET_VARIANTS


# region GET_GENE_SPECIFIC_VARIANTS


def get_promoter_variants(
    mt_path: str,  # checkpoint from function above
    ht_path: str,
    gene_dict: dict[str, str],  # ouput of make_gene_loc_dict
    window_size: int,
    plink_output_prefix: str,  # 'tob_wgs_rv/pseudobulk_rv_association/{celltype}/'
):
    """Subset hail matrix table

    Input:
    mt_path: path to already subsetted hail matrix table
    gene_file: path to a tsv file with information on
    the chrom, start and end location of genes (rows)

    Output:
    For retained variants, that are: 1) in promoter regions and
    2) within 50kb up or down-stream of the gene body (or in the gene body itself)
    (on top of all filters done above)

    returns:
    - Path to plink prefix for genotype files (plink format: .bed, .bim, .fam)
    - Path to hail table with variant (rows) stats for downstream analyses
    """
    # read hail matrix table object (pre-filtered)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    gene_name = gene_dict["gene_name"]

    # get relevant chromosome
    chrom = gene_dict["chrom"]

    # subset to window
    # get gene body position (start and end) and build interval
    left_boundary = max(1, int(gene_dict["gene-start"]) - window_size)
    right_boundary = min(
        int(gene_dict["gene-end"]) + window_size,
        hl.get_reference("GRCh38").lengths[f"chr{chrom}"],
    )
    # get gene-specific genomic interval
    gene_interval = f"{chrom}:{left_boundary}-{right_boundary}"
    logging.info(f"Interval considered: {gene_interval}")  # 'chr22:23219960-23348287'

    # include variants up to {window size} up- and downstream
    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval(gene_interval, reference_genome="GRCh38")]
    )
    mt_path = output_path("in_window.mt", "tmp")
    mt = mt.checkpoint(
        mt_path, overwrite=True
    )  # add checkpoint to avoid repeat evaluation
    logging.info(f"Number of variants within interval: {mt.count()[0]}")

    # annotate using VEP
    vep_ht = hl.read_table(ht_path)
    mt = mt.annotate_rows(vep=vep_ht[mt.row_key].vep)

    # filter variants found to be in promoter regions
    mt = mt.filter_rows(
        mt.vep.regulatory_feature_consequences["biotype"].contains("promoter")
    )
    promoter_path = output_path(f"{gene_name}promoter_variants.mt", "tmp")
    mt = mt.checkpoint(promoter_path, overwrite=True)  # checkpoint
    logging.info(
        f"Number of rare (freq<5%) QC-passing, biallelic SNPs in promoter regions: {mt.count()[0]}"
    )

    # export this as a Hail table for downstream analysis
    ht_path = output_path(f"{gene_name}_rare_promoter_summary.ht")
    ht = mt.rows()
    ht.write(ht_path)

    # export MT object to PLINK (promoter variants)
    plink_path = output_path(f"{plink_output_prefix}/{gene_name}")
    export_plink(mt, plink_path, ind_id=mt.s)

    return plink_path


# endregion GET_GENE_SPECIFIC_VARIANTS

# region PREPARE_INPUT_FILES


def prepare_input_files(  # consider splitting into multiple functions
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
    expression_filename = AnyPath(output_path(f"{gene_name}_{cell_type}.csv"))
    genotype_filename = AnyPath(output_path(f"{gene_name}_rare_regulatory.csv"))
    kinship_filename = AnyPath(output_path("kinship_common_samples.csv"))

    # read in phenotype file (tsv)
    phenotype = pd.read_csv(phenotype_file, sep="\t", index_col=0)

    phenotype = xr.DataArray(
        phenotype.values,
        dims=["sample", "gene"],
        coords={"sample": phenotype.index.values, "gene": phenotype.columns.values},
    )

    # read in genotype file (plink format)
    to_path(genotype_file_bed).copy("temp.bed")  # bed
    to_path(genotype_file_bim).copy("temp.bim")  # bim
    to_path(genotype_file_fam).copy("temp.fam")  # fam
    geno = read_plink1_bin("temp.bed")

    # read in GRM (genotype relationship matrix; kinship matrix)
    kinship = pd.read_csv(kinship_file, index_col=0)
    kinship.index = kinship.index.astype("str")
    assert all(kinship.columns == kinship.index)  # symmetric matrix, donors x donors
    kinship = xr.DataArray(
        kinship.values,
        dims=["sample_0", "sample_1"],
        coords={"sample_0": kinship.columns, "sample_1": kinship.index},
    )
    kinship = kinship.sortby("sample_0").sortby("sample_1")

    # this file will map different IDs (and OneK1K ID to CPG ID)
    sample_mapping = pd.read_csv(sample_mapping_file, sep="\t")

    # ensure samples are the same and in the same order across input files
    # samples with expression data
    donors_exprs = set(phenotype.sample.values).intersection(
        set(sample_mapping["OneK1K_ID"].unique())
    )

    logging.info(f"Number of unique donors with expression data: {len(donors_exprs)}")

    # samples with genotype data
    donors_geno = set(geno.sample.values).intersection(
        set(sample_mapping["InternalID"].unique())
    )
    logging.info(f"Number of unique donors with genotype data: {len(donors_geno)}")

    # samples with both (can this be done in one step?)
    sample_mapping1 = sample_mapping.loc[sample_mapping["OneK1K_ID"].isin(donors_exprs)]
    sample_mapping_both = sample_mapping1.loc[
        sample_mapping1["InternalID"].isin(donors_geno)
    ]
    donors_e = sample_mapping_both["OneK1K_ID"].unique()
    donors_g = sample_mapping_both["InternalID"].unique()
    assert len(donors_e) == len(donors_g)

    # samples in kinship
    donors_e_short = [re.sub(".*_", "", donor) for donor in donors_e]
    donors_k = sorted(set(list(kinship.sample_0.values)).intersection(donors_e_short))

    logging.info(f"Number of unique common donors: {len(donors_g)}")

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
    with expression_filename.open("w") as ef:
        y_df.to_csv(ef, index=False)

    with genotype_filename.open("w") as gf:
        geno_df.to_csv(gf, index=False)

    with kinship_filename.open("w") as kf:
        kinship_df.to_csv(kf, index=False)

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
        y=pheno, G=genotypes, W=covs, E=contexts, mask="mask.max"
    )[0]
    pv2 = run_burden_association(
        y=pheno, G=genotypes, W=covs, E=contexts, mask="mask.sum"
    )[0]
    pv3 = run_burden_association(
        y=pheno, G=genotypes, W=covs, E=contexts, mask="mask.comphet"
    )[0]
    pv4 = omnibus_set_association(np.array([pv0, pv1]))
    pv5 = omnibus_set_association(np.array([pv0, pv2]))
    pv6 = omnibus_set_association(np.array([pv0, pv3]))

    return np.array([pv_norm, pv0, pv1, pv2, pv3, pv4, pv5, pv6])


# endregion GET_CRM_PVALUES

# region RUN_ASSOCIATION


def run_gene_association(
    gene_name: str,  # 'VPREB3'
    genotype_mat_path: str,  # 'VPREB3_50K_window/SNVs.csv'
    phenotype_vec_path: str,  # 'Bnaive/VPREB3_pseudocounts.csv'
):
    """Run gene-set association test

    Input:
    input files

    Output:
    table with p-values
    """
    # get genotypes
    # these are variants in and around gene VPREB3 on chrom 22
    g_file = AnyPath(output_path(genotype_mat_path))
    g_df = pd.read_csv(g_file)
    # because the current matrix is counting the copies of the reference allele
    # while we are interested in the alternative allele, flip the genotypes
    genotypes = 2 - g_df

    # get phenptypes
    p_file = AnyPath(output_path(phenotype_vec_path))
    p_df = pd.read_csv(p_file)
    pheno = p_df.values

    # covs
    covs = ones((genotypes.shape[0], 1))  # intercept of ones as covariates

    # contexts
    contexts = eye(genotypes.shape[0])

    # TODO: kinship

    cols = [
        "P_shapiro",
        "P_CRM_VC",
        "P_CRM_burden_max",
        "P_CRM_burden_sum",
        "P_CRM_burden_comphet",
        "P_CRM_omnibus_max",
        "P_CRM_omnibus_sum",
        "P_CRM_omnibus_comphet",
    ]

    # create p-values data frame
    pv_df = pd.DataFrame(
        data=get_crm_pvs(pheno, covs, genotypes, contexts),
        columns=cols,
        index=gene_name,
    )

    pv_filename = AnyPath(
        output_path(
            "simulations/CRM/1000samples_10causal_singletons/10tested_samebeta.csv"
        )
    )
    with pv_filename.open("w") as pf:
        pv_df.to_csv(pf, index=False)

    return pv_filename


# endregion RUN_ASSOCIATION


# region AGGREGATE_RESULTS


def summarise_association_results(
    pv_dfs: list[str],
):
    """Summarise results

    Input:
    p-values from all association tests

    Ouput:
    one matrix (what format??) per cell type,
    combining results across all genes in a single file
    """
    pv_all_df = pd.concat([pv_dfs])  # test syntax

    # run qvalues for all tests
    pv_all_df["Q_CRM_VC"] = qvalue(pv_all_df["P_CRM_VC"])
    pv_all_df["Q_CRM_burden_max"] = qvalue(pv_all_df["P_CRM_burden_max"])
    pv_all_df["Q_CRM_burden_sum"] = qvalue(pv_all_df["P_CRM_burden_sum"])
    pv_all_df["Q_CRM_burden_comphet"] = qvalue(pv_all_df["P_CRM_burden_comphet"])
    pv_all_df["Q_CRM_omnibus_max"] = qvalue(pv_all_df["P_CRM_omnibus_max"])
    pv_all_df["Q_CRM_omnibus_sum"] = qvalue(pv_all_df["P_CRM_omnibus_sum"])
    pv_all_df["Q_CRM_omnibus_comphet"] = qvalue(pv_all_df["P_CRM_omnibus_comphet"])

    return pv_all_df


# endregion AGGREGATE_RESULTS

# region MISCELLANEOUS


def make_gene_loc_dict(file) -> dict[str, dict]:
    from csv import DictReader

    gene_dict = {}

    with open(file) as handle:
        reader = DictReader(handle, delimiter="\t")

        for row in reader:
            gene_dict[row["gene_name"]] = row
            # gene_dict[row["gene_name"]] = row["gene_name"]
            # gene_dict[row["chrom"]] = row["chr"]
            # gene_dict[row["gene_start"]] = row["start"]
            # gene_dict[row["gene_end"]] = row["end"]

    return gene_dict


# copied from https://github.com/populationgenomics/tob-wgs/blob/main/scripts/eqtl_hail_batch/launch_eqtl_spearman.py
# check whether it needs modifying
def get_genes_for_chromosome(*, expression_tsv_path, geneloc_tsv_path) -> list[str]:
    """get index of total number of genes
    Input:
    expression_df: a data frame with samples as rows and genes as columns,
    containing normalised expression values (i.e., the average number of molecules
    for each gene detected in each person).
    geneloc_df: a data frame with the gene location (chromosome, start, and
    end locations as columns) specified for each gene (rows).
    Returns:
    The number of genes (as an int) after filtering for lowly expressed genes.
    This integer number gets fed into the number of scatters to run.
    """
    expression_df = pd.read_csv(AnyPath(expression_tsv_path), sep="\t")
    geneloc_df = pd.read_csv(AnyPath(geneloc_tsv_path), sep="\t")

    # expression_df = filter_lowly_expressed_genes(expression_df) # I might add this in but not for now
    gene_ids = set(list(expression_df.columns.values)[1:])

    genes = set(geneloc_df.gene_name).intersection(gene_ids)
    return list(sorted(genes))


# copied from https://github.com/populationgenomics/tob-wgs/blob/main/scripts/eqtl_hail_batch/launch_eqtl_spearman.py
# check whether it needs modifying
def remove_sc_outliers(df):
    """
    Remove outlier samples, as identified by sc analysis
    """
    outliers = ["966_967", "88_88"]
    df = df[-df.sampleid.isin(outliers)]

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

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
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

    assert pv.min() >= 0 and pv.max() <= 1, "p-values should be between 0 and 1"

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
            # print("qvalues pi0=%.3f, estimated proportion of null features " % pi0)
            logging.info(
                "qvalues pi0=%.3f, estimated proportion of null features " % pi0
            )

        if pi0 > 1:
            if verbose:
                # print(
                #     "got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0
                # )
                logging.info(
                    "got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0
                )
            pi0 = 1.0

    assert pi0 >= 0 and pi0 <= 1, "pi0 is not between 0 and 1: %f" % pi0

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
@click.option("--sc-samples")
@click.option("--chromosomes")
@click.option("--genes")
@click.option("--celltypes")
@click.option("--expression-file-prefix")
@click.option("--sample-mapping-file")
@click.option("--mt-path")
@click.option("--anno-ht-path")
# @click.option("--fdr-threshold")
def main(
    sc_samples: list["str"],
    chromosomes: list[str],
    genes: list[str],
    celltypes: list[str],
    expression_files_prefix: str,  # 'scrna-seq/grch38_association_files'
    sample_mapping_file: str,
    mt_path: str = DEFAULT_JOINT_CALL_MT,  # 'mt/v7.mt'
    anno_ht_path: str = DEFAULT_ANNOTATION_HT,  # 'tob_wgs_vep/104/vep104.3_GRCh38.ht'
    # fdr_threshold: float = 0.05,
):

    sb = hb.ServiceBackend(
        billing_project=config["hail"]["billing_project"],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch("CellRegMap pipeline", backend=sb)

    mt = filter_variants(mt_path=mt_path, samples=sc_samples)

    # step1: for each gene, extract relevant variants
    genotype_jobs = []

    # loop over chromosomes
    for chromosome in chromosomes:
        geneloc_tsv_path = os.path.join(
            expression_files_prefix,
            "gene_location_files",
            f"GRCh38_geneloc_chr{chromosome}.tsv",
        )
        geneloc_dict = make_gene_loc_dict(geneloc_tsv_path)

        # find genes in that chromosome
        # loop over genes
        _genes = genes or get_genes_for_chromosome(
            expression_tsv_path=expression_tsv_path,
            geneloc_tsv_path=geneloc_tsv_path,
        )

        # submit a job for each gene (export genotypes to plink)
        plink_paths = {}
        for gene in _genes:
            plink_job = batch.new_python_job(f"Create plink files for: {gene}")
            genotype_jobs.append(plink_job)
            plink_path = plink_job.call(
                get_promoter_variants,
                mt_input_path=mt,  # ouput of filter_variants
                anno_ht_path=anno_ht_path,  # 'tob_wgs_vep/104/vep104.3_GRCh38.ht'
                gene_dict=geneloc_dict[gene],
                window_size=50000,
                # plink_output_prefix=plink_output_prefix[gene],
            )
            plink_paths[gene] = plink_path  # syntax??

    for celltype in celltypes:
        expression_tsv_path = os.path.join(
            expression_files_prefix, "expression_files", f"{celltype}_expression.tsv"
        )

        genes = extract_genes(expression_files_prefix)
        pv_file_paths = []
        # maybe this makes no sense (to loop over genes twice)
        # but I also didn't want to re-select variants for the same gene repeatedly
        # for every new cell type?
        gene_prepare_jobs = []
        gene_run_jobs = []
        for gene in genes:
            plink_output_prefix = plink_paths[gene]
            # prepare input files
            prepare_input_job = batch.new_python_job(f"Prepare inputs for: {gene}")
            gene_prepare_jobs.append(prepare_input_job)
            pheno_path, geno_path, _ = prepare_input_job.call(
                prepare_input_files,
                gene_name=gene,
                cell_type=celltype,
                genotype_file_bed=plink_output_prefix + ".bed",
                genotype_file_bim=plink_output_prefix + ".bim",
                genotype_file_fam=plink_output_prefix + ".fam",
                phenotype_file=expression_tsv_path,
                kinship_file=None,
                sample_mapping_file=sample_mapping_file,
            )
            # run association
            run_job = batch.new_python_job(f"Run association for: {gene}")
            run_job.image(CELLREGMAP_IMAGE)
            gene_run_jobs.append(run_job)
            pv_file = run_job.call(
                run_gene_association,
                gene_name=gene,
                genotype_mat_path=geno_path,
                phenotype_vec_path=pheno_path,
            )
            pv_file_paths.append(pv_file)
        # combine all p-values across all chromosomes, genes (per cell type)
        summarise_job = batch.new_python_job(f"Summarise all results for {celltype}")
        pv_all = summarise_job.call(
            summarise_association_results, pv_dfs=pv_file_paths
        )  # no idea how do to this (get previous job's dataframes and add them in a list)

        pv_filename = AnyPath(output_path(f"{celltype}_all_pvalues.csv"))
        with pv_filename.open("w") as pf:
            pv_all.to_csv(pf, index=False)

    # dependencies between jobs
    prepare_input_job.depends_on(*genotype_jobs)

    # this only depends the corresponding job (for that same gene)??
    # do this inside loop??
    run_job.depends_on(prepare_input_job)

    summarise_job.depends_on(*gene_run_jobs)

    # determine for each chromosome

    # below is pseudocode from @illusional I don't understand but don't want to delete

    # plink_job_reference_and_outputs_for_gene: dict[str, tuple[hb.Job, str]] = {}
    # if you know the output, you need to wait for the actual job that processes it
    # so store the reference to the job that's producing the file

    # for chr_loc_file in gene_loc_files:
    #     with AnyPath(chr_loc_file).open() as f:
    #         for line in f:
    #             # i have no idea, I'm just randomly guessing
    #             gene = line.split("\t")[0]
    #             job = batch.new_python_job(f"Preprocess: {gene}")

    #             plink_output_prefix = output_path(f"plink_files/{gene}_rare_promoter")
    #             plink_job_reference_and_outputs_for_gene[gene] = (
    #                 job,
    #                 plink_output_prefix,
    #             )
    #             result = job.call(
    #                 get_promoter_variants,
    #                 gene_file=chr_loc_file,
    #                 gene_name=gene,
    #                 window_size=50_000,
    #                 plink_output_prefix=plink_output_prefix,
    #             )

    # _gene = "BRCA1"
    # dependent_job, output_prefix = plink_job_reference_and_outputs_for_gene[_gene]
    # downstream_job = batch.new_job("Downstream job")
    # downstream_job.depends_on(dependent_job)

    # # run that function


if __name__ == "__main__":
    main()  # rename
