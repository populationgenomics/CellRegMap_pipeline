#!/usr/bin/env python3
# pylint: disable=wrong-import-position,import-error

"""
Hail Batch workflow for the rare-variant association analysis, including:

- get relevant variants around a gene and export genotypes as plink files
- generate input files for association tests 
- run association tests
"""

# import python modules # is there a specific order, rationale for grouping imports??
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
import xarray as xr  # what am I doing here?? install this in a new image? build an image within this??

from limix.qc import quantile_gaussianize
from numpy import eye, ones
from pandas_plink import read_plink1_bin
from scipy.stats import shapiro

import hail as hl
import hailtop.batch as hb
from hail.methods import export_plink

from cellregmap import (  # figure out how to import this (in image?)
    run_gene_set_association,
    run_burden_association,
    omnibus_set_association,
)

# use logging to print statements, display at info level
logging.basicConfig(
    format="%(levelname)s:%(message)s", level=logging.INFO
)  # consider removing some print statements?

DEFAULT_JOINT_CALL_MT = dataset_path("mt/v7.mt")
DEFAULT_ANNOTATION_HT = dataset_path("tob_wgs_vep/104/vep104.3_GRCh38.ht")  # atm VEP

# region GET_RELEVANT_VARIANTS


def get_promoter_variants(
    mt_path: str,
    ht_path: str,
    gene_file: str,
    gene_name: str,
    window_size: int,
    samples: list[str],
    # plink_output_prefix: str, # figure out how to deal with this
):
    """Subset hail matrix table

    Input:
    gene_file: path to a tsv file with information on
    the chrom, start and end location of genes (rows)

    Output:
    For retained variants, that are: 1) not ref-only, 2) biallelic,
    3) meet QC filters, 4) rare (MAF<5%), 5) in promoter regions and
    6) within 50kb up or down-stream of the gene body (or in the gene body itself)
    returns:
    - Path to plink prefix for genotype files (plink format: .bed, .bim, .fam)
    - Path to hail table with variant (rows) stats for downstream analyses
    """
    # TODO: subset to relevant samples (at least before frequency filter)
    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    # get relevant chromosome
    gene_df = pd.read_csv(AnyPath(dataset_path(gene_file)), sep="\t", index_col=0)
    chrom = gene_df[gene_df["gene_name"] == gene_name]["chr"]
    # subset to chromosome
    mt = mt.filter_rows(mt.locus.contig == ("chr" + chrom))

    # densify
    mt = hl.experimental.densify(mt)  # never know when i should do this step

    # subset to window
    # get gene body position (start and end) and build interval
    interval_start = int(gene_df[gene_df["gene_name"] == gene_name]["start"]) - int(
        window_size
    )
    interval_end = int(gene_df[gene_df["gene_name"] == gene_name]["end"]) + int(
        window_size
    )
    # clip to chromosome boundaries
    left_boundary = max(1, interval_start)
    right_boundary = min(
        interval_end, hl.get_reference("GRCh38").lengths[f"chr{chrom}"]
    )
    # get gene-specific genomic interval
    gene_interval = f"chr{chrom}:{left_boundary}-{right_boundary}"
    logging.info(f"Interval considered: {gene_interval}")  # 'chr22:23219960-23348287'

    # include variants up to {window size} up- and downstream
    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval(gene_interval, reference_genome="GRCh38")]
    )
    logging.info(f"Number of variants within interval: {mt.count()[0]}")
    logging.info(f"Number of variants in window: {mt.count()[0]}")

    # filter out low quality variants and consider biallelic variants only (no multi-allelic, no ref-only)
    mt = mt.filter_rows(
        (hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)  # QC
        & (hl.len(mt.alleles) == 2)  # remove hom-ref
        & (mt.n_unsplit_alleles == 2)  # biallelic
        & (hl.is_snp(mt.alleles[0], mt.alleles[1]))  # SNVs
    )

    # annotate using VEP
    vep_ht = hl.read_table(ht_path)
    mt = mt.annotate_rows(vep=vep_ht[mt.row_key].vep)
    mt_path = output_path("densified_qced_snps_only_vep_annotated.mt", "tmp")
    mt = mt.checkpoint(
        mt_path, overwrite=True
    )  # add checkpoint to avoid repeat evaluation
    logging.info(f"Number of QC-passing, biallelic SNPs: {mt.count()[0]}")

    # filter variants found to be in promoter regions
    mt = mt.filter_rows(
        mt.vep.regulatory_feature_consequences["biotype"][0] == "promoter"
    )
    mt_path = output_path("promoter_variants.mt", "tmp")
    mt = mt.checkpoint(mt_path, overwrite=True)  # checkpoint
    logging.info(
        f"Number of rare variants (freq<5%) in promoter regions: {mt.count()[0]}"
    )

    # filter rare variants only (MAF < 5%)
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] < 0.05) & (mt.variant_qc.AF[1] > 0)
        | (mt.variant_qc.AF[1] > 0.95) & (mt.variant_qc.AF[1] < 1)
    )
    mt_path = output_path("rare_variants.mt", "tmp")
    mt = mt.checkpoint(mt_path, overwrite=True)  # checkpoint
    logging.info(f"Number of rare variants (freq<5%): {mt.count()[0]}")

    # export this as a Hail table for downstream analysis
    ht_filename = output_path(f"{gene_name}_rare_promoter_summary.ht")
    ht = mt.rows()
    ht.write(ht_filename)

    # export MT object to PLINK (promoter variants)
    mt_path = output_path(f"plink_files/{gene_name}_rare_promoter")
    export_plink(mt, mt_path, ind_id=mt.s)

    return [mt_path, ht_filename]


# endregion GET_RELEVANT_VARIANTS

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
    expression_filename = AnyPath(output_path(f'{gene_name}_{cell_type}.csv'))
    genotype_filename = AnyPath(output_path(f'{gene_name}_rare_regulatory.csv'))
    kinship_filename = AnyPath(output_path('kinship_common_samples.csv'))

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
    donors_e_short = [re.sub('.*_', '', donor) for donor in donors_e]
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
    with expression_filename.open('w') as ef:  
        y_df.to_csv(ef, index=False)

    with genotype_filename.open('w') as gf:  
        geno_df.to_csv(gf, index=False)

    with kinship_filename.open('w') as kf: 
        kinship_df.to_csv(kf, index=False)
    
    return [y_df, geno_df, kinship_df]

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
    return [pv_norm, pv0, pv1, pv2, pv3, pv4, pv5, pv6]  # do I need brackets?


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
    fdrThreshold: int = 1,
):
    """Summarise results

    Input:
    p-values from all association tests

    Ouput:
    one matrix (what format??) per cell type,
    combining results across all genes in a single file
    """
    for pv_df in pv_dfs:
        pv_all_df = pd.concat(pv_df)

    # run qvalues for all tests


# endregion AGGREGATE_RESULTS


config = get_config()


@click.command()
@click.option("--gene-list")
@click.option("--celltype-list")
@click.option("--mt-path")
@click.option("--anno_ht_path")
@click.option("--fdr-threshold")
def main(
    gene_list: list[str],
    celltype_list: list[str],
    mt_path: str = DEFAULT_JOINT_CALL_MT,  # 'mt/v7.mt'
    anno_ht_path: str = DEFAULT_ANNOTATION_HT,  # 'tob_wgs_vep/104/vep104.3_GRCh38.ht'
    fdr_threshold: float = 0.05,
):

    sb = hb.ServiceBackend(
        billing_project=config["hail"]["billing_project"],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch("CellRegMap pipeline", backend=sb)

    # determine for each chromosome

    plink_job_reference_and_outputs_for_gene: dict[str, tuple[hb.Job, str]] = {}
    # if you know the output, you need to wait for the actual job that processes it
    # so store the reference to the job that's producing the file

    for chr_loc_file in gene_loc_files:
        with AnyPath(chr_loc_file).open() as f:
            for line in f:
                # i have no idea, I'm just randomly guessing
                gene = line.split("\t")[0]
                job = batch.new_python_job(f"Preprocess: {gene}")

                plink_output_prefix = output_path(f"plink_files/{gene}_rare_promoter")
                plink_job_reference_and_outputs_for_gene[gene] = (
                    job,
                    plink_output_prefix,
                )
                result = job.call(
                    get_promoter_variants,
                    gene_file=chr_loc_file,
                    gene_name=gene,
                    window_size=50_000,
                    plink_output_prefix=plink_output_prefix,
                )

    _gene = "BRCA1"
    dependent_job, output_prefix = plink_job_reference_and_outputs_for_gene[_gene]
    downstream_job = batch.new_job("Downstream job")
    downstream_job.depends_on(dependent_job)

    # run that function
