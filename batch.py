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

# make gene loc dict import
from csv import DictReader

import hail as hl
import hailtop.batch as hb
from hail.methods import export_plink

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
    subset hail matrix table, containing only variants that:
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
    gene_details: dict[str, str],  # ouput of make_gene_loc_dict
    window_size: int,
    plink_output_prefix: str,  # 'tob_wgs_rv/pseudobulk_rv_association/plink_files/'
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

    gene_name = gene_details["gene_name"]

    # get relevant chromosome
    chrom = gene_details["chr"]

    # subset to window
    # get gene body position (start and end) and build interval
    left_boundary = max(1, int(gene_details["start"]) - window_size)
    right_boundary = min(
        int(gene_details["end"]) + window_size,
        hl.get_reference("GRCh38").lengths[f"chr{chrom}"],
    )
    # get gene-specific genomic interval
    gene_interval = f"{chrom}:{left_boundary}-{right_boundary}"
    logging.info(f"Interval considered: {gene_interval}")  # 'chr22:23219960-23348287'

    # include variants up to {window size} up- and downstream
    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval(gene_interval, reference_genome="GRCh38")]
    )
    mt_path = output_path(f"{gene_name}_in_window.mt", "tmp")
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

# region MISCELLANEOUS


def make_gene_loc_dict(file) -> dict[str, dict]:
    """
    Turn gene information into a dictionary
    to avoid opening this file for every gene
    """
    gene_dict = {}

    with open(file) as handle:
        reader = DictReader(handle, delimiter="\t")

        for row in reader:
            gene_dict[row["gene_name"]] = row

    return gene_dict


# copied from https://github.com/populationgenomics/tob-wgs/blob/main/scripts/eqtl_hail_batch/launch_eqtl_spearman.py
# check whether it needs modifying
def remove_sc_outliers(df):
    """
    Remove outlier samples, as identified by sc analysis
    """
    outliers = ["966_967", "88_88"]
    df = df[-df.sampleid.isin(outliers)]

    return df


# endregion MISCELLANEOUS



@click.command()
@click.option("--chromosomes")
@click.option("--genes")
@click.option("--celltypes")
@click.option("--expression-file-prefix")
@click.option("--sample-mapping-file")
@click.option("--output-path")
@click.option("--mt-path")
@click.option("--anno-ht-path")
def crm_pipeline(
    chromosomes: list[str],
    genes: list[str],
    celltypes: list[str],
    expression_files_prefix: str,  # 'scrna-seq/grch38_association_files'
    sample_mapping_file: str,
    output_path: str,  # 'tob_wgs_rv/pseudobulk_rv_association'
    mt_path: str = DEFAULT_JOINT_CALL_MT,  # 'mt/v7.mt'
    anno_ht_path: str = DEFAULT_ANNOTATION_HT,  # 'tob_wgs_vep/104/vep104.3_GRCh38.ht'
    window_size: int = 50000,
):

    sb = hb.ServiceBackend(
        billing_project=get_config()["hail"]["billing_project"],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch("CellRegMap pipeline", backend=sb)

    # extract samples for which we have single-cell (sc) data
    sample_mapping_file = remove_sc_outliers(sample_mapping_file)
    sc_samples = sample_mapping_file["InternalID"].unique()

    # filter to QC-passing, rare, biallelic variants
    filter_job = batch.new_python_job(name="MT filter job")
    output_mt_path = output_path("tob_wgs_rv/densified_rv_only.mt")
    filter_job.image(CELLREGMAP_IMAGE)
    filter_job.call(
        filter_variants,
        mt_path=mt_path,
        samples=sc_samples,
        output_mt_path=output_mt_path,
    )

    # grab all relevant genes across all chromosomes
    # simpler if gene details are condensed to one file
    gene_dict: dict[str, dict] = {}
    for chromosome in chromosomes:
        geneloc_tsv_path = os.path.join(
            expression_files_prefix,
            "gene_location_files",
            f"GRCh38_geneloc_chr{chromosome}.tsv",
        )

        # concatenating across chromosomes to have a single dict
        gene_dict.update(make_gene_loc_dict(geneloc_tsv_path))

    # isolate to the genes we are interested in
    genes_of_interest = genes or list[gene_dict.keys()]

    # for each gene, extract relevant variants (in window + with some annotation)
    # submit a job for each gene (export genotypes to plink)
    genotype_jobs = []
    for gene in genes_of_interest:

        # final path for this gene - generate first (check syntax)
        plink_file = output_path(f"{output_path}/plink_files/{gene}")
        gene_dict[gene]["plink"] = plink_file

        # if the plink output exists, do not re-generate it
        if to_path(f'{plink_file}.bim').exists():
            continue

        plink_job = batch.new_python_job(f"Create plink files for: {gene}")
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
        # concatenate jobs so they can be depended on
        genotype_jobs.append(plink_job)

    # set jobs running
    batch.run(wait=False)


if __name__ == "__main__":
    crm_pipeline()
