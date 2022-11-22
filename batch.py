#!/usr/bin/env python3
# pylint: disable=wrong-import-position,import-error 

"""
Hail Batch workflow for the rare-variant association analysis, including:

- get relevant variants around a gene and export genotypes as plink files
- generate input files for association tests 
- run association tests
"""

# import python modules # is there a specific order, rationale for grouping imports??
import logging
import click
import hail as hl
import hailtop.batch as hb
import pandas as pd
from hail.methods import export_plink
from cpg_utils.hail_batch import get_config, remote_tmpdir
from cloudpathlib import AnyPath
from cpg_utils.hail_batch import dataset_path, init_batch, output_path

# use logging to print statements, display at info level
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

DEFAULT_JOINT_CALL_MT = dataset_path('mt/v7.mt')
DEFAULT_ANNOTATION_HT = dataset_path('tob_wgs_vep/104/vep104.3_GRCh38.ht')  # atm VEP

# region GET_RELEVANT_VARIANTS

def get_promoter_variants(
    mt_path: str,
    ht_path: str,
    gene_file: str,
    gene_name: str,
    window_size: int,
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
    # read hail matrix table object (WGS data)
    init_batch()
    mt = hl.read_matrix_table(mt_path)

    # get relevant chromosome
    gene_df = pd.read_csv(AnyPath(dataset_path(gene_file)), sep='\t', index_col=0)
    chrom = gene_df[gene_df['gene_name'] == gene_name]['chr']
    # subset to chromosome
    mt = mt.filter_rows(mt.locus.contig == ('chr' + chrom))

    # densify
    mt = hl.experimental.densify(mt)  # never know when i should do this step

    # subset to window
    # get gene body position (start and end) and build interval
    interval_start = int(gene_df[gene_df['gene_name'] == gene_name]['start']) - int(
        window_size
    )
    interval_end = int(gene_df[gene_df['gene_name'] == gene_name]['end']) + int(
        window_size
    )
    # clip to chromosome boundaries
    left_boundary = max(1, interval_start)
    right_boundary = min(
        interval_end, hl.get_reference('GRCh38').lengths[f'chr{chrom}']
    )
    # get gene-specific genomic interval
    gene_interval = f'chr{chrom}:{left_boundary}-{right_boundary}'
    logging.info(f'Interval considered: {gene_interval}')  # 'chr22:23219960-23348287'

    # include variants up to {window size} up- and downstream
    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval(gene_interval, reference_genome='GRCh38')]
    )
    logging.info(f'Number of variants within interval: {mt.count()[0]}')
    logging.info(f'Number of variants in window: {mt.count()[0]}')

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
    mt_path = output_path('densified_qced_snps_only_vep_annotated.mt', 'tmp')
    mt = mt.checkpoint(
        mt_path, overwrite=True
    )  # add checkpoint to avoid repeat evaluation
    logging.info(f'Number of QC-passing, biallelic SNPs: {mt.count()[0]}')

    # filter variants found to be in promoter regions
    mt = mt.filter_rows(
        mt.vep.regulatory_feature_consequences['biotype'][0] == 'promoter'
    )
    mt_path = output_path('promoter_variants.mt', 'tmp')
    mt = mt.checkpoint(mt_path, overwrite=True)  # checkpoint
    logging.info(
        f'Number of rare variants (freq<5%) in promoter regions: {mt.count()[0]}'
    )

    # filter rare variants only (MAF < 5%)
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(
        (mt.variant_qc.AF[1] < 0.05) & (mt.variant_qc.AF[1] > 0)
        | (mt.variant_qc.AF[1] > 0.95) & (mt.variant_qc.AF[1] < 1)
    )
    mt_path = output_path('rare_variants.mt', 'tmp')
    mt = mt.checkpoint(mt_path, overwrite=True)  # checkpoint
    logging.info(f'Number of rare variants (freq<5%): {mt.count()[0]}')

    # export this as a Hail table for downstream analysis
    ht_filename = output_path(f'{gene_name}_rare_promoter_summary.ht')
    ht = mt.rows()
    ht.write(ht_filename)

    # export MT object to PLINK (promoter variants)
    mt_path = output_path(f'plink_files/{gene_name}_rare_promoter')
    export_plink(mt, mt_path, ind_id=mt.s)

    return [mt_path, ht_filename]

# endregion GET_RELEVANT_VARIANTS

# region PREPARE_INPUT_FILES

def prepare_input_files():
    """Prepare association test input files

    Input:
    genotype: plink files
    phenotype: gene expression

    Output:
    genotype matrix
    phenotype vector
    """

# endregion PREPARE_INPUT_FILES

# region RUN_ASSOCIATION

def run_gene_association():
    """Run gene-set association test

    Input:
    input files

    Output:
    list of p-values from the three tests
    """

# endregion RUN_ASSOCIATION


config = get_config()


@click.command()
@click.option('--gene-list')
def main(
    # genotypeFiles: dict[str, str], # one file per chromosome
    # genotypeFilesBims: dict[str, str],
    # genotypeFilesFams: dict[str, str],
    # phenotypeFiles: dict[str, str],
    # mafFiles: dict[str, str],
    # contextFile: str,
    # kinshipFile: str,
    # sampleMappingFile: str,
    # featureVariantFile: str,
    # nContexts: int=10,
    # fdrThreshold: int=1,
    mt_path: str,               # 'mt/v7.mt'
    vep_ht_path: str,           # 'tob_wgs_vep/104/vep104.3_GRCh38.ht'
    gene_loc_files: list[str],
):

    

    sb = hb.ServiceBackend(
    billing_project=config['hail']['billing_project'],
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

                plink_output_prefix = output_path(f'plink_files/{gene}_rare_promoter')
                plink_job_reference_and_outputs_for_gene[gene] = (job, plink_output_prefix)
                result = job.call(get_promoter_variants, gene_file=chr_loc_file, gene_name=gene, window_size=50_000, plink_output_prefix=plink_output_prefix)



    _gene = "BRCA1"
    dependent_job, output_prefix = plink_job_reference_and_outputs_for_gene[_gene]
    downstream_job = batch.new_job("Downstream job")
    downstream_job.depends_on(dependent_job)

        # run that function
            




