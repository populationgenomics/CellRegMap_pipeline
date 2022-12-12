#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,import-outside-toplevel,consider-using-enumerate,chained-comparison,consider-using-dict-items

import os
import logging
from collections import defaultdict
from typing import Any, Dict, List

import click
import numpy as np
import pandas as pd

from multipy.fdr import qvalue

# use logging to print statements, display at info level
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


@click.command()
@click.option('--file-with-filenames', required=True)
@click.option('--fdr-threshold', required=False, default=0.05)
@click.option(
    '--output-folder', required=False, default=''
)  # by default current directory, where you are running your script from
def main(file_with_filenames: str, fdr_threshold: float, output_folder: str):

    table: Dict[str, List[Any]] = defaultdict(list)

    with open(file_with_filenames, encoding='utf-8') as f:
        list_of_files = [line.strip() for line in f.readlines() if line.strip()]

    for file in list_of_files:
        df = pd.read_csv(file, index_col=0)
        nsnps = int(len(df))
        if nsnps == 0:
            continue
        gene = os.path.splitext(os.path.basename(file))[0]
        chrom = df['chrom'].values[0]
        for i in range(nsnps):
            temp = {}
            temp['chrom'] = chrom
            temp['gene'] = gene
            temp['n_snps'] = nsnps
            temp['snp_id'] = df['variant'].values[i]
            temp['pv_raw'] = df['pv'].values[i]
            temp['pv_Bonf'] = nsnps * temp['pv_raw']
            if temp['pv_Bonf'] > 1:
                temp['pv_Bonf'] = 1
            if temp['pv_Bonf'] < 0:
                temp['pv_Bonf'] = 0

        for key in temp:
            table[key].append(temp[key])

    for key in table.keys():
        table[key] = np.array(table[key])

    df = pd.DataFrame.from_dict(table)
    outfile = 'summary.csv'
    myp = os.path.join(output_folder, outfile)
    df.to_csv(myp)

    # apply multiple testing correction (q-value)
    _, qvals = qvalue(df['pv_Bonf'])
    df['qv'] = list(qvals)
    # select only significant results (at given FDR threshold)
    df_sign = df[df['qv'] <= fdr_threshold]
    outfile = 'significant_results.csv'
    myp = os.path.join(output_folder, outfile)
    df_sign.to_csv(myp)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
