import os
import click
import pandas as pd
import numpy as np

def smartAppend(table,name,value):
    """ helper function for appending in a dictionary """
    if name not in table.keys():
        table[name] = []
    table[name].append(value)

@click.command()
@click.option('--file-with-filenames-1', required=True)
@click.option('--file-with-filenames-2', required=True)

def main(file_with_filenames_1, file_with_filenames_2):

    # betaG
    x = 0
    table = {}

    with open(file_with_filenames_1, encoding='utf-8') as f:
        list_of_files1 = f.readlines()

    for file in list_of_files1:
        x += 1
        if x % 500 == 0: 
            print(x)
        df = pd.read_csv(file, index_col=0)
        nsnps = int(len(df))
        if nsnps==0:
            continue
        gene = os.path.splitext(os.path.basename(file))[0] 
        # print(gene)
        chrom = df['chrom'].values[0]
        # print(chrom)
        for i in range(nsnps):
            temp = {}
            temp['gene'] = gene
            temp['n_snps'] = nsnps
            temp['snp_id'] = df['variant'].values[i]
            temp['betaG'] = df['betaG'].values[i]

        for key in temp.keys():
            smartAppend(table, key, temp[key])

    print(x)
    for key in table.keys():
        table[key] = np.array(table[key])

    df = pd.DataFrame.from_dict(table)
    outfile = "summary_betaG.csv" 
    myp = os.path.join(path_results, outfile)
    df.to_csv(myp)

    # betaGxC
    x = 0
    table = {}
    
    with open(file_with_filenames_2, encoding='utf-8') as f:
        list_of_files2 = f.readlines()

    for file in list_of_files2:
        x += 1
        if x % 500 == 0: 
            print(x)
        df = pd.read_csv(file, index_col=0)
        nsnps = int(len(df))
        if nsnps==0:
            continue
        gene = os.path.splitext(os.path.basename(file))[0] 
        # print(gene)
        chrom = df['chrom'].values[0]
        # print(chrom)
        for i in range(nsnps):
            temp = {}
            temp['gene'] = gene
            temp['n_snps'] = nsnps
            temp['snp_id'] = df['variant'].values[i]
            temp['pv_raw'] = df['pv'].values[i]
            temp['pv_Bonf'] = nsnps * temp['pv_raw']
            if temp['pv_Bonf']>1: temp['pv_Bonf'] = 1
            if temp['pv_Bonf']<0: temp['pv_Bonf'] = 0

        for key in temp.keys():
            smartAppend(table, key, temp[key])

    print(x)
    for key in table.keys():
        table[key] = np.array(table[key])

    df = pd.DataFrame.from_dict(table)
    outfile = "summary_betaGxC.csv" 
    myp = os.path.join(path_results, outfile)
    df.to_csv(myp)

    if __name__ == '__main__':
    main()