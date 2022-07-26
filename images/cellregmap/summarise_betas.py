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
@click.option(
    "--output-folder", required=False, default=""
)  # by default current directory, where you are running your script from

def main(file_with_filenames_1, file_with_filenames_2, output_folder):

    # betaG
    table = {}

    with open(file_with_filenames_1, encoding='utf-8') as f:
        list_of_files1 = [line.strip() for line in f.readlines() if line.strip()]

    for file in list_of_files1:
        df = pd.read_csv(file, index_col=0)
        nsnps = int(len(df))
        if nsnps==0:
            continue
        filename = os.path.splitext(os.path.basename(file))[0]
        gene = filename.replace("_betaG", "")
        chrom = df['chrom'].values[0]
        for i in range(nsnps):
            temp = {}
            temp['chrom'] = chrom
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
    myp = os.path.join(output_folder, outfile)
    df.to_csv(myp)

    # betaGxC
    df_all = pd.DataFrame()
    
    with open(file_with_filenames_2, encoding='utf-8') as f:
        list_of_files2 = [line.strip() for line in f.readlines() if line.strip()]

    for file in list_of_files2:
        df = pd.read_csv(file, index_col=0)
        nsnps = int(len(df))
        if nsnps==0:
            continue
        filename = os.path.splitext(os.path.basename(file))[0]
        gene = filename.replace("_betaGxC", "")
        df.columns = gene + "_" + df.columns
        df_all = pd.concat([df_all, df], axis=1)

    outfile = "summary_betaGxC.csv" 
    myp = os.path.join(output_folder, outfile)
    df_all.to_csv(myp)

if __name__ == '__main__':
    main()
