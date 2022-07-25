import os
import click
import pandas as pd
import numpy as np

def smartAppend(table,name,value):
    """ helper function for appending in a dictionary """
    if name not in table.keys():
        table[name] = []
    table[name].append(value)

def qv_estimate(pv, m=None, verbose=False, lowmem=False, pi0=None):
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

    import scipy as sp
    from scipy import interpolate

    assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"

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
            pi0.append(counts[l]/(m*(1-lam[l])))

        pi0 = sp.array(pi0)

        # fit natural cubic spline
        tck = interpolate.splrep(lam, pi0, k=3)
        pi0 = interpolate.splev(lam[-1], tck)
        if verbose:
            print("qvalues pi0=%.3f, estimated proportion of null features " % pi0)

        if pi0 > 1:
            if verbose:
                print("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)
            pi0 = 1.0

    assert(pi0 >= 0 and pi0 <= 1), "pi0 is not between 0 and 1: %f" % pi0

    if lowmem:
        # low memory version, only uses 1 pv and 1 qv matrices
        qv = sp.zeros((len(pv),))
        last_pv = pv.argmax()
        qv[last_pv] = (pi0*pv[last_pv]*m)/float(m)
        pv[last_pv] = -sp.inf
        prev_qv = last_pv
        for i in range(int(len(pv))-2, -1, -1):
            cur_max = pv.argmax()
            qv_i = (pi0*m*pv[cur_max]/float(i+1))
            pv[cur_max] = -sp.inf
            qv_i1 = prev_qv
            qv[cur_max] = min(qv_i, qv_i1)
            prev_qv = qv[cur_max]

    else:
        p_ordered = sp.argsort(pv)
        pv = pv[p_ordered]
        qv = pi0 * m/len(pv) * pv
        qv[-1] = min(qv[-1], 1.0)

        for i in range(len(pv)-2, -1, -1):
            qv[i] = min(pi0*m*pv[i]/(i+1.0), qv[i+1])

        # reorder qvalues
        qv_temp = qv.copy()
        qv = sp.zeros_like(qv)
        qv[p_ordered] = qv_temp

    # reshape qvalues
    qv = qv.reshape(original_shape)

    return qv

@click.command()
@click.option('--file-with-filenames', required=True)
@click.option('--fdr-threshold', required=False, default=0.05)
@click.option(
    "--output-folder", required=False, default=""
)  # by default current directory, where you are running your script from

def main(file_with_filenames, fdr_threshold, output_folder):

    x = 0
    table = {}

    with open(file_with_filenames, encoding='utf-8') as f:
        list_of_files = [line.strip() for line in f.readlines() if line.strip()]

    for file in list_of_files:
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
            temp['chrom'] = chrom
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
    outfile = "summary.csv" 
    myp = os.path.join(output_folder, outfile)
    df.to_csv(myp)

    # apply multiple testing correction (q-value)
    df["qv"] = qv_estimate(df["pv_Bonf"])
    # select only significant results (at given threshold)
    df_sign = df[df["qv"]<=fdr_threshold]
    outfile = "significant_results.csv" 
    myp = os.path.join(output_folder, outfile)
    df_sign.to_csv(myp)



if __name__ == '__main__':
    main()
