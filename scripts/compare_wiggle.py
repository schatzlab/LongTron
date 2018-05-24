#!/bin/env python2.7
import sys
import pandas as p
import time

#static columns
cols={}
#long read set
#cols['lr']=[0,1]
#short
cols['short']=[2,6]
#refseq
cols['refseq']=[7,11]
#gtex
cols['gtex']=[12,16]
#sra
cols['sra']=[17,21]

filters = [1,2,5,7,10,15,20,50,100,500,1000]

def main():
    counts_file = sys.argv[1]
    c = p.read_csv(counts_file, sep='\t')
    lr = p.read_csv('counts/lr.filtered', header=None)
    lr_key='lr_nr'
    header = ["lr_min_reads","target_min_reads","matches","all_lr_jxs","all_target_jxs","lr_recall","target_recall"]

    for d in cols.keys():
        target = p.read_csv('counts/%s.filtered' % d, header=None)
        key = '%s_nr' % d
        (c1,c2) = cols[d]
        fout = None
        if d != 'refseq':
            fout = open("%s_exact_filtered/wiggle_compare.tsv" % d, "wb")
        uout = open("%s_exact_unfiltered/wiggle_compare.tsv" % d, "wb")

        for (f, of) in [['filtered',fout],['unfiltered',uout]]:
            if f == 'filtered' and d == 'refseq':
                continue
            fcounts = p.DataFrame([c.loc[(c[lr_key] >= n) & (c[key] >= n)].count()[1] for n in filters])
            target_filters = filters
            if f == 'unfiltered':
                target_filters = [1]*len(filters)
                fcounts = p.DataFrame([c.loc[(c[lr_key] >= n) & (c[key] >= 1)].count()[1] for n in filters])
            lr_ratios = fcounts/lr
            target_ratios = fcounts/target
            out = p.DataFrame([filters, target_filters, fcounts[0].tolist(), lr[0].tolist(), target[0].tolist(), lr_ratios[0].tolist(), target_ratios[0].tolist()])
            out = out.transpose()
            out.to_csv(of, sep='\t', index=False, header=(header))
       
        if d != 'refseq':
            fout.close()
        uout.close()
        

if __name__ == '__main__':
    main()

