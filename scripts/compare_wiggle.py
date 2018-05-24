#!/bin/env python2.7
import sys
import pandas as p
import time

#static columns
cols={}
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
        target_filtered = p.read_csv('counts/%s.filtered' % d, header=None)
        target_unfiltered = p.DataFrame([target_filtered[0:1]]*len(target_filtered))
        key = '%s_nr' % d
        (c1,c2) = cols[d]
        fout = None
        if d != 'refseq':
            fout = open("%s_wiggle_filtered_compare.tsv" % d, "wb")
        uout = open("%s_wiggle_unfiltered_compare.tsv" % d, "wb")

        for (f, of) in [['filtered',fout],['unfiltered',uout]]:
            if f == 'filtered' and d == 'refseq':
                continue
            fcounts = p.DataFrame([c.loc[(c[lr_key] >= n) & (c[key] >= n)].count()[1] for n in filters])
            target_filters = filters
            target = target_filtered
            if f == 'unfiltered':
                target_filters = [1]*len(filters)
                fcounts = p.DataFrame([c.loc[(c[lr_key] >= n) & (c[key] >= 1)].count()[1] for n in filters])
                target = target_unfiltered
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

