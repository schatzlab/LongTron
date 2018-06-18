#!/bin/env python2.7
#takes long-read-derived "isoforms", one per line and sorted by chr,start,end
#and clusters them by overlap, then collapses any isoforms which share the 
#same exact junctions in the same order (ignores ends)
#and outputs a GTF of "transcripts" with scores based on read support
import sys
import re
from collections import defaultdict

#1) for clusters, just keep track of lowest/highest start/end
#from set of overlapping isoforms

#2) create a concatenated string of all raw isoforms (id:junction_list) in a cluster
#then loop through list sorted by shorted to longest to find 
#junction lists contained/matching others, collapase found containments/matches
#into one isoform keeping number of reads supporting/rids


#UPDATE: new criteria
#1) right 3' ends have to match exactly
#2) internal splices have to match (already known)
#3) left 5' ends can be contained

JX_COL=4
ID_COL=0
#semi-abitrary min. intron length requirement
#from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5435141/
MIN_INTRON_LENGTH=50

def process_cluster(c, s, e, o, cluster):
    #one read in the cluster
    if len(cluster) == 1:
        sys.stdout.write('%s\tBAM\ttranscript\t%s\t%s\t%d\t%s\t.\ttranscript_id "%s";\n' % (c,str(s),str(e),1,o,cluster[0][ID_COL]))
    #2 or more reads in the cluster
    else:
        #sort by jx length
        cluster = sorted(cluster, key=lambda x: len(x[JX_COL]))
        search_str = ';'+';'.join([x[ID_COL]+':'+x[JX_COL] for x in cluster if x[JX_COL] != ''])+';'
        rcounts = {x[ID_COL]:[1,x[3],x[5],[x[ID_COL]]] for x in cluster if x[JX_COL] != '']}

        singletons = defaultdict(lambda: [0,[]])
       
        #order of shortest to longest is important here
        for r in cluster:
            rid = r[ID_COL]
            start = r[3]
            end = r[5]
            if r[JX_COL] == '':
            #    key = str(start)+':'+str(end)
            #    singletons[key][0] += 1
            #    singletons[key][1].append(rid)
                continue
            jx_patt = re.compile(r[JX_COL])
            matching = False
            #find all other rids which this matches to (even as a substring)
            for match in re.finditer(jx_patt, search_str):
                #determine the rid of this match by parsing backwards up the search string
                i = match.start(0) - 1
                rid_ = []
                rid_started = False
                #idx 0 should be ';'
                while i > 0 and search_str[i] != ';':
                    if search_str[i] == ':':
                        rid_started = True
                        i -= 1
                        continue
                    if rid_started:
                        rid_.append(search_str[i])
                    i -= 1
                #only add 1 even if this had a larger count
                #since substrings will have already been added to this
                #one's superstring count
                rid_.reverse()
                rid_ = ''.join(rid_)
                #if not a rid, or it's invalid, or it's equal to the current one or it's already been processed
                if len(rid_) === 0 or rid_ not in rcounts or rid_ != rid or rcounts[rid_][0] == -1:
                    continue
                rcounts[rid_][0] += 1
                rcounts[rid_][3].append(rid)
                if start < rcounts[rid_][1]:
                    rcounts[rid_][1] = start
                if end > rcounts[rid_][2]:
                    rcounts[rid_][2] = end
                matching = True
            #remove this one from consideration
            #if this one didn't match any others, we're fine printing it out
            if not matching:
                (count, start, end, rids) = rcounts[rid]
                sys.stdout.write('%s\tBAM\ttranscript\t%d\t%d\t%d\t%s\t.\ttranscript_id "%s";\n' % (c,start,end,count,o,';'.join(rids)))
            rcounts[rid][0] = -1
        #TODO currently skipping singleton exons, do we need to fold singletons into existing isoforms if they're compatible?
        #probably not going to have too many of these
        #for (k,v) in singletons.iteritems():
        #    (start, end) = k.split(':')
        #    sys.stdout.write('%s\tBAM\ttranscript\t%s\t%s\t%d\t%s\t.\ttranscript_id "%s";\n' % (c,start,end,v[0],o,';'.join(v[1])))


#removes short junctions which the aligner reports as introns but are probability alignment artifacts            
def remove_short_introns(jxs):
    filtered = []
    jxs = jxs.split(',')
    for jx in jxs:
        (st,en)=jx.split('-')
        if (int(en) - int(st))+1 >= MIN_INTRON_LENGTH:
            filtered.append(jx)
    return ','.join(filtered)


def main():
    (pc,ps,pe,po) = (None, None, None, None)
    cluster = []
    idx = 1
    #input needs to be sorted on chrm + end + start:
    #sort -k1,1 -k3,3n -k2,2n
    fin = open("h1k.plus","rb")
    for line in fin:
        (rid, c, o, s, jxs, e) = line.rstrip().split('\t')
        jxs = remove_short_introns(jxs)
        rid = str(idx)
        idx += 1
        s = int(s)
        e = int(e)
        #if reverse, swap coordinates
        #since we're sorting on the 3' end
        if o == '-':
            t = s
            s = e
            e = t
        #overlaps this cluster
        if pc is not None:
            #same chrm, same strand, and same end
            if pc == c and o == po and e == pe:
                if s < ps:
                    ps = s
                cluster.append([rid, c, o, s, jxs, e])
                continue
            #doesn't overlap, process old cluster
            else:
                process_cluster(pc, ps, pe, po, cluster)

        #new cluster
        pc = c
        ps = s
        pe = e
        po = o
        cluster = [[rid, c, o, s, jxs, e]]
    fin.close()
    if pc is not None:
        process_cluster(pc, ps, pe, po, cluster)


if __name__ == '__main__':
    main()
        
#ahocorasick cruft
###old thoughts which didn't pan out
#2) for collapsing, use ahocorasick:
#   a) concatenate rid1:intron_str1|rid2:intron_str2|...
#   b) sort list of jx strings by length (shortist first)
#   c) loop through sorted list tracking collapsed read counts
#       in a hash, unless the read was collapsed into another read already

#just reduce isoform to ordered set of junction coordinates string
#and do a hash keyed by the string but containing longest start/ends and read count
'''
        import ahocorasick
        acs = ahocorasick.Automaton()
        #search_str = '|'.join([x[0]+":"+x[JX_COL] for x in cluster])
        [acs.add_word(x[JX_COL], x[ID_COL]) for x in cluster]
        matches = {r[1]:r[0] for r in acs.iter(search_str)}
        #search_patt = re.compile(search_str)
        #for (pos, rid) in acs.iter(search_str):
        #    #skip over ':'
                

        for r in cluster:
            rid = r[ID_COL]
            if rid in matches:
                pos = matches[rid]
                i = pos + 1
                rid_ = ""
                while i < len(search_str) and search_str[i] != ';':
                    rid_ += search_str[i]
                if rid_ != rid:

            
            #matching = acs.iter('|'+r[JX_COL]+':'+x[ID_COL])
            #for (pos, rid) in acs.iter(r[JX_COL]):


            #search automaton for string ocurrences
            #update rcounts for any which match other than this one
            #and where they're still present
'''
