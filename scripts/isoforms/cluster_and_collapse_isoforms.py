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

def process_cluster_by_junction(c, s, e, o, cluster):
    #one read in the cluster
    if len(cluster) == 1:
        #swap coords back
        if o == '-':
            t = s
            s = e
            e = t
        sys.stdout.write('%s\tBAM\ttranscript\t%s\t%s\t%d\t%s\t.\ttranscript_id "%s";\n' % (c,str(s),str(e),1,o,cluster[0][ID_COL]))
    #2 or more reads in the cluster
    else:
        #sort by jx length
        cluster = sorted(cluster, key=lambda x: len(x[JX_COL]))
        search_str = ';'+';'.join([x[ID_COL]+':'+x[JX_COL] for x in cluster if x[JX_COL] != ''])+';'
        rcounts = {x[ID_COL]:[1,x[3],x[5],[x[ID_COL]]] for x in cluster if x[JX_COL] != ''}

        singletons = defaultdict(lambda: [0,[]])
       
        #order of shortest to longest is important here
        for r in cluster:
            rid = r[ID_COL]
            start = r[3]
            end = r[5]
            if r[JX_COL] == '':
                key = str(start)+':'+str(end)
                singletons[key][0] += 1
                singletons[key][1].append(rid)
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
                if len(rid_) == 0 or rid_ not in rcounts or rid_ == rid or rcounts[rid_][0] == -1:
                    continue
                rcounts[rid_][0] += 1
                rcounts[rid_][3].append(rid)
                if (start < rcounts[rid_][1] and rcounts[rid_][2] == '+') or \
                    (start > rcounts[rid_][1] and rcounts[rid_][2] == '-'):
                    rcounts[rid_][1] = start
                #if end > rcounts[rid_][2]:
                #    rcounts[rid_][2] = end
                matching = True
            #remove this one from consideration
            #if this one didn't match any others, we're fine printing it out
            if not matching:
                (count, start, end, rids) = rcounts[rid]
                #swap coords back
                if o == '-':
                    t = start
                    start = end
                    end = t
                sys.stdout.write('%s\tBAM\ttranscript\t%d\t%d\t%d\t%s\t.\ttranscript_id "%s";\n' % (c,start,end,count,o,';'.join(rids)))
            rcounts[rid][0] = -1
        #TODO currently skipping singleton exons, do we need to fold singletons into existing isoforms if they're compatible?
        #probably not going to have too many of these
        #for (k,v) in singletons.iteritems():
        #    (start, end) = k.split(':')
        #    sys.stdout.write('%s\tBAM\ttranscript\t%s\t%s\t%d\t%s\t.\ttranscript_id "%s";\n' % (c,start,end,v[0],o,';'.join(v[1])))

def process_cluster_simple(c, s, e, o, cluster):
    if o == '-':
        t = s
        s = e
        e = t
    all_jxs = cluster[0][JX_COL]
    rids = cluster[0][ID_COL]
    #one read in the cluster
    score = len(cluster)
    if score > 1:
        rids = [x[ID_COL] for x in cluster]
    #2 or more reads in the cluster
        #get all jxs sorted by start across all reads
        all_jxs = [[rids[idx],[str(y[3]),str(y[3])]] for (idx,y) in enumerate(cluster)]
        all_jxs.extend([[rids[idx],x.split('-')] for (idx,y) in enumerate(cluster) for x in y[JX_COL].split(',') if len(y[JX_COL]) > 0])
        #all_jxs = sorted([[idx,x.split('-')] for (idx,y) in enumerate(cluster) for x in y[JX_COL].split(',') if len(y[JX_COL]) > 0], key=lambda x: x[1][0])
        all_jxs = sorted(all_jxs, key=lambda x: int(x[1][0]))
        all_jxs = ','.join([':'.join([str(z[0]),'-'.join(z[1])]) for z in all_jxs])
        rids = ':'.join(rids)
    #sys.stdout.write('%s\tBAM\ttranscript\t%s\t%s\t%d\t%s\t.\ttranscript_id "%s";\t%s\n' % (c,str(s),str(e),1,o,cluster[0][ID_COL],all_jxs))
    sys.stdout.write('%s\tBAM\ttranscript\t%s\t%s\t%d\t%s\t.\ttranscript_id "%s";\t%s\n' % (c,str(s),str(e),score,o,rids,all_jxs))

process_cluster = process_cluster_simple

#removes short junctions which the aligner reports as introns but are probably alignment artifacts
def remove_short_introns(jxs):
    if len(jxs) == 0:
        return ''
    filtered = []
    #sys.stderr.write(jxs+"\n")
    jxs = jxs.split(',')
    jxs.pop()
    for jx in jxs:
        (st,en)=jx.split('-')
        if (int(en) - int(st))+1 >= MIN_INTRON_LENGTH:
            filtered.append(jx)
    return ','.join(filtered)


def main():
    (pc,ps,pe,po) = (None, None, None, None)
    cluster = []
    idx = 0
    #input needs to be sorted on chrm + end + start:
    #sort -k1,1 -k3,3n -k2,2n
    fin = open(sys.argv[1],"rb")
    #fin = open("h1k.plus","rb")
    #fin = open("h1k.minus","rb")
    for line in fin:
        (rid, c, o, s, jxs, e) = line.rstrip().split('\t')
        idx += 1
        jxs = remove_short_introns(jxs)
        if len(jxs) == 0:
            continue
        rid = str(idx)
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
            #if pc == c and o == po and s <= pe:
                if (s < ps and o == '+') or (s > ps and o == '-'):
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
