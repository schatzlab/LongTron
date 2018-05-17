#!/bin/env python2.7
#takes long-read-derived "isoforms", one per line and sorted by chr,start,end
#and clusters them by overlap, then collapses any isoforms which share the 
#same exact junctions in the same order (ignores ends)
#and outputs a GTF of "transcripts" with scores based on read support
import sys
import ahocorasick

#1) for clusters, just keep track of lowest/highest start/end
#from set of overlapping isoforms

#2) for collapsing, use ahocorasick:
#   a) concatenate rid1:intron_str1|rid2:intron_str2|...
#   b) sort list of jx strings by length (shortist first)
#   c) loop through sorted list tracking collapsed read counts
#       in a hash, unless the read was collapsed into another read already

#just reduce isoform to ordered set of junction coordinates string
#and do a hash keyed by the string but containing longest start/ends and read count

JX_COL=4
ID_COL=0

def process_cluster(c, s, e, o, cluster):
    #one read in the cluster
    if len(cluster) == 1:
        sys.stdout.write('\t'.join(cluster[0])+"\n")
    #two reads in the cluster
    elif len(cluster) == 2:
        short_one = cluster[0][JX_COL]
        long_one = cluster[1][JX_COL]
        if len(cluster[0][JX_COL]) > len(cluster[1][JX_COL]):
            short_one = cluster[1]
        if short_one[JX_COL] in long_one[JX_COL]:
            sys.stdout.write('%s\tBAM\ttranscript\t%d\t%d\t%d\t%s\t.\ttranscript_id "%s.%s";\n' % (long_one[1],s,e,2,o,short_one[ID_COL],long_one[ID_COL]))
        else:
            sys.stdout.write('%s\tBAM\ttranscript\t%d\t%d\t%d\t%s\t.\ttranscript_id "%s";\n' % (short_one[1],short_one[2],short_one[3],1,o,short_one[ID_COL]))
            sys.stdout.write('%s\tBAM\ttranscript\t%d\t%d\t%d\t%s\t.\ttranscript_id "%s";\n' % (long_one[1],long_one[2],long_one[3],1,o,long_one[ID_COL]))
    #3 or more reads in the cluster
    else:
        #sort by jx length
        cluster = sorted(cluster, key=lambda x: len(x[JX_COL]))
        search_str = ';'.join(x[ID_COL]+':'+[x[JX_COL] for x in cluster])
        rcounts = {x[ID_COL]:[1,x[3],x[5],[x[ID_COL]]] for x in cluster}
        
        for r in cluster:
            rid = r[ID_COL]
            start = r[3]
            end = r[5]
            jx_patt = re.compile(r[JX_COL])
            matching = False
            #find all other rids which this matches to (even as a substring)
            for match in re.finditer(jx_patt, search_str):
                #determine the rid of this match
                i = match.start(0) - 1
                rid_ = None
                while i > 0:
                    if search_str[i] == ':':
                        rid_ = search_str[i-1]
                        break
                    i -= 1
                #only add 1 even if this had a larger count
                #singe substrings will have already been added to this
                #one's superstring's count
                if rid_ != rid and rcounts[rid_][0] != -1:
                    rcounts[rid_][0] += 1
                    rcounts[rid_][3].append(rid)
                    if start < rcounts[rid_][1]:
                        rcounts[rid_][1] = start
                    if end > rcounts[rid_][2]:
                        rcounts[rid_][2] = end
                    matching = True
            #remove this one from consideration
            rcounts[rid][0] = -1
            if not matching:
                (count, start, end, rids) = rcounts[rid]
                sys.stdout.write('%s\tBAM\ttranscript\t%d\t%d\t%d\t%s\t.\ttranscript_id "%s";\n' % (c,start,end,count,o,';'.join(rids)))



def main():
    (pc,ps,pe) = (None, None, None)
    cluster = []
    idx = 1
    for line in sys.stdin:
        (rid, c, o, s, jxs, e) = line.rstrip().split('\t')
        rid = idx
        idx += 1
        #overlaps this cluster
        if pc is not None:
            if s <= pe and o == po:
                if e > pe:
                    pe = e
                cluster.add([rid, c, o, s, jxs, e])
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

    if pc is not None:
        process_cluster(pc, ps, pe, po, cluster)

if __name__ == '__main__':
    main()
        
#ahocorasick cruft
'''
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
