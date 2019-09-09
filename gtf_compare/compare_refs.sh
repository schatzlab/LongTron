
cut -f 2 *.refmap | sort -u > refs.in.refmap.u
cut -f 2 *.tmap | sort -u > refs.in.tmap.u

comm -1 -2 refs.in.refmap.u refs.in.tmap.u > refs.in.both_maps
comm -2 -3 refs.in.refmap.u refs.in.tmap.u > refs.only_in.refmap
comm -1 -3 refs.in.refmap.u refs.in.tmap.u > refs.only_in.tmap
