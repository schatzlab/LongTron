cat tests/test_disjoin.tsv | perl scripts/disjoin.pl > test_disjoin.tsv.out
diff test_disjoin.tsv.out tests/test_disjoin.tsv.expected
