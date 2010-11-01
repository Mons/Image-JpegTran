#!/usr/bin/env bash

MODULE=`perl -ne 'print($1),exit if m{(?:version|all)_from.+?([\w/.]+)}i' Makefile.PL`;
perl=perl
$perl -v

rm -rf MANIFEST.bak Makefile.old && \
pod2text $MODULE > README && \
#$perl -i -lpne 's{^\s+$}{};s{^    ((?: {8})+)}{" "x(4+length($1)/2)}se;' README && \
$perl Makefile.PL && \
make manifest && \
make && \
TEST_AUTHOR=1 make test && \
#perl tmp/test.pl bench && \
#perl -MDevel::Leak tmp/test.pl && \
#perl -MTest::Valgrind tmp/test.pl && \
#DEBUG_MEM=1 perl tmp/test.pl && \
#TEST_AUTHOR=1 runprove 'xt/*.t' && \
make disttest && \
make dist && \
cp -f *.tar.gz dist/ && \
make clean && \
rm -rf MANIFEST.bak Makefile.old && \
echo "All is OK"
