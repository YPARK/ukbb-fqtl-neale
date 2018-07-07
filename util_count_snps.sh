#!/bin/bash -f
ld=$1
rdir=$2
out=$3
ff=${rdir}/${ld}.snp-factor.gz
[ -f $ff ] || exit 1
zcat $ff | tail -n+2 | cut -f10 | awk '{ n[$1]++ } END { print length(n) }' > ${out}
exit 0

