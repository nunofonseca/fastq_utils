#!/bin/env bash
file=$1
ncols=$(expr `head -n 1 $file|wc -w` - 1)
nrows=$(wc -l $file|cut -f 1 -d\ )
uniq_rows=$(cut -f 1-$ncols $file|sort -u|wc -l)
if [ "$nrows-" == "$uniq_rows-" ]; then
    exit 0
else
    exit 1
fi

