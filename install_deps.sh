#!/bin/bash
set -e
############################################
## zlib
wget http://zlib.net/zlib-1.2.11.tar.gz
tar xvzf zlib-1.2.11.tar.gz
pushd zlib-1.2.11
./configure --prefix `pwd`/..
make
make install
popd
#############################################
#
wget https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download -O samtools-0.1.19.tar.bz2
tar -jxvf samtools-0.1.19.tar.bz2 
pushd samtools-0.1.19
make
make razip
popd
