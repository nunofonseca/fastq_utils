#!/bin/bash
set -e
TOPLEVEL_DIR=$PWD
mkdir -p deps

############################################
## zlib
wget -c https://zlib.net/fossils/zlib-1.2.13.tar.gz
tar xvzf zlib-1.2.13.tar.gz
pushd zlib-1.2.13
./configure --prefix `pwd`/..
make
make install
popd
#############################################
#
if [ ! -e deps/samtools-0.1.19.tar.bz2  ]; then
    wget -c https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download -O deps/samtools-0.1.19.tar.bz2
fi
tar jxvf deps/samtools-0.1.19.tar.bz2 
pushd samtools-0.1.19
make
make razip
popd

###############################################
# Downloads were failing often so keep a local copy to speedup the installation
echo "Installing samtools 1.x" 
if [ ! -e deps/samtools-1.9.tar.bz2 ] ; then
    wget -c https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2/download -O deps/samtools-1.9.tar.bz2
fi

tar xvjf deps/samtools-1.9.tar.bz2
pushd samtools-1.9
./configure  prefix=$TOPLEVEL_DIR
make -j $J
#prefix=$TOPLEVEL_DIR
make install prefix:=$TOPLEVEL_DIR
