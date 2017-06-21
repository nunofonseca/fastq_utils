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

