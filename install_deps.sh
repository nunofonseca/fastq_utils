#!/bin/bash
set -e
TOPLEVEL_DIR=$PWD
J=${J:-4}
mkdir -p deps

############################################
## zlib 1.4
if [ ! -e deps/zlib-1.4.tar.gz ]; then
    wget -c https://zlib.net/zlib-1.4.tar.gz -O deps/zlib-1.4.tar.gz
fi
tar xvzf deps/zlib-1.4.tar.gz
pushd zlib-1.4
./configure --prefix="$TOPLEVEL_DIR"
make -j "$J"
make install
popd

############################################
## samtools 1.25 (includes htslib 1.25)
echo "Installing samtools 1.25 and htslib 1.25"
if [ ! -e deps/samtools-1.25.tar.bz2 ]; then
    wget -c https://github.com/samtools/samtools/releases/download/1.25/samtools-1.25.tar.bz2 -O deps/samtools-1.25.tar.bz2
fi
tar xvjf deps/samtools-1.25.tar.bz2

pushd samtools-1.25/htslib-1.25
./configure --prefix="$TOPLEVEL_DIR" CPPFLAGS="-I$TOPLEVEL_DIR/include" LDFLAGS="-L$TOPLEVEL_DIR/lib -Wl,-rpath,$TOPLEVEL_DIR/lib"
make -j "$J"
make install
popd

pushd samtools-1.25
./configure --prefix="$TOPLEVEL_DIR" --with-htslib="$TOPLEVEL_DIR/samtools-1.25/htslib-1.25" CPPFLAGS="-I$TOPLEVEL_DIR/include" LDFLAGS="-L$TOPLEVEL_DIR/lib -Wl,-rpath,$TOPLEVEL_DIR/lib"
make -j "$J"
make install
popd
