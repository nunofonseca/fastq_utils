# 
## fastq utils  ![Build Status](https://travis-ci.org/nunofonseca/fastq_utils.svg?branch=master) [![License](http://img.shields.io/badge/license-GPL%203-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
Set of Linux utilities to validate and manipulate fastq files.

### Building

#### Dependencies

zlib (http://zlib.net) version 1.2.11 or latest is required to compile fastq_utils. 

#### Getting sources

Option 1: to use git to download the repository  with the entire code history, type:

    git clone https://github.com/nunofonseca/fastq_utils.git
    cd fastq_utils

Option 2: download the latest source release tarball from https://github.com/nunofonseca/fastq_utils/releases, and then from your download directory type:

    tar xzf fastq_utils-x.x.x.tar.gz
    cd fastq_utils-x.x.x

#### Generating the executables

To compile fastq_utils, type:

    make && make install

If something goes wrong, then remove the whole build subdirectory with make clean and start new with make. The executables will be installed in the bin folder.


### Fastq utilities

#### fastq_info - validates and collects information from single or paired fastq files.

Usage: fastq_info fastq_file1 [fastq_file2|pe]

Check a single fastq file

    fastq_info file1.fastq.gz

Check a interleaved fastq file

    fastq_info file.fastq.gz pe
    
Check two paired fastq files
    
     fastq_info file_1.fastq.gz file_2.fastq.gz    

#### fastq_filterpair - sorts and keeps the reads with a mate in two paired fastq files.

Usage: fastq_filterpair fastq_file1 fastq_file2 out_fastq_file1.fastq.gz out_fastq_file2.fastq.gz out_fastq_sing.fastq.gz

The reads with a mate in fastq_file1 and fastq_file2 are written, respectively, to out_fastq_file1.fastq.gz out_fastq_file2.fastq.gz. Reads without a mate (singleton) are kept in out_fastq_sing.fastq.gz.

Example

    fastq_filterpair file_1.fastq.gz file_1.fastq.gz file_matched_1.fastq.gz file_matched_1.fastq.gz file_singletons.fastq.gz


#### fastq_filter_n - discards reads with more than x% of uncalled bases (N).

Usage: fastq_filter_n [-n threshold] fastq_file

Outputs a gziped fastq content where reads with more than the maximum number of allowed uncalled bases are not included.
Threshold is the maximum percentage (ranging from 0 to 100) of bases with uncalled bases that a read can have. default value is 0, which means that a read with a single uncalled base would be discarded. 

#### fastq_num_reads - prints the number of reads in a fastq file

Usage: fastq_num_reads fastq_file

