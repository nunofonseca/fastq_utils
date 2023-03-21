# 
## fastq utils   
![Build Status](https://travis-ci.org/nunofonseca/fastq_utils.svg?branch=master) [![codecov](https://codecov.io/gh/nunofonseca/fastq_utils/branch/master/graph/badge.svg)](https://codecov.io/gh/nunofonseca/fastq_utils) [![Dockerhub](https://img.shields.io/docker/automated/jrottenberg/ffmpeg.svg)](https://hub.docker.com/r/nunofonseca/fastq_utils/tags/) [![License](http://img.shields.io/badge/license-GPL%203-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html) [![DOI](https://zenodo.org/badge/94895209.svg)](https://zenodo.org/badge/latestdoi/94895209) [![Conda](https://img.shields.io/conda/pn/bioconda/fastq_utils.svg?style=flat)](https://anaconda.org/bioconda/fastq_utils)

![conda version](https://anaconda.org/bioconda/fastq_utils/badges/version.svg) ![conda downloads](https://anaconda.org/bioconda/fastq_utils/badges/downloads.svg)



Set of Linux utilities to validate and manipulate fastq files. It also includes a set of programs to preprocess barcodes (namely UMIs, cells and samples), add the barcodes as tags in BAM files and count UMIs.

Programs:
   1. [fastq_info](#fastq_info---validates-and-collects-information-from-single-or-paired-fastq-files)
   2. [fastq_filterpair](#fastq_filterpair---sorts-and-keeps-the-reads-with-a-mate-in-two-paired-fastq-files)
   3. [fastq_filter_n](#fastq_filter_n---discards-reads-with-more-than-x-of-uncalled-bases-n)
   4. [fastq_num_reads](#fastq_num_reads---prints-the-number-of-reads-in-a-fastq-file)
   5. [fastq_trim_poly_at](#fastq_trim_poly_at---trims-poly-a-stretches-at-the-3-end-and-poly-t-at-5-end-of-each-read-optionally-discarding-reads-with-a-length-below-the-given-threshold)
   6. [fastq_pre_barcodes](#fastq_pre_barcodes---preprocess-the-reads-to-move-the-barcodes-umi-cell--to-the-respective-readname-optionally-discarding-reads-with-bases-in-the-barcode-regions-below-a-given-threshold)
   7. [bam_add_tags](#bam_add_tags---companion-program-to-fastq_pre_barcodes)
   8. [bam_umi_count](#bam_umi_count---count-the-number-of-unique-umis)
   9. [fastq2bam](#fastq2bam---lossless-fastq-to-bam-convertor)
   10. [bam2fastq](#bam2fastq---bam-to-fastq-convertor)	

### Installation

#### Conda

An easy way to install fastq_utils is via Conda. See the the [Conda documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to get started if you don't already have a working installation.

Installation of fastq_utils is then as simple as:

```
conda install -c bioconda fastq_utils
```

This handles software dependencies and provides the pre-built binaries for fastq_utils itself.

#### Building from source

Alternatively you can build the package yourself from source. 

##### Dependencies

samtools (version 0.1.19) and zlib (http://zlib.net) version 1.2.11 or latest are required to compile fastq_utils. 
The [install_deps.sh](https://github.com/nunofonseca/fastq_utils/blob/master/install_deps.sh) script in the toplevel folder tries to download and compile the dependencies. The bam_annotate.sh script requires samtools (version 1.5 or higher).

##### Getting sources

Option 1: download the latest source release tarball from https://github.com/nunofonseca/fastq_utils/releases, and then from your download directory type:

    tar xzf fastq_utils-x.x.x.tar.gz
    cd fastq_utils-x.x.x

Option 2: to use git to download the repository  with the entire code history, type:

    git clone https://github.com/nunofonseca/fastq_utils.git
    cd fastq_utils


##### Generating the executables

To compile fastq_utils, type:

    make && make install

If something goes wrong, then remove the whole build subdirectory with make clean and start new with make. The executables will be installed in the bin folder.


### Programs



#### fastq_info - validates and collects information from single or paired fastq files.

Usage: fastq_info [-s -r ] fastq_file1 [fastq_file2|pe]

If the fastq file(s) pass the checks then the program will exit with an exit status of 0 otherwise an error message is printed describing the error and the exit status will be different from 0. Further details about the checks are available in the [wiki](https://github.com/nunofonseca/fastq_utils/wiki/FASTQ-validation) page.

By using the -r option no checks are made to determine if the read names/identifiers are unique (fastq_info will run faster and use less memory). The -s option can be used when the reads are sorted in the same way in two paired fastq files. This option combined with -r for paired fastq files will make the validation checks less strict but fastq_info will run faster and use a fraction of the memory.

##### Examples

Check a single fastq file

    fastq_info file1.fastq.gz

Check a interleaved fastq file

    fastq_info file.fastq.gz pe
    
Check two paired fastq files
    
     fastq_info file_1.fastq.gz file_2.fastq.gz    

There is a script available (fastq_validator.sh) that also accepts a BAM file or bzip2 compressed .FASTQ files as input besides of .FASTQ or .FASTQ.gz files.

     fastq_validator.sh file_1.fastq.bzip2 file_2.fastq.bzip2    


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

#### fastq_trim_poly_at - trims poly-A stretches at the 3'-end and poly-T at 5'-end of each read, optionally discarding reads with a length below the given threshold.

Usage: fastq_trim_poly_at --file input_fastq_file --outfile output_fastq_file --min_poly_at_len integer --min_len integer

Example:

    fastq_trim_poly_at --file my.fastq.gz --outfile my.trimmed.fastq.gz --min_poly_at_len 10 --min_len 20

#### fastq_pre_barcodes - preprocess the reads to move the barcodes (UMI, Cell, ...) to the respective readname, optionally discarding reads with bases in the barcode regions below a given threshold.

Run fastq_pre_barcodes --help to get the full list of options.

Example:

    fastq_pre_barcodes  --read1 my.umi.fastq.gz   --outfile1 tmp.fastq.gz --phred_encoding 33 --read1_offset 22 --read1_size -1 --umi_read read1 --umi_size=8 --umi_offset 12

In the above command, the UMIs (starting in the base 12 and with a length of 8 bases) are extracted from the sequences and inserted in the respective read name. The read sequences in the output file includes the bases starting in position 22 until the end of the sequence. The modified readname will have the following format

@STAGS_CELL=[cell]_UMI=[umi]_SAMPLE=[sample]\_ETAGS\_[ORIGINAL READ NAME]

where [cell], [umi], and [sample] will have the value of the barcode (if available) and [ORIGINAL_READ_NAME] is, as the name suggest, the read name found in the input fastq file.

#### bam_add_tags - companion program to fastq_pre_barcodes. 

Given a bam file generated from fastq files preprocessed by fastq_pre_barcodes, bam_add_tags will add UM (UMI), CR (Cell), and BC (sample) tags to each alignment in the BAM file based on the information found in the respective readnames.

Usage: bam_add_tags --inbam input.bam --outbam output.bam [--tx  --tx2gx gene2trans_map.tsv --help]

#### bam_annotate.sh - annotate alignments with GX/TX tags

Given a BAM file annotates the alignments with the features defined in the bed files (genes->GX tag, transcritps->TX tag). Alignments can also be annotated as exonic or intronic if they overlap the regions defined in the exon or intron bed file respectively.

Usage: bam_annotate.sh  -b bam_file [ -e exon_bed -i intron_bed -g gene_bed -t transcript_bed ]

#### bam_umi_count - count the number of unique UMIs

Given a BAM file with the UM, CR, and BC tags (as produced by bam_add_tags) together with some extra tag. By default the bam_umi_count will count unique UMIs associated to uniquely mapped reads overlapping annotated genes. The GX tag is expected to contain the gene id. If an alignment overlaps multiple features then the UMI count will be partially (1/y) assigned to each feature. The output file (--ucounts) will contain two or more columns (tab-separated): the feature id (gene id by default); cell (if found in the BAM); sample (if found in the bam); and the respective number of unique UMIs (with at least x number of reads, where x is passed in the parameter --min_reads). Alternatively, a Matrix Market file (mtx) file is generated if --ucounts_MM option is passed. A white list of known UMIs can provided using the --known_umi option and a white list of cells provided with the --known_cells option. This is a simpler and faster approach to count UMIs - as an alternative you may want to consider the `umis count` command available in the [umis package](https://github.com/vals/umis) which will try to correct the barcodes.
 
Usage: bam_umi_count --bam in.bam --ucounts output_filename.tsv [--min_reads 0] [--uniq_mapped|--multi_mapped]  [--dump file.tsv] [--tag GX|TX]  [--known_umi file_one_umi_per_line]  [--known_cells file_one_cell_per_line] [--ucounts_MM]

#### fastq2bam - lossless fastq to bam convertor

Usage: fastq2bam -s schema -b out_bam -1 fastqfile1 [ -2 fastqfile2 -3 fastqfile3 -4 fastqfile4]

Currently supported schemas are: drop-seq, 10xV1, 10xV1i, 10xV1a, 10xV2, 10xV3, or none. This will be extended in the near future. To get more information type `fastq2bam -s selected_schema`

Examples:

    fastq2bam -s none -b test_se.bam  -1 single_end.fastq
    fastq2bam -s none -b test_pe.bam  -1 x_1.fastq.gz -2 x_2.fastq.gz
    fastq2bam -s 10xV1 -b test_10xv1.bam  -1 x_1.fastq.gz -2 x_2.fastq.gz 
    fastq2bam -s 10xV1 -b test_10xv1.bam  -1 x_1.fastq.gz -2 x_2.fastq.gz -3 x_I5.fastq.gz
    fastq2bam -s 10xV1i -b test_10xv1a.bam  -1 x_RA.fastq.gz -2 x_I1.fastq.gz
    fastq2bam -s 10xV1a -b test_10xv1a.bam  -1 x_RA_1.fastq.gz -2 x_RA_2.fastq.gz -3 x_I1.fastq.gz
    fastq2bam -s 10xV1a -b test_10xv1a.bam  -1 x_RA_1.fastq.gz -2 x_RA_2.fastq.gz -3 x_I1.fastq.gz -4 x_I2.fastq.gz		
    fastq2bam -s 10xV2 -b test_10xv2.bam  -1 x_R1.fastq.gz -2 x_R2.fastq.gz -3 x_I5.fastq.gz
    fastq2bam -s 10xV2 -b test_10xv2.bam  -1 x_R1.fastq.gz -2 x_R2.fastq.gz -3 x_I5.fastq.gz -4 x_I1.fastq.gz
    fastq2bam -s drop-seq -b test_ds.bam  -1 x_1.fastq.gz -2 x_2.fastq.gz

#### bam2fastq - bam to fastq convertor

Usage: bam2fastq --bam in.bam --out fastq_prefix [--10xV2 --10xV3]

Converts a bam file generated by fastq2bam into fastq format -
the following fastq files may be generated depending on the content of the
BAM file: fastq_prefix_1.fastq.gz, fastq_prefix_2.fastq.gz,
fastq_prefix_cell.fastq.gz, fastq_prefix_umi.fastq.gz,
fastq_prefix_sample.fastq.gz.

If `--10xV2` or `--10xV3` is passed to bam2fastq then the following files will be generated: fastq_prefix_R1.fastq.gz,  fastq_prefix_R2.fastq.gz, and, depending on contents of the BAM file, fastq_prefix_I1.fastq.gz. 


