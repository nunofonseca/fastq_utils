# Note: avoid having boost in  -I : it may result on including the wrong regexp.h
# zlib is required
#CFLAGS+= -g -Wunused-result  
CFLAGS+= -Wall -O3
ifdef DO_COV_TEST
CFLAGS+= -fprofile-arcs -ftest-coverage
endif
#
SHELL=bash

SAMTOOLS_PATH=../samtools-0.1.19


ifeq ($(DEBUG),1)
CFLAGS+= -g -DDEBUG=1
endif

TARGETS=fastq_truncate fastq_filterpair  fastq_filter_n fastq_num_reads fastq_info fastq_not_empty fastq_pre_barcodes fastq_trim_poly_at fastq_split_interleaved

ifdef SAMTOOLS_PATH
TARGETS+= bam_add_tags bam_umi_count
endif

############################################################################
all: $(TARGETS)

install: $(TARGETS)
	cp $^ ../bin

fastq_filterpair: hash.o fastq_filterpair.o fastq.o 
	gcc  $(CFLAGS) $^ -o $@ -lz

fastq_info:  hash.o fastq_info.o fastq.o
	gcc  $(CFLAGS) $^ -o $@ -lz

fastq_filter_n: fastq_filter_n.o fastq.o hash.o
	gcc  $(CFLAGS) $^ -o $@  -lz

fastq_num_reads: fastq_num_reads.o hash.o fastq.o
	gcc  $(CFLAGS) $^  -o $@  -lz

fastq_not_empty: fastq_not_empty.o hash.o fastq.o
	gcc  $(CFLAGS) $^  -o $@  -lz

fastq_truncate:  fastq_truncate.o  hash.o fastq.o
	gcc  $(CFLAGS) $^  -o $@  -lz

fastq_split_interleaved: fastq_split_interleaved.o fastq.o hash.o
	gcc  $(CFLAGS) $^ -o $@  -lz

# deprecated
#fastq_validator:  hash.o fastq_validator.o
#	gcc  $(CFLAGS) $^ -o $@

fastq_trim_poly_at: fastq_trim_poly_at.o hash.o fastq.o
	gcc  $(CFLAGS) $^  -o $@  -lz

fastq_pre_barcodes: fastq.o fastq.h fastq_pre_barcodes.o hash.o
	gcc  $(CFLAGS) -g  $^ -o $@ -lz


# Companion of fastq preprocess barcodes fastq_pre_barcodes
bam_add_tags: hash.o bam_add_tags.o
	gcc  $(CFLAGS) $^ -L $(SAMTOOLS_PATH) -I $(SAMTOOLS_PATH) -lpthread  -lbam -lm -lz -pthread -o $@

bam_umi_count:   hash.o  bam_umi_count.o 
	gcc  $(CFLAGS)  $^ -L $(SAMTOOLS_PATH) -I $(SAMTOOLS_PATH) -lpthread  -lbam -lm -lz -pthread -o $@


##################################################


fastq.o: fastq.c fastq.h hash.h
	gcc $(CFLAGS) -c $< -lz

hash.o: hash.c hash.h
	gcc $(CFLAGS) -c $<

%.o: %.c
	gcc -I $(SAMTOOLS_PATH) $(CFLAGS)  -c $<

clean:
	rm -f *.o $(TARGETS) *~

###########################################################################

tests:
	../run_tests.sh

