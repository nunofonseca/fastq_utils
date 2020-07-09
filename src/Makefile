# Note: avoid having boost in  -I : it may result on including the wrong regexp.h
# zlib is required
#CFLAGS+= -g -Wunused-result  -pg -Wall
CFLAGS+=  -O3 -Wunused-result  -Wall

ifdef PROFILE
CFLAGS+= -pg
endif
ifdef DO_COV_TEST
CFLAGS+= -fprofile-arcs -ftest-coverage -O0
endif
#
SHELL=bash

SAMTOOLS_PATH=../samtools-0.1.19
ZLIB_PATH=../zlib-1.2.11

ifeq ($(DEBUG),1)
CFLAGS+= -g -DDEBUG=1 -O1
endif

TARGETS=fastq_truncate fastq_filterpair  fastq_filter_n fastq_num_reads fastq_not_empty fastq_pre_barcodes fastq_trim_poly_at fastq_split_interleaved fastq_tests

ifdef SAMTOOLS_PATH
TARGETS+= bam_add_tags bam_umi_count bam2fastq fastq_info
endif

############################################################################
all: $(TARGETS)

install: $(TARGETS)
	cp $^ ../bin


fastq_filterpair: hash.o fastq_filterpair.o fastq.o 
	gcc  $(CFLAGS) $^ -lz -o $@

fastq_info:  hash.o fastq_info.o fastq.o
	gcc  $(CFLAGS) $^ -lz -o $@

fastq_filter_n: fastq_filter_n.o fastq.o hash.o
	gcc  $(CFLAGS) $^ -lz -o $@ 

fastq_num_reads: fastq_num_reads.o hash.o fastq.o
	gcc  $(CFLAGS) $^ -lz -o $@ 

fastq_not_empty: fastq_not_empty.o hash.o fastq.o
	gcc  $(CFLAGS) $^ -lz -o $@ 

fastq_truncate:  fastq_truncate.o  hash.o fastq.o
	gcc  $(CFLAGS) $^ -lz -o $@  

fastq_split_interleaved: fastq_split_interleaved.o fastq.o hash.o
	gcc  $(CFLAGS) $^ -lz -o $@ 

fastq_tests: fastq_tests.o hash.o fastq.o range_list.o
	gcc  $(CFLAGS) $^ -lz -o $@


# deprecated
#fastq_validator:  hash.o fastq_validator.o
#	gcc  $(CFLAGS) $^ -o $@

fastq_trim_poly_at: fastq_trim_poly_at.o hash.o fastq.o
	gcc  $(CFLAGS) $^ -lz -o $@

##fastq_trim_poly_at: fastq_sanger2phred.o hash.o fastq.o
##	gcc  $(CFLAGS) $^ -lz -o $@


fastq_pre_barcodes: fastq.o fastq.h fastq_pre_barcodes.o hash.o
	gcc  $(CFLAGS) $(patsubst %.h,,$^) -lz -o $@ 


# Companion of fastq preprocess barcodes fastq_pre_barcodes
bam_add_tags: hash.o bam_add_tags.o
	gcc  $(CFLAGS) $^ -L $(ZLIB_PATH) -I $(ZLIB_PATH) -L $(SAMTOOLS_PATH) -I $(SAMTOOLS_PATH) -lpthread  -lbam -lm -lz -pthread -o $@

bam_umi_count: range_list.o  hash.o  bam_umi_count.o 
	gcc  $(CFLAGS)  $^  -L $(ZLIB_PATH) -I $(ZLIB_PATH) -L $(SAMTOOLS_PATH) -I $(SAMTOOLS_PATH) -lpthread  -lbam -lm -lz -pthread -o $@


bam_umi_count_old:   hash.o  bam_umi_count_old.o 
	gcc  $(CFLAGS)  $^  -L $(ZLIB_PATH) -I $(ZLIB_PATH) -L $(SAMTOOLS_PATH) -I $(SAMTOOLS_PATH) -lpthread  -lbam -lm -lz -pthread -o $@


bam2fastq:   bam2fastq.o fastq.o hash.o
	gcc  $(CFLAGS)  $^  -L $(ZLIB_PATH) -I $(ZLIB_PATH) -L $(SAMTOOLS_PATH) -I $(SAMTOOLS_PATH) -lpthread  -lbam -lm -lz -pthread -o $@


##################################################


fastq.o: fastq.c fastq.h hash.h
	gcc $(CFLAGS) -I $(ZLIB_PATH) -lz -c $< 

hash.o: hash.c hash.h
	gcc $(CFLAGS) -c $<

range_list.o: range_list.c range_list.h
	gcc $(CFLAGS) -c $<

fastq_tests.o: fastq_tests.c
	gcc $(CFLAGS) -c $<

%.o: %.c
	gcc -I $(SAMTOOLS_PATH) $(CFLAGS)   -c $<

fastq_pre_barcodes.o: fastq_pre_barcodes.c sam_tags.h hash.h fastq.h
	gcc -I $(SAMTOOLS_PATH) $(CFLAGS)   -c $<
clean:
	rm -f *.o $(TARGETS) *~

###########################################################################

tests:
	../run_tests.sh


gcov: 
	gcov $(TARGETS) hash.o range_list.o fastq.o


