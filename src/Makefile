# Note: avoid having boost in -I : it may result on including the wrong regexp.h
# zlib and htslib are required
CFLAGS+= -O3 -Wunused-result -Wall

ifdef PROFILE
CFLAGS+= -pg
endif
ifdef DO_COV_TEST
CFLAGS+= -fprofile-arcs -ftest-coverage -O0
endif

SHELL=bash
CC?=gcc

PREFIX?=..
HTS_PATH?=$(PREFIX)
ZLIB_PATH?=$(PREFIX)

ZLIB_INC=-I$(ZLIB_PATH)/include
ZLIB_LIB=-L$(ZLIB_PATH)/lib -Wl,-rpath,'$$ORIGIN/../lib' -lz
HTS_INC=-I$(HTS_PATH)/include
HTS_LIB=-L$(HTS_PATH)/lib -Wl,-rpath,'$$ORIGIN/../lib' -lhts

ifeq ($(DEBUG),1)
CFLAGS+= -g -DDEBUG=1 -O1
endif

TARGETS=fastq_truncate fastq_filterpair fastq_filter_n fastq_num_reads fastq_not_empty fastq_pre_barcodes fastq_trim_poly_at fastq_split_interleaved fastq_tests
TARGETS+= bam_add_tags bam_umi_count bam2fastq fastq_info

############################################################################
all: $(TARGETS)

install: $(TARGETS)
	cp $^ ../bin

fastq_filterpair: hash.o fastq_filterpair.o fastq.o
	$(CC) $(CFLAGS) $^ $(ZLIB_LIB) -o $@

fastq_info: hash.o fastq_info.o fastq.o
	$(CC) $(CFLAGS) $^ $(ZLIB_LIB) -o $@

fastq_filter_n: fastq_filter_n.o fastq.o hash.o
	$(CC) $(CFLAGS) $^ $(ZLIB_LIB) -o $@

fastq_num_reads: fastq_num_reads.o hash.o fastq.o
	$(CC) $(CFLAGS) $^ $(ZLIB_LIB) -o $@

fastq_not_empty: fastq_not_empty.o hash.o fastq.o
	$(CC) $(CFLAGS) $^ $(ZLIB_LIB) -o $@

fastq_truncate: fastq_truncate.o hash.o fastq.o
	$(CC) $(CFLAGS) $^ $(ZLIB_LIB) -o $@

fastq_split_interleaved: fastq_split_interleaved.o fastq.o hash.o
	$(CC) $(CFLAGS) $^ $(ZLIB_LIB) -o $@

fastq_tests: fastq_tests.o hash.o fastq.o range_list.o
	$(CC) $(CFLAGS) $^ $(ZLIB_LIB) -o $@

fastq_trim_poly_at: fastq_trim_poly_at.o hash.o fastq.o
	$(CC) $(CFLAGS) $^ $(ZLIB_LIB) -o $@

fastq_pre_barcodes: fastq.o fastq.h fastq_pre_barcodes.o hash.o
	$(CC) $(CFLAGS) $(patsubst %.h,,$^) $(ZLIB_LIB) -o $@

# Companion of fastq preprocess barcodes fastq_pre_barcodes
bam_add_tags: hash.o bam_add_tags.o
	$(CC) $(CFLAGS) $^ $(HTS_LIB) $(ZLIB_LIB) -lpthread -lm -pthread -o $@

bam_umi_count: range_list.o hash.o bam_umi_count.o
	$(CC) $(CFLAGS) $^ $(HTS_LIB) $(ZLIB_LIB) -lpthread -lm -pthread -o $@

bam2fastq: bam2fastq.o fastq.o hash.o
	$(CC) $(CFLAGS) $^ $(HTS_LIB) $(ZLIB_LIB) -lpthread -lm -pthread -o $@

##################################################

fastq.o: fastq.c fastq.h hash.h
	$(CC) $(CFLAGS) $(ZLIB_INC) -c $<

hash.o: hash.c hash.h
	$(CC) $(CFLAGS) -c $<

range_list.o: range_list.c range_list.h
	$(CC) $(CFLAGS) -c $<

fastq_tests.o: fastq_tests.c
	$(CC) $(CFLAGS) $(ZLIB_INC) -c $<

fastq_pre_barcodes.o: fastq_pre_barcodes.c sam_tags.h hash.h fastq.h
	$(CC) $(CFLAGS) $(HTS_INC) $(ZLIB_INC) -c $<

%.o: %.c
	$(CC) $(CFLAGS) $(HTS_INC) $(ZLIB_INC) -c $<

clean:
	rm -f *.o $(TARGETS) *~

###########################################################################

tests:
	../run_tests.sh

gcov:
	gcov $(TARGETS) hash.o range_list.o fastq.o
