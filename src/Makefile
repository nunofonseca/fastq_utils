# Note: avoid having boost in  -I : it may result on including the wrong regexp.h
# zlib is required
CFLAGS=  -Wall -O3
#-g -Wunused-result 
#
SHELL=bash

ifdef DEBUG
CFLAGS+= -g -DDEBUG=1
endif

TARGETS=fastq_truncate fastq_filterpair fastq_randsample fastq_filter_n fastq_num_reads fastq_info fastq_not_empty

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


# deprecated
fastq_validator:  hash.o fastq_validator.o
	gcc  $(CFLAGS) $^ -o $@


fastq.o: fastq.c fastq.h hash.h
	gcc $(CFLAGS) -c $< -lz

hash.o: hash.c hash.h
	gcc $(CFLAGS) -c $<

%.o: %.c
	gcc $(CFLAGS) -c $<


fastq_randsample:  fastq_randsample.o 
	gcc  $(CFLAGS) -o $@ $^

clean:
	rm -f *.o $(TARGETS) *~

###########################################################################

tests:
	../run_tests.sh

