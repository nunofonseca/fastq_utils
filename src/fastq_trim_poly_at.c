/*
# =========================================================
# Copyright 2017,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
# This file is part of fastq_utils.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# if not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================
*/
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <regex.h> 
#include <zlib.h> 
#include <getopt.h>

#include "fastq.h"

struct params_s {
  char *file;
  char *outfile;
  int min_poly_at_len;
  long min_len;
  FASTQ_BOOLEAN verbose;
};
typedef struct params_s Params;

#define PRINT_VERBOSE(p,s...) {if (p->verbose) fprintf(stderr,##s );  }


void validate_options(Params* options) {
  if (options->file==NULL) {
    PRINT_ERROR("missing input file (-file)");
    exit(PARAMS_ERROR_EXIT_STATUS);
  }

  if (options->outfile==NULL) {
    PRINT_ERROR("missing output file name (-outfile)");
    exit(PARAMS_ERROR_EXIT_STATUS);
  }

  return;
}

Params* Params_new(void) {
  Params* new=(Params*)malloc(sizeof(Params));
  if (new==NULL) {
    PRINT_ERROR("unable to allocate %ld bytes of memory",sizeof(Params));
    exit(SYS_INT_ERROR_EXIT_STATUS);
  }

  new->file=NULL;
  new->outfile=NULL;
  new->min_poly_at_len=10;
  new->min_len=10;
  return(new);
}

// return 1 if the read was trimmed, 0 otherwise
short trim_poly_at(FASTQ_ENTRY *m1,
		   const int min_poly_at_len) {
  if ( min_poly_at_len <=0 ) return(0);
  FASTQ_READ_OFFSET x,matched1=0,matched2=0;

  // 3 end
  // seq includes the \n
  for (x=get_elength(m1);x>=0;--x) {
    //PRINT_ERROR("%ld %c!!!!!!!!!!!!\n",x,m1->seq[x]);
    if ( m1->seq[x]!='N' && m1->seq[x]!='A' &&
	 m1->seq[x]!='n' && m1->seq[x]!='a' )
      break;
    ++matched1;
  }
  if ( matched1 >= min_poly_at_len ) {
    m1->read_len=m1->read_len-matched1;
    m1->seq[x+1]='\n';
    m1->seq[x+2]='\0';
    m1->qual[x+1]='\n';
    m1->qual[x+2]='\0';
    return(1);
  }
  //PRINT_ERROR("!!!!!!!!!!!!%s\n>%ld/%d\n",m1->seq,matched1,min_poly_at_len);
  // 5 end
  for (x=0;x<m1->read_len;++x) {
    if ( m1->seq[x]!='N' && m1->seq[x]!='T' &&
	 m1->seq[x]!='n' && m1->seq[x]!='t' )
      break;
    ++matched2;
  }
  //
  if ( matched2 >=min_poly_at_len ) {
    for(x=0; x<=m1->read_len-matched2; ++x) {
      m1->seq[x]=m1->seq[x+matched2];
      m1->qual[x]=m1->qual[x+matched2];
    }
    m1->read_len=m1->read_len-matched2;
    return(1);
  }
  return(0);
}

void print_usage(void) {

  char msg[]="\n\
  --help       :print the usage\n\
  --file <filename> :fastq (optional gzipped) file name \n\
  --ofile <filename> : fastq file name where the processed reads will be written \n\
  --min_poly_at_len integer     : minimum length of poly-A|T sequence to remove.\n\
  --min_len integer     : minimum read length.\n\
";
  fprintf(stdout,"usage: fastq_trim_poly_at --file1 fastq_file --outfile1 out_file [optional parameters]");
  fprintf(stdout,"%s",msg);
}
//
int main(int argc, char **argv ) {
  //
  Params *p=Params_new();
  int c;
  static int help=FALSE;
  opterr = 0;

  fastq_print_version();
  
  static struct option long_options[] = {
    /* These options set a flag. */
    {"help",   no_argument, &help, TRUE},
    /* These options don’t set a flag.
       We distinguish them by their indices. */
    {"min_poly_at_len",  required_argument, 0, 'a'},
    {"file",  required_argument, 0, 'b'},
    {"outfile",  required_argument, 0, 'c'},
    {"min_len",  required_argument, 0, 'd'}
  };

  // process arguments
  while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    c = getopt_long (argc, argv, "a:b:c:",
		     long_options, &option_index);
    if (c == -1) // no more options
      break;
      
    switch (c) {
      
    case 'a':
      p->min_poly_at_len=atol(optarg);
      break;

    case 'b':
      p->file=optarg;
      break;

    case 'c': 
      p->outfile=optarg;
      break;

    case 'd':
      p->min_len=atol(optarg);
      break;

    default:  // ignore
      break;
    }
  }
  if (help ) {
    print_usage();
    exit(0);
  }
  /* Print any remaining command line arguments (not options). */
  // 
  PRINT_INFO("Validating options...");
  validate_options(p);
  PRINT_INFO("Options OK.");
  PRINT_ERROR(">%d",p->min_poly_at_len);
  // Assumptions:
  // 1) the fastq files have been validated and
  // 2) reads have the same order in the multiple files
  unsigned long trimmed_reads=0;
  unsigned long discarded_reads=0;
  unsigned long processed_reads=0;

  FASTQ_ENTRY* m=NULL;
  FASTQ_FILE* fdw=NULL;
  FASTQ_FILE* fdi=NULL;
  
  fdi=fastq_new(p->file,FALSE,"r");
  m=fastq_new_entry();
  fdw=fastq_new(p->outfile,FALSE,"w4"); 

  //
  while(!gzeof(fdi->fd) ) {
    if (fastq_read_next_entry(fdi,m)==0) break; // EOF
    //PRINT_ERROR(">%s\n",m->seq);
    ++processed_reads;
    // initialize
    if (trim_poly_at(m,p->min_poly_at_len) ) {
      	// other errors would have resulted in the exit of the program
	PRINT_VERBOSE(p,"Trimmed %s\n",m->seq);
	++trimmed_reads;
    }
    if ( m->read_len >= p->min_len) 
      fastq_write_entry(fdw,m);
    else
      ++discarded_reads;
    PRINT_READS_PROCESSED(fdi->cline/4,100000);
  }
  //   extract the info, change read name, trim the read, write
  PRINT_INFO("Reads processed: %ld",processed_reads);
  PRINT_INFO("Reads trimmed: %ld",trimmed_reads);
  PRINT_INFO("Reads discarded: %ld",discarded_reads);
  fastq_destroy(fdw);  
  exit(0);
}

