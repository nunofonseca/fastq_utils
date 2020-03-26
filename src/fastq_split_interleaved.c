/*
# =========================================================
# Copyright 2012-2020,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
# This file is part of iRAP.
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

#include "fastq.h"

FASTQ_FILE* split_interleaved(char *f,char *out_prefix) {
  //unsigned long cline=1;

  fprintf(stderr,"Paired-end interleaved\n");
    
  FASTQ_FILE* fd1=fastq_new(f,FALSE,"r");
  fastq_is_pe(fd1);

  char outfile1[1024];
  char outfile2[1024];
  sprintf(&outfile1[0],"%s_1.fastq.gz",out_prefix);
  sprintf(&outfile2[0],"%s_2.fastq.gz",out_prefix);
  FASTQ_FILE* fdw1=fastq_new(&outfile1[0],FALSE,"w4");
  FASTQ_FILE* fdw2=fastq_new(&outfile2[0],FALSE,"w4"); 

  
  FASTQ_ENTRY *m1=fastq_new_entry(),
    *m2=fastq_new_entry();

  char rname1[MAX_LABEL_LENGTH];
  char rname2[MAX_LABEL_LENGTH];
  unsigned long nreads1=0;
  unsigned long len=0;    

  while(!gzeof(fd1->fd)) {
    // read 1
    if (fastq_read_entry(fd1,m1)==0) break;
    // read 2
    if (fastq_read_entry(fd1,m2)==0) {
      PRINT_ERROR("Error in file %s: line %lu: file truncated?",f,fd1->cline);
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }
    // match
    char *readname1=fastq_get_readname(fd1,m1,&rname1[0],&len,TRUE);
    char *readname2=fastq_get_readname(fd1,m2,&rname2[0],&len,TRUE);


    if ( strcmp(readname1,readname2) ) {
      PRINT_ERROR("Error in file %s: line %lu: unpaired read - %s",f,fd1->cline,readname1);
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    } 

    if (fastq_validate_entry(fd1,m1)) {
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }
    if (fastq_validate_entry(fd1,m2)) {
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }
    PRINT_READS_PROCESSED(fd1->cline/4,100000);
    nreads1+=2;
    // write to the out files
    fastq_write_entry(fdw1,m1);
    fastq_write_entry(fdw2,m2);

  }
  printf("\n");
  fastq_destroy(fdw1);
  fastq_destroy(fdw2);  
  return(fd1);
}


int main(int argc, char **argv ) {
  fastq_print_version();
    
  if (argc!=3 ) {
    PRINT_ERROR("Usage: fastq_split_interleaved interleaved_fastq out_prefix");
    exit(PARAMS_ERROR_EXIT_STATUS);
  }
  
  // interleaved    
  split_interleaved(argv[1],argv[2]);

  exit(0);
}



