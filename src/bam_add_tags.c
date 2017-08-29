/*
 * =========================================================
 * Copyright 2017,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
 *
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================
 */
#include <stdio.h>  
#include <math.h>
#include <bam.h>
#include <sam.h>
#include <kstring.h>      
#include <getopt.h>

#include <string.h>
#include <assert.h>

#include "fastq.h"

#define MAX_FEAT_LEN 50
typedef struct trans_gene_map {
  char tx[MAX_FEAT_LEN];
  char gx[MAX_FEAT_LEN];
} TGM;
/*
 * Returns 1 on sucess, 0 otherwise.
 */
short get_barcodes(const char *s,char *sample,char *umi,char *cell,
		   int* sample_len, int* umi_len, int* cell_len) {
  cell[0]=sample[0]=umi[0]='\0';

  *cell_len=*sample_len=*umi_len=0;

  if ( s[0]!='S'  || s[1]!='T' || s[2]!='A' ||
       s[3]!='G' || s[4]!='S' || s[5]!='_' ) {
    //fprintf(stderr,"failed cond1\n");
    return(0);
  }
  //
  int idx=6;
  int z=0;
  if ( s[idx]!='C' || s[idx+1]!='E'  || s[idx+2]!='L' || s[idx+3]!='L' ||
       s[idx+4]!='=' ) {
    //fprintf(stderr,"failed cond2\n");
    return(0);
  }
  // cell
  idx=idx+5;
  while (s[idx]!='_') {
    cell[z++]=s[idx++];
  }
  cell[z]='\0';
  *cell_len=z;
  idx++;
  // umi
  if ( s[idx]!='U' || s[idx+1]!='M'  || s[idx+2]!='I' || s[idx+3]!='=') {
    //fprintf(stderr,"failed cond3\n");
    return(0);
  }
  idx=idx+4;
  z=0;
  while (s[idx]!='_') {
    umi[z++]=s[idx++];
  }
  umi[z]='\0';
  ++idx;
  *umi_len=z;
  // sample
  if ( s[idx]!='S' || s[idx+1]!='A'  || s[idx+2]!='M' || s[idx+3]!='P' ||
       s[idx+4]!='L' || s[idx+5]!='E' || s[idx+6]!='=' ) {
    return(0);
  }
  idx=idx+7;
  z=0;
  while (s[idx]!='_') {
    sample[z++]=s[idx++];
  }
  sample[z]='\0';
  *sample_len=z;
  return(1);
}

void print_usage(int error) {
  char msg[]="Usage: bam_add_tags --inbam <in.bam> --outbam <out.bam or - for stdout> [--tx] [--tx2gx map_file_gene_2_trans.tsv]";
  if (error>0 ) {
    PRINT_ERROR(msg); exit(error);
  } else {
    fprintf(stderr,"%s\n",msg);
  }
}

static ulong char2ulong(char *str) {

  ulong hash = 0;
  int c;
  while ((c = *str++))
    hash = c + (hash << 6) + (hash << 16) - hash;  
  return(hash);
}


char* get_gene( char* trans,hashtable ht) {

  // lookup
  ulong key=char2ulong(trans);
  TGM* e=(TGM*)get_object(ht,key);
  while (e!=NULL) {     
    if ( !strcmp(trans,e->tx) ) return(e->gx);
    e=(TGM*)get_next_object(ht,key);
  }
  return(NULL);
}


int main(int argc, char *argv[])   {  
  short out2stdout=0;
  bamFile in; 
  bamFile out; 
  // inbam
  // outbam
  // --tx_to_gx
  // --tx field (seq)
  char *inbam_file=NULL;
  char *outbam_file=NULL;
  char *map_file=NULL;
  hashtable t2g_ht=NULL;
  
  static int verbose=0;  
  static int help=FALSE;
  static int tx_tag;
  static struct option long_options[] = {
    {"verbose", no_argument,       &verbose, TRUE},
    {"tx", no_argument,      &tx_tag, TRUE},
    {"help",   no_argument, &help, TRUE},
    {"inbam",  required_argument, 0, 'i'},
    {"outbam",  required_argument, 0, 'o'},
    {"tx_2_gx",  required_argument, 0, 'm'},
    {0,0,0,0}
  };


  // process arguments
  while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    int c = getopt_long (argc, argv, "i:o:m:",
			 long_options, &option_index);      
    if (c == -1) // no more options
      break;
    
    switch (c) {      
    case 'i': inbam_file=optarg; break;
    case 'o': outbam_file=optarg; break;
    case 'm': map_file=optarg; break;
    default: break;
    }
  }
  if ( inbam_file == NULL ) print_usage(1);
  if ( outbam_file == NULL ) print_usage(1);
  if ( !tx_tag && map_file != NULL ) {
    PRINT_ERROR("missing  --tx when --tx_2_gx is provided\n");
    print_usage(PARAMS_ERROR_EXIT_STATUS);
  }
  if (help ) {
    print_usage(0);
    exit(0);
  }

  
  in = strcmp(inbam_file, "-")? bam_open(inbam_file, "rb") : bam_dopen(fileno(stdin), "rb"); 
  out2stdout = strcmp(outbam_file, "-")? 0 : 1; 
  out = strcmp(outbam_file, "-")? bam_open(outbam_file, "w") : bam_dopen(fileno(stdout), "w"); 
  if (in == 0 ) {  
    PRINT_ERROR("Failed to open BAM file %s", inbam_file);  
    return(PARAMS_ERROR_EXIT_STATUS);  
  }  

  if (out == 0) {  
    PRINT_ERROR("Failed to open BAM file %s", outbam_file);  
    return(PARAMS_ERROR_EXIT_STATUS);  
  }
  
  unsigned long num_map_entries=0;   
  if ( map_file!=NULL ) {
    FILE *map_fd;
    if ((map_fd=fopen(map_file,"r"))==NULL) {
      PRINT_ERROR("Failed to open file %s", map_file);  
      exit(1);
    }
    // map transcript to gene 
    t2g_ht=new_hashtable(100001);
    char buf[1000];

    while (!feof(map_fd) ) {
      char *s=fgets(&buf[0],1000,map_fd);
      if (s==NULL || s[0]=='\0') continue;
      
      TGM *e=(TGM*)malloc(sizeof(TGM));
      if ( e==NULL ) {
       PRINT_ERROR("Failed to allocate memory\n");
       exit(1);
      }
      // gene trans
      char *gx=strtok(s,"\t\n");
      char *tx=strtok(NULL,"\t\n");
      if (gx==NULL || tx==NULL ) {
	PRINT_ERROR("Failed to find the gene and transcript ids in %s\n",s);
	exit(1);
      }
      strcpy(e->tx,tx);
      strcpy(e->gx,gx);
      insere(t2g_ht,char2ulong(e->tx),e);// create a dummy entry (key is unique)
      ++num_map_entries;
    }
    fclose(map_fd);
    fprintf(stderr,"unique gene/transcript pairs %lu\n",t2g_ht->n_entries);
  }

  unsigned long num_alns=0;  
  // ***********
  // Copy header
  bam_header_t *header;
  header = bam_header_read(in);
  bam_header_write(out,header);

  // 
  bam1_t *aln=bam_init1();

  if (!out2stdout) {
    fprintf(stderr,"bam_add_tags version %s\n",VERSION);
    fprintf(stderr,"Processing %s\n",inbam_file);
  }
  //
  // place holders to keep the information
  char sample[MAX_BARCODE_LENGTH];
  char umi[MAX_BARCODE_LENGTH];
  char cell[MAX_BARCODE_LENGTH];
  int sample_len, umi_len,cell_len;
  num_alns=0;
  while(bam_read1(in,aln)>=0) {
    //if (aln->core.tid < 0) continue;//ignore unaligned reads
    //if (aln->core.flag & BAM_FUNMAP) continue; // the mate is unmapped
    ++num_alns;
    //assert(r!=NULL);
    char *qn=bam1_qname(aln);
    //fprintf(stderr,"-->%s\n",qn);
    if (get_barcodes(qn,&sample[0],&umi[0],&cell[0],&sample_len,&umi_len,&cell_len)) {

      if ( umi_len > 0 )   // UMI
	bam_aux_append(aln, "UM", 'Z', umi_len+1, (uint8_t *)umi);
      if ( cell_len > 0 )  // cellular barcode sequence as reported by the sequencer
	bam_aux_append(aln, "CR", 'Z', cell_len+1, (uint8_t *)cell); 
      if ( sample_len > 0 ) // sample index
	bam_aux_append(aln, "BC", 'Z', sample_len+1, (uint8_t *)sample);
      if (tx_tag ) {
	char *tx=NULL;
	if ( aln->core.tid >= 0 ) {
	  tx=header->target_name[aln->core.tid];
	  //fprintf(stderr,">>>>TX: %s\n",tx);
	  bam_aux_append(aln,"TX",'Z',strlen(tx)+1,(uint8_t *)tx);
	  //
	  if ( map_file != NULL ) {
	    char* gene=(char*)get_gene(tx,t2g_ht);
	    if ( gene!=NULL ) {
	      bam_aux_append(aln,"GX",'Z',strlen(gene)+1,(uint8_t *)gene);
	    } // else { exit(3); }  should this fail now?
	  }
	}
      }
    }
    bam_write1(out,aln);
  }

  //bam_close(in); 
  bam_close(out);
  bam_destroy1(aln);
  if (!out2stdout) {
    fprintf(stderr,"Processing %s complete\n",inbam_file);
  }

  return(0);
}

