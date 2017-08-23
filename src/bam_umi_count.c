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

#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <limits.h>
#include <float.h>

#define HASHSIZE 1000001
#define uint_64 unsigned long long

#include "hash.h"
#include "fastq.h"

#define FEAT_ID_MAX_LEN 21
#define MAX_BARCODE_LEN 25


typedef struct count_ENTRY {
  ulong key;
  char feat_id[FEAT_ID_MAX_LEN];
  uint_64 cell;
  uint_64 umi;
  uint_64 sample;
  float umi_obs;
  float reads_obs;
} COUNT_ENTRY;

typedef struct uniq_keys {
  char *feat_id;
  uint_64 cell;
  uint_64 sample;
  struct uniq_keys* next;
} UNIQ_KEYS;

// this can be optimized
const char INT2NT[]={' ','A','C','G','T','N','.'};

short base2int(char b) {
  switch(b) {
  case 'A': return(1);
  case 'C': return(2);
  case 'G': return(3);
  case 'T': return(4);
  case 'N': return(5);
  case 'a': return(1);
  case 'c': return(2);
  case 'g': return(3);
  case 't': return(4);
  case 'n': return(5);
  default:
    return(0);
  }
}


char* uint_642char(const uint_64 i,char *s) {
  
  uint_64 i2=i;
  unsigned short x=0;
  while ( i2>0 ) {
    char c=INT2NT[i2%10];
    if ( c == 0 ) 
      break;
    i2=i2/10;
    s[x++]=c;
  }
  s[x]='\0';
  //fprintf(stderr,":::: %lu (%d)->%s\n",i,x,s);
  return(s);
}

// convert a barcode (maximum 20 bases long) to integer
uint_64 char2uint_64(const char* s) {
  uint_64 i=0;
  int pos=0;
  
  if ( s==NULL ) return(i);
  //fprintf(stderr,"s=%s\n",s);
  while ( s[pos] != '\0' && s[pos]!='\n' ) ++pos;
  // reverse so that the conversion is faster
  if ( !pos ) return(i); // nothing todo
  --pos;
  while ( pos>=0 ) {
    short base=base2int(s[pos]);
    if ( !base  ) break;
    i=i*10+base;
    --pos;
  }
#ifdef DEBUG  
  char test[30];
  char *res=uint_642char(i,&test[0]);
  fprintf(stderr,"%s<-%llu<<\n",res,i);
  if ( strcmp(res,s) ) {
    fprintf(stderr,"Internal error\nline 120: mismatch %s->%llu<<\n",res,i);
    exit(1);
  }
#endif
  return(i);
}

//
static ulong hash_uniq_umi_key(const char* feat_id,const uint_64 cell,const uint_64 sample) {

  ulong hash = 0L;
  int c;
  const char *str;
  str=feat_id;
  while ((c = *str++))
    hash = c + (hash << 6) + (hash << 16) - hash;

  hash+=(ulong)cell+(ulong)sample;
  return(hash);
}


static ulong hash_countkey(const char* feat_id,const uint_64 umi,const uint_64 cell,const uint_64 sample) {

  ulong hash = 0L;
  int c;
  const char *str;
  str=feat_id;
  while ((c = *str++))
    hash = c + (hash << 6) + (hash << 16) - hash;

  hash+=(ulong)umi+(ulong)cell+(ulong)sample;
  return(hash);
}
//
void print_count_entry(const COUNT_ENTRY *e,const char* sep,FILE *stream,short print_header,unsigned int min_num_reads) {
  char umi[MAX_BARCODE_LEN+1];
  char cell[MAX_BARCODE_LEN+1];
  char sample[MAX_BARCODE_LEN+1];
  char *labels[]={"UMI","CELL","SAMPLE"};
  char *v[3]={NULL,NULL,NULL};
  char *vl[3]={NULL,NULL,NULL};

  int x=0;
  uint_642char(e->umi,&umi[0]);
  uint_642char(e->cell,&cell[0]);
  uint_642char(e->sample,&sample[0]);

  if ( umi[0]!='\0') { v[x]=&umi[0];vl[x]=labels[0];++x;}  
  if ( cell[0]!='\0') { v[x]=&cell[0];vl[x]=labels[1];++x;}
  if ( sample[0]!='\0') { v[x]=&sample[0];vl[x]=labels[2];++x;}

    
  if (print_header) {
    if ( x==0 )
      fprintf(stream,"%s%s%s%s%s\n","Feature",sep,"COUNT",sep,"NUM READS");
    else if (x==1) 
      fprintf(stream,"%s%s%s%s%s%s%s\n","Feature",sep,vl[0],sep,"COUNT",sep,"NUM READS");    
    else if (x==2)
      fprintf(stream,"%s%s%s%s%s%s%s%s%s\n","Feature",sep,vl[0],sep,vl[1],sep,"COUNT",sep,"NUM READS");
    else
      fprintf(stream,"%s%s%s%s%s%s%s%s%s%s%s\n","Feature",sep,vl[0],sep,vl[1],sep,vl[2],sep,"COUNT",sep,"NUM READS");
  }
  if ( e->reads_obs <min_num_reads) return;


  if ( x==0 )
    fprintf(stream,"%s%s%.1f%s%1.f\n",e->feat_id,sep,e->umi_obs,sep,e->reads_obs);
  else if (x==1) 
    fprintf(stream,"%s%s%s%s%.1f%s%.1f\n",e->feat_id,sep,v[0],sep,e->umi_obs,sep,e->reads_obs);
  else if (x==2)
    fprintf(stream,"%s%s%s%s%s%s%.1f%s%.1f\n",e->feat_id,sep,v[0],sep,v[1],sep,e->umi_obs,sep,e->reads_obs);
  else   
    fprintf(stream,"%s%s%s%s%s%s%s%s%.1f%s%.1f\n",e->feat_id,sep,v[0],sep,v[1],sep,v[2],sep,e->umi_obs,sep,e->reads_obs);  
}

void print_ukey(const UNIQ_KEYS *e) {
  char cell[MAX_BARCODE_LEN+1];
  char sample[MAX_BARCODE_LEN+1];

  uint_642char(e->cell,&cell[0]);
  uint_642char(e->sample,&sample[0]);
  fprintf(stdout,"%s %s %s\n",e->feat_id,cell,sample);
}

void print_ucount(const UNIQ_KEYS *e, const ulong n,const char* sep,FILE *stream,const short print_header) {
  char cell[MAX_BARCODE_LEN+1];
  char sample[MAX_BARCODE_LEN+1];
  char *labels[]={"CELL","SAMPLE"};
  char *v[2]={NULL,NULL};
  char *vl[2]={NULL,NULL};

  uint_642char(e->cell,&cell[0]);
  uint_642char(e->sample,&sample[0]);

  int x=0;
  if ( cell[0]!='\0') { v[x]=&cell[0];vl[x]=labels[0];++x;}
  if ( sample[0]!='\0') { v[x]=&sample[0];vl[x]=labels[1];++x;}

    
  if (print_header) {
    if (x==1) 
      fprintf(stream,"%s%s%s%s%s\n","Feature",sep,vl[0],sep,"COUNT");    
    else
      fprintf(stream,"%s%s%s%s%s%s%s\n","Feature",sep,vl[0],sep,vl[1],sep,"COUNT");
  }

  if (x==1) 
    fprintf(stream,"%s%s%s%s%lu\n",e->feat_id,sep,v[0],sep,n);
  else 
    fprintf(stream,"%s%s%s%s%s%s%lu\n",e->feat_id,sep,v[0],sep,v[1],sep,n);
}


UNIQ_KEYS* add_to_list(UNIQ_KEYS* keys, COUNT_ENTRY *e) {
  UNIQ_KEYS *new=(UNIQ_KEYS*)malloc(sizeof(UNIQ_KEYS));
  if (new==NULL) { return(NULL);}
  new->cell=e->cell;
  new->sample=e->sample;
  new->feat_id=&e->feat_id[0]; // optimization only valid if memory is not released
  new->next=keys;
  return(new);
}

COUNT_ENTRY* new_count_entry(const char* feat_id,const char* umi,const char* cell,const char* sample) {

  int len1=strlen(feat_id);
  assert ( len1+1 < FEAT_ID_MAX_LEN );

  // 64
  //fprintf(stderr,"size=%d\n",sizeof(COUNT_ENTRY));
  COUNT_ENTRY *ck=(COUNT_ENTRY*)malloc(sizeof(COUNT_ENTRY));
  if (ck==NULL) { return(NULL);}  

  strncpy(&ck->feat_id[0],feat_id,len1+1);
  ck->cell=char2uint_64(cell);
  ck->umi=char2uint_64(umi);
  ck->sample=char2uint_64(sample);
  ck->key=hash_countkey(ck->feat_id,ck->umi,ck->cell,ck->sample);
  ck->umi_obs=0;
  ck->reads_obs=0;
  return(ck);
}

// 
void add_count_entry_to_uht(hashtable uniq_ht,COUNT_ENTRY *e) {

  ulong key=hash_uniq_umi_key(e->feat_id,e->cell,e->sample);
  COUNT_ENTRY* obj=(COUNT_ENTRY*)get_object(uniq_ht,key);
  while (obj!=NULL) {     
    if (
	obj->key==e->key && 
	obj->cell==e->cell &&
	obj->umi==e->umi &&
	obj->sample==e->sample &&
	!strcmp(obj->feat_id,e->feat_id)
	) break;
    obj=(COUNT_ENTRY*)get_next_object(uniq_ht,key);
  }
  if ( obj == NULL ) {
    // add object
    if(insere(uniq_ht,key,e)<0) {
      fprintf(stderr,"ERROR: unable to add entry to hash table");
      exit(1);
    }
  }
  return;
}

ulong count_uniq_entries(hashtable uniq_ht,UNIQ_KEYS *key,ulong min_num_reads) {

  ulong hkey=hash_uniq_umi_key(key->feat_id,key->cell,key->sample);
  ulong uniq_entries=0;
  //print_ukey(key);
  COUNT_ENTRY* e=(COUNT_ENTRY*)get_object(uniq_ht,hkey);
  while (e!=NULL) {     
    if (
	key->cell==e->cell &&
	key->sample==e->sample &&
	!strcmp(key->feat_id,e->feat_id) &&
	e->reads_obs >= min_num_reads
	) ++uniq_entries;
    e=(COUNT_ENTRY*)get_next_object(uniq_ht,hkey);
  }
  return(uniq_entries);
}

// check if an entry exists, returns a pointer for the the existing COUNT_ENTRY if available
COUNT_ENTRY* get_count_key(const char* feat_id,const char* umi,const char* cell,const char* sample,hashtable ht) {
  uint_64 celli=char2uint_64(cell);
  uint_64 samplei=char2uint_64(sample);
  uint_64 umii=char2uint_64(umi);

  // lookup
  ulong key=hash_countkey(feat_id,umii,celli,samplei);
  COUNT_ENTRY* e=(COUNT_ENTRY*)get_object(ht,key);
  while (e!=NULL) {     
    if ( 
	 celli==e->cell &&
	 umii==e->umi &&
	 samplei==e->sample &&
	 !strcmp(feat_id,e->feat_id)
	 ) break;
    e=(COUNT_ENTRY*)get_next_object(ht,key);
  }
  if ( e!=NULL) return(e);
  // new entry
  
  e=new_count_entry(feat_id,umi,cell,sample);
  if(insere(ht,e->key,e)<0) {
    fprintf(stderr,"ERROR: unable to add entry to hash table");
    exit(1);
  }

  return(e);
}

char EMPTY_STRING[]="";
char *get_tag(bam1_t *aln,const char tagname[2]) {

  uint8_t *s=bam_aux_get(aln,tagname);
  if (s==0) return(EMPTY_STRING);

  char *s2=bam_aux2Z(s);
  if (s2==0) return(EMPTY_STRING);

  return(s2);
}

int valid_umi(hashtable ht,char *umi_s) {
  if (ht==NULL) return(TRUE); // by default all UMIs are valid
  uint_64 umi=char2uint_64(umi_s);
  ulong key;
  if ( umi> 4294967296) {
    key=umi/1000000000;// convert to ulong
  } else {
    key=umi;
  }
  uint_64 *ptr=(uint_64*)get_object(ht,key);
  while ( ptr!=NULL ) {
    if ( *ptr==umi) return(TRUE);
    ptr=(uint_64*)get_next_object(ht,key);
  }
  //fprintf(stderr,"not valid %s\n",umi_s);
  return(FALSE);
}


void print_usage(int exit_status) {
    PRINT_ERROR("Usage: bam_umi_count --bam in.bam --ucounts output_filename [--min_reads 0] [--uniq_mapped|--multi_mapped]  [--dump filename] [--tag GX|TX] [--known_umi file_one_umi_per_line]");
    if ( exit_status>=0) exit(exit_status);
}


int main(int argc, char *argv[])  
{  
  bamFile in; 
  uint min_num_reads=0;
  static int uniq_mapped_only=FALSE;
  char feat_tag[]="GX";
  unsigned long long num_alns=0;
  unsigned long long num_tags_found=0;
  UNIQ_KEYS* key_list=NULL;

  char *bam_file=NULL;
  char *ucounts_file=NULL;
  char *dump_file=NULL;
  char *known_umi_file=NULL;
  static int verbose=0;  
  static int help=FALSE;
  static struct option long_options[] = {
    {"verbose", no_argument,       &verbose, TRUE},
    {"multi_mapped", no_argument,       &uniq_mapped_only, FALSE},
    {"uniq_mapped", no_argument,       &uniq_mapped_only, TRUE},
    {"help",   no_argument, &help, TRUE},
    {"bam",  required_argument, 0, 'b'},
    {"known_umi",  required_argument, 0, 'k'},
    {"ucounts",  required_argument, 0, 'u'},
    {"dump",  required_argument, 0, 'd'},
    {"tag",  required_argument, 0, 'x'},
    {"min_reads",  required_argument, 0, 't'},
    {0,0,0,0}
  };

  fprintf(stderr,"bam_umi_count version %s\n",VERSION);
  // process arguments
  while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    int c = getopt_long (argc, argv, "b:u:d:t:x:",
		     long_options, &option_index);      
    if (c == -1) // no more options
      break;

    switch (c) {      
    case 'b':
      bam_file=optarg;
      break;
    case 'u':
      ucounts_file=optarg;
      break;
    case 'd':
      dump_file=optarg;
      break;
    case 'k':
      known_umi_file=optarg;
      break;
    case 'x':
      strncpy(feat_tag,optarg,3);
      break;
    case 't':
      min_num_reads=atol(optarg);
      break;
    default:
      print_usage(1);
      break;
    }
  }
  if ( bam_file == NULL ) print_usage(1);
  if ( ucounts_file == NULL ) print_usage(1);
  
  if (help ) {
    print_usage(0);
  }
  
  // fprintf(stderr,"HASHSIZE=%u\n",HASHSIZE);
  hashtable ht=new_hashtable(HASHSIZE);
  // uniq=feature|cell|sample 
  hashtable uniq_ht=new_hashtable(HASHSIZE);
  hashtable kumi_ht=NULL;
  
  if ( known_umi_file!=NULL ) {
    FILE *kumi_fd;
    if ((kumi_fd=fopen(known_umi_file,"r"))==NULL) {
      PRINT_ERROR("Failed to open file %s", known_umi_file);  
      exit(1);
    }
    // known UMIs
    kumi_ht=new_hashtable(20341);
    char buf[200];
    unsigned long num_read_umis=0;
    while (!feof(kumi_fd) ) {
      char *umi=fgets(&buf[0],200,kumi_fd);
      if (umi==NULL || umi[0]=='\0') continue;
      ++num_read_umis;
      uint_64 umi_num=char2uint_64(umi);
      //fprintf(stderr,"%s->%llu\n",umi,umi_num);
      ulong key;
      if ( umi_num> 4294967296) {
	key=umi_num/1000000000;// convert to ulong
      } else {
	key=umi_num;
      }
      uint_64 *ptr=(uint_64*)get_object(kumi_ht,key);
      while ( ptr!=NULL ) {
	if ( *ptr==umi_num) continue;
	ptr=(uint_64*)get_next_object(kumi_ht,key);
      }
      if ( ptr==NULL ) {
	uint_64 *e=(uint_64*)malloc(sizeof(uint_64));
	if ( e==NULL ) {
	  PRINT_ERROR("Failed to allocate memory\n");
	  exit(1);
	}
	*e=umi_num;
	insere(kumi_ht,key,e);// create a dummy entry (key is unique)
      }
    }
    fclose(kumi_fd);
    fprintf(stderr,"unique UMIs %lu\n",kumi_ht->n_entries);
  }
  
  // Open file and exit if error
  in = strcmp(bam_file, "-")? bam_open(bam_file, "rb") : bam_dopen(fileno(stdin), "rb"); 
  
  if (in == 0 ) {  
    PRINT_ERROR("Failed to open BAM file %s", bam_file);  
    return(PARAMS_ERROR_EXIT_STATUS);  
  }  

  FILE *ucounts_fd=NULL,*dump_fd=NULL;
  if ((ucounts_fd=fopen(ucounts_file,"w+"))==NULL) {
    PRINT_ERROR("Failed to open file %s", ucounts_file);  
    exit(1);
  }

  if ( dump_file!=NULL) 
    if ((dump_fd=fopen(dump_file,"w+"))==NULL) {
      PRINT_ERROR("Failed to open file %s", dump_file);  
      exit(1);
    }

  fprintf(stderr,"@min_num_reads=%u\n",min_num_reads);
  fprintf(stderr,"@uniq mapped reads=%u\n",uniq_mapped_only);
  fprintf(stderr,"@tag=%s\n",feat_tag);
  fprintf(stderr,"@unique counts file=%s\n",ucounts_file);

  //
  // 
  bam1_t *aln=bam_init1();

  fprintf(stderr,"Processing %s\n",bam_file);

  // read header
  bam_header_read(in);

  uint8_t *nh;
  char *feat,*umi,*cell,*sample;
  //
  num_alns=0;
  while(bam_read1(in,aln)>=0) { // read alignment
    if ( num_alns == ULLONG_MAX ) {
      PRINT_ERROR("counter overflow (number of alignments) - %llu\n",num_alns);
      exit(3);
    }
    ++num_alns;
    if (num_alns%100000==0) { fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%llu",num_alns); fflush(stderr); }

    if (aln->core.tid < 0) continue;//ignore unaligned reads
    if (aln->core.flag & BAM_FUNMAP) continue;
    if (aln->core.flag & BAM_FSECONDARY) continue; // use only primary alignments
    if (aln->core.flag & BAM_FPAIRED & BAM_FPROPER_PAIR & BAM_FREAD2 ) continue; // avoid double counting

    //assert(r!=NULL);
    nh = bam_aux_get(aln, "NH");
    // exclude multimaps?
    if (nh!=NULL && uniq_mapped_only) {
      int x=bam_aux2i(nh);
      if ( x > 1 ) continue;
    }
    //int is_paired=aln->core.flag & BAM_FPAIRED;

    //fprintf(stderr,"aaaa1 %s\n",feat_tag);
    feat=get_tag(aln,feat_tag);
    if (feat[0]!='\0') {
      num_tags_found++;
      umi=get_tag(aln,"UM");
      cell=get_tag(aln,"CR");
      sample=get_tag(aln,"BC");

#ifdef DEBUG      
      fprintf(stderr,"umi2-->%s %s %s %s",umi,cell,sample,feat);
#endif
      // skip if the UMI is not valid
      if ( !valid_umi(kumi_ht,umi) ) {
        fprintf(stderr,"skipping %s\n",umi);
        continue;
      }
      //
      // feature id
      char *f=strtok(feat,",");
      char *prev_f=NULL;
      int n_feat=0;
      while (f !=NULL ) {
	if (prev_f==NULL || !strcmp(f,prev_f)) 
	  ++n_feat;
	prev_f=f;
	f=strtok(NULL,",");
      }
      
      f=strtok(feat,",");
      prev_f=NULL;
      while (f !=NULL ) {
	if (prev_f==NULL || !strcmp(f,prev_f)) {
	  COUNT_ENTRY* lookup=get_count_key(f,umi,cell,sample,ht);
	  lookup->reads_obs+=1/n_feat; // for now it will be 0 if multiple feat. overlap
	  lookup->umi_obs+=1/n_feat; // for now it will be 0 if multiple feat. overlap
	  // check for overflow
	  if ( lookup->reads_obs + 1 >= FLT_MAX ) {
	    PRINT_ERROR("Counter overflow (umi_obs) %.2lf\n",lookup->reads_obs);	    
	    exit(1);
	  }
	  if ( lookup->umi_obs + 1 >= FLT_MAX ) {
	    PRINT_ERROR("Counter overflow (umi_obs) %.2lf\n",lookup->umi_obs);	    
	    exit(1);
	  }
	  // add the entry to the other unique table
	  add_count_entry_to_uht(uniq_ht,lookup);
	  key_list=add_to_list(key_list,lookup);
	  if ( key_list==NULL ) {
	    PRINT_ERROR("Error: failed to allocate memory\n");
	    exit(1);
	  }
	}
	prev_f=f;
	f=strtok(NULL,",");
      }
    }
  }
  fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\n");fflush(stderr);
  
  bam_destroy1(aln);
  // write output
  fprintf(stderr,"Alignments processed: %llu\n",num_alns);
  fprintf(stderr,"%s encountered  %llu times\n",feat_tag,num_tags_found);
  if ( !num_tags_found ) {
    fprintf(stderr,"ERROR: no valid alignments tagged with %s were found in %s.\n",feat_tag,bam_file);
    exit(1);
  }
  // uniq UMIs that overlap each gene per cell (and optionally per sample)
  // traverse the hashtable and count how many distinct UMIs are assigned to each gene (with the number of reads above minimum number of reads threshold)
  if ( ucounts_file !=NULL) {
    UNIQ_KEYS *ukey=key_list;
    int pheader=TRUE;  
    while ( ukey!=NULL ) {
      ulong n=count_uniq_entries(uniq_ht,ukey,min_num_reads);
      if ( n>=min_num_reads)
	print_ucount(ukey,n,"\t",ucounts_fd,pheader);
      pheader&=FALSE;
      ukey=ukey->next;
    }
    fclose(ucounts_fd);
  }
  // dump the counts
  if ( dump_file != NULL ) {
    init_hash_traversal(ht);
    COUNT_ENTRY* e;
    int pheader=TRUE;  
    while((e=(COUNT_ENTRY*)next_hash_object(ht))!=NULL) {
      print_count_entry(e,"\t",dump_fd,pheader,min_num_reads);
      pheader&=FALSE;
    }
    fclose(dump_fd);
  }
  return(0);
}
