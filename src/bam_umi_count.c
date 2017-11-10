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

#define HASHSIZE 50000000
#define uint_64 unsigned long long

#include "hash.h"
#include "fastq.h"

#define FEAT_ID_MAX_LEN 25
#define MAX_BARCODE_LEN 19


typedef struct count_ENTRY {
  uint feat_id;
  uint_64 cell;
  uint_64 umi;
  uint_64 sample;
  float umi_obs;
  float reads_obs;
} COUNT_ENTRY;

typedef struct uniq_keys {
  uint feat_id;
  uint_64 cell;
  uint_64 sample;
  struct uniq_keys* next;
} UNIQ_KEYS;


typedef struct uniq_keys_ht {
  hashtable ht;
  UNIQ_KEYS* keys;
} UNIQ_KEYS_HT;

// single label => feature
typedef struct label2id {
  uint id;
  char *label;
  struct label2id* next;
} LABEL2ID;

typedef struct labels {
  LABEL2ID *first;
  LABEL2ID *last;
  uint  ctr;
  hashtable ht;
} LABELS;

// possibly two barcode labels
typedef struct blabel2id {
  uint id;
  uint_64 label;
  uint_64 label2;
  struct blabel2id* next;
} BLABEL2ID;

typedef struct blabels {
  BLABEL2ID *first;
  BLABEL2ID *last;
  uint  ctr;
  hashtable ht;
} BLABELS;


static ulong hash_str(const char *str);
static ulong hash_uniq_umi_key(const uint feat_id,const uint_64 cell,const uint_64 sample);
char* uint_642char(const uint_64 i,char *s);

// Map labels (e.g. genes) to ids (1,...N)
LABELS* init_labels(uint hashsize) {
  LABELS *new=(LABELS*)malloc(sizeof(LABELS));
  new->first=new->last=NULL;
  new->ctr=0;
  new->ht=new_hashtable(hashsize);
  return(new);
}
// feature to id
uint_64 label_str2id(const char* lab,LABELS* lm ) {
  //
  assert(lm!=NULL);
  // lookup
  ulong ikey=hash_str(lab);
  LABEL2ID *e=(LABEL2ID*)get_object(lm->ht,ikey);
  
  // look for the match
  while (e!=NULL) {
    if ( !strcmp(e->label,lab)) break;
    e=e->next;
  }

  if ( e==NULL ) {
    uint len=strlen(lab);
    LABEL2ID *new=(LABEL2ID*)malloc(sizeof(LABEL2ID)+len+1);
    lm->ctr++;
    if ( lm->last==NULL ) {
      lm->last=new; lm->first=new;
    } else {
      lm->last->next=new; lm->last=new;
    }
    new->id=lm->ctr;
    new->next=NULL;
    new->label=(char*)&(*new)+sizeof(LABEL2ID);
    strncpy(new->label,lab,len+1);
    if(insere(lm->ht,ikey,new)<0) {
      fprintf(stderr,"ERROR: unable to add entry to hash table");
      exit(1);
    }
    e=new;
  }
  return(e->id);
}

char* label_id2str(const uint id, const LABELS *lm) {
  assert(lm!=NULL);
  LABEL2ID* l2id=lm->first;
  while(l2id!=NULL) {
    if ( l2id->id==id) return(l2id->label);
    l2id=l2id->next;
  }
  return(NULL);
}
uint label_entries(const LABELS *lm) {
  return(lm->ctr);
}

// barcode labels
// Map labels (e.g. genes) to ids (1,...N)
BLABELS* init_blabels(uint hashsize) {
  BLABELS *new=(BLABELS*)malloc(sizeof(BLABELS));
  new->first=new->last=NULL;
  new->ctr=0;
  new->ht=new_hashtable(hashsize);
  return(new);
}
// barcode based label to id
uint_64 blabel2id(const uint_64 lab,const uint_64 lab2,BLABELS* lm ) {
  //
  assert(lm!=NULL);
  // lookup
  ulong ikey=hash_uniq_umi_key(0,lab,lab2);
  BLABEL2ID *e=(BLABEL2ID*)get_object(lm->ht,ikey);
  
  // look for the match
  while (e!=NULL) {
    if ( e->label==lab && e->label2==lab2) break;
    e=e->next;
  }

  if ( e==NULL ) {
    BLABEL2ID *new=(BLABEL2ID*)malloc(sizeof(BLABEL2ID));
    lm->ctr++;
    if ( lm->last==NULL ) {
      lm->last=new; lm->first=new;
    } else {
      lm->last->next=new; lm->last=new;
    }
    new->id=lm->ctr;
    new->next=NULL;
    new->label=lab;
    new->label2=lab2;
    if(insere(lm->ht,ikey,new)<0) {
      fprintf(stderr,"ERROR: unable to add entry to hash table");
      exit(1);
    }
    e=new;
  }
  return(e->id);
}

char buf[MAX_BARCODE_LEN*2+2+1];
const char* blabel_id2str(const BLABEL2ID *e) {
  assert(e!=NULL);
  if (e->label2==0) {
    uint_642char(e->label,&buf[0]);	      
  } else {
    char buf2[MAX_BARCODE_LEN+1];
    char buf3[MAX_BARCODE_LEN+1];
    uint_642char(e->label,&buf2[0]);
    uint_642char(e->label2,&buf3[0]);	
    sprintf(&buf[0],"%s::%s",buf2,buf3);
  }
  return(&buf[0]);
}

uint blabel_entries(const BLABELS *lm) {
  return(lm->ctr);
}


#define MAPSEP "\t"
void write_map2fileL(const char* file,const char *labname,LABELS* map) {
  FILE *fd;
  char buf[300];
  sprintf(&buf[0],"%s_%s",file,labname);
  file=buf;
  if ((fd=fopen(buf,"w+"))==NULL) {
    PRINT_ERROR("Failed to open file %s for writing", buf);  
    exit(1);
  }
  if ( map->first!=NULL) {
    uint_64 ctr=0;
    LABEL2ID *e=map->first;    
    while ( e!=NULL ) {
      fprintf(fd,"%lu%s%s\n",e->id,MAPSEP,e->label);
      e=e->next;
      ++ctr;
    }
  }
  fclose(fd);
}

void write_map2fileB(const char* file,const char *labname,BLABELS* map) {
  FILE *fd;
  char buf[300];
  sprintf(&buf[0],"%s_%s",file,labname);
  file=buf;
  if ((fd=fopen(buf,"w+"))==NULL) {
    PRINT_ERROR("Failed to open file %s for writing", buf);  
    exit(1);
  }
  if ( map->first!=NULL) {
    uint_64 ctr=0;
    BLABEL2ID *e=map->first;    
    while ( e!=NULL ) {
      fprintf(fd,"%lu%s%s\n",e->id,MAPSEP,blabel_id2str(e));
      e=e->next;
      ++ctr;
    }
  }
  fclose(fd);
}


// this can be optimized
const char INT2NT[]={' ','A','C','G','T','N','.'};

// TODO: use a 3 bit encoding?
static inline   short base2int(char b) {
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
  if ( x > MAX_BARCODE_LEN ) {
    fprintf(stderr,"ERROR: barcode should be at most %u bases\n",MAX_BARCODE_LEN);
    exit(1);
  }
  //fprintf(stderr,":::: %lu (%d)->%s\n",i,x,s);
  return(s);
}

// convert a barcode (maximum 20 bases long) to integer
uint_64 char2uint_64(const char* s) {
  uint_64 i=0;
  int pos=0;
  
  if ( s==NULL || *s=='\0' ) return(i);
  //fprintf(stderr,"s=%d %s\n",s[1],s);
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

static ulong hash_str(const char *str) {
  ulong hash = 0L;
  int c;
  while ((c = *str++))
    hash = c + (hash << 6) + (hash << 16) - hash;

  return(hash);
}

//
static ulong hash_uniq_umi_key(const uint feat_id,const uint_64 cell,const uint_64 sample) {

  ulong hash = 0L;
  //hash+=feat_id+cell+sample;
  //srand((unsigned)cell*feat_id);
  //hash=rand()+(sample*10000)+feat_id;
  hash=sample+feat_id+10000+cell+1000000;
  /*  hash=cell*feat_id;
  if ( sample>0 ) {   hash*=sample;}
  hash=2166136261;
  hash=(hash*16777619) ^ cell+rand();
  hash=(hash*16777619) ^ feat_id+rand();
  if ( sample>0 ) {   hash=(hash*16777619) ^ sample+rand(); }*/
  //fprintf(stdout,"%llu %llu %lu %lu %lu\n",hash%5000000,hash,feat_id,cell,sample);
  //hash=cell;
  /* hash+=feat_id;  hash+=(hash<<10); hash^=(hash>>6); */
  /* hash+=cell;  hash+=(hash<<10); hash^=(hash>>6); */
  /* if ( sample>0 ) { hash+=sample;  hash+=(hash<<10); hash^=(hash>>6);} */
  /* hash+=(hash<<3); */
  /* hash^=(hash>>11); */
  /* hash+=(hash<<15); */
  /*  hash+=(hash<<10);
  hash^=(hash>>6);
  hash^=cell;
  hash+=(hash<<10);
  hash^=(hash>>6);
  hash+=hash ^ (ulong)sample;*/
  return(hash);
}



static ulong hash_countkey(const uint feat_id,const uint_64 umi,const uint_64 cell,const uint_64 sample) {

  ulong hash = 0L;
  hash=umi+hash_uniq_umi_key(feat_id,cell,sample);
  return(hash);
}
//
void print_count_entry(const COUNT_ENTRY *e,const char* sep,FILE *stream,LABELS* feat_map, short print_header,unsigned int min_num_reads) {
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
      fprintf(stream,"%s%s%s%s%s%s%s%s%s%s%s\n","Feature",sep,vl[0],"::",vl[1],sep,vl[2],sep,"COUNT",sep,"NUM READS");
  }
  if ( e->reads_obs <min_num_reads) return;


  if ( x==0 )
    fprintf(stream,"%s%s%.1f%s%1.f\n",label_id2str(e->feat_id,feat_map),sep,e->umi_obs,sep,e->reads_obs);
  else if (x==1) 
    fprintf(stream,"%s%s%s%s%.1f%s%.1f\n",label_id2str(e->feat_id,feat_map),sep,v[0],sep,e->umi_obs,sep,e->reads_obs);
  else if (x==2)
    fprintf(stream,"%s%s%s%s%s%s%.1f%s%.1f\n",label_id2str(e->feat_id,feat_map),sep,v[0],sep,v[1],sep,e->umi_obs,sep,e->reads_obs);
  else   
    fprintf(stream,"%s%s%s%s%s%s%s%s%.1f%s%.1f\n",label_id2str(e->feat_id,feat_map),sep,v[0],sep,v[1],"::",v[2],sep,e->umi_obs,sep,e->reads_obs);  
}

void print_ukey(const UNIQ_KEYS *e) {
  char cell[MAX_BARCODE_LEN+1];
  char sample[MAX_BARCODE_LEN+1];

  uint_642char(e->cell,&cell[0]);
  uint_642char(e->sample,&sample[0]);
  fprintf(stdout,"%lu %s %s\n",e->feat_id,cell,sample);
}

void print_ucount(const UNIQ_KEYS *e, const ulong n,const char* sep,FILE *stream,LABELS* feat_map,const short print_header) {
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

  if (print_header==TRUE) {
    if (x==0)
      fprintf(stream,"%s%s%s\n","Feature",sep,"COUNT");    
    else if (x==1) 
      fprintf(stream,"%s%s%s%s%s\n","Feature",sep,vl[0],sep,"COUNT");    
    else
      fprintf(stream,"%s%s%s%s%s%s%s\n","Feature",sep,vl[0],"::",vl[1],sep,"COUNT");
  }

  if (x==0)
    fprintf(stream,"%s%s%lu\n",label_id2str(e->feat_id,feat_map),sep,n);
  else if (x==1) 
    fprintf(stream,"%s%s%s%s%lu\n",label_id2str(e->feat_id,feat_map),sep,v[0],sep,n);
  else 
    fprintf(stream,"%s%s%s%s%s%s%lu\n",label_id2str(e->feat_id,feat_map),sep,v[0],"::",v[1],sep,n);
}


UNIQ_KEYS_HT* add_to_list(UNIQ_KEYS_HT* keys, COUNT_ENTRY *e) {
  if (keys==NULL) {
    UNIQ_KEYS_HT *new=(UNIQ_KEYS_HT*)malloc(sizeof(UNIQ_KEYS_HT));
    new->ht=new_hashtable(1000001);
    new->keys=NULL;
    keys=new;
  }
  ulong ikey=hash_countkey(e->feat_id,0,e->cell,e->sample);
  UNIQ_KEYS *f=(UNIQ_KEYS*)get_object(keys->ht,ikey);
  while (f!=NULL) {
    if ( 
	 f->cell==e->cell &&
	 f->sample==e->sample &&
	 f->feat_id==e->feat_id
	 ) return(keys);
    f=f->next;
  }
  UNIQ_KEYS *new=(UNIQ_KEYS*)malloc(sizeof(UNIQ_KEYS));
  if (new==NULL) { return(NULL);}
  new->cell=e->cell;
  new->sample=e->sample;
  new->feat_id=e->feat_id;

  new->next=keys->keys;
  keys->keys=new;
  if(insere(keys->ht,ikey,new)<0) {
    fprintf(stderr,"ERROR: unable to add entry to hash table");
    exit(1);
  }
  return(keys);
}

COUNT_ENTRY* new_count_entry(const uint feat_id,const uint_64 umi,const uint_64 cell,const uint_64 sample) {

  // 64
  //fprintf(stderr,"size=%d\n",sizeof(COUNT_ENTRY));
  COUNT_ENTRY *ck=(COUNT_ENTRY*)malloc(sizeof(COUNT_ENTRY));
  if (ck==NULL) { return(NULL);}  

  ck->feat_id=feat_id;
  ck->cell=cell;
  ck->umi=umi;
  ck->sample=sample;
  ck->umi_obs=0;
  ck->reads_obs=0;
  return(ck);
}

// return 1 if entry exists, 0 if a new entry is created
hashtable add_count_entry_to_uht(hashtable uniq_ht,COUNT_ENTRY *e) {

  ulong key=hash_uniq_umi_key(e->feat_id,e->cell,e->sample);
  COUNT_ENTRY* obj=(COUNT_ENTRY*)get_object(uniq_ht,key);
  while (obj!=NULL) {     
    if (
	obj->umi==e->umi &&
	obj->cell==e->cell &&
	obj->sample==e->sample &&
	obj->feat_id==e->feat_id
	) return(uniq_ht);
    obj=(COUNT_ENTRY*)get_next_object(uniq_ht,key);
  }

  if ( obj == NULL ) {
    // add object
    if(insere(uniq_ht,key,e)<0) {
      fprintf(stderr,"ERROR: unable to add entry to hash table");
      exit(1);
    }
    //bin/fprintf(stderr,"miss\n");
    return(uniq_ht);
  }
  //fprintf(stderr,"hit\n");
  return(uniq_ht);
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
	key->feat_id==e->feat_id &&
	e->reads_obs >= min_num_reads
	) ++uniq_entries;
    e=(COUNT_ENTRY*)get_next_object(uniq_ht,hkey);
  }
  return(uniq_entries);
}

// check if an entry exists, returns a pointer for the the existing COUNT_ENTRY if available
COUNT_ENTRY* get_count_key(const uint feat_id,const uint_64 umii,const uint_64 celli,const uint_64 samplei,hashtable ht) {

  uint feat_i=feat_id;

  // lookup
  ulong key=hash_countkey(feat_id,umii,celli,samplei);
  COUNT_ENTRY* e=(COUNT_ENTRY*)get_object(ht,key);
  while (e!=NULL) {     
    if (
	umii==e->umi &&
	celli==e->cell &&
	samplei==e->sample &&
	feat_id == e->feat_id 
	) break;
    e=(COUNT_ENTRY*)get_next_object(ht,key);
  }
  if ( e!=NULL) return(e);
  // new entry
  
  e=new_count_entry(feat_id,umii,celli,samplei);
  if(insere(ht,key,e)<0) {
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

int valid_barcode(hashtable ht,const uint_64 barcode) {
  if (ht==NULL) return(TRUE); // by default all BARCODEs are valid

  ulong key;
  //if ( barcode> 4294967296) {
  //  key=barcode/1000000000;// convert to ulong
  //} else {
    key=barcode;
    //}
  uint_64 *ptr=(uint_64*)get_object(ht,key);
  while ( ptr!=NULL ) {
    if ( *ptr==barcode) return(TRUE);
    ptr=(uint_64*)get_next_object(ht,key);
  }
  return(FALSE);
}

// get a column id - based on the cell/sample name
uint_64 get_col_id(const uint_64 cell, const uint_64 sample,BLABELS *map) {
  return(blabel2id(cell,sample,map));
}


hashtable load_whitelist(const char* file,uint hashsize) {

  FILE *fd;
  if ((fd=fopen(file,"r"))==NULL) {
    PRINT_ERROR("Failed to open file %s", file);  
    exit(1);
  }
  fprintf(stderr,"Loading whitelist from %s\n",file);
  // known barcodes
  hashtable ht=new_hashtable(hashsize);
  char buf[200];
  unsigned long num_read_b=0;
  while (!feof(fd) ) {
    char *l=fgets(&buf[0],200,fd);
    if (l==NULL || l[0]=='\0') continue;
    ++num_read_b;
    uint_64 num=char2uint_64(l);
    //fprintf(stderr,"%s->%llu\n",umi,umi_num);
    ulong key;
    //if ( num> 4294967296) {
    //  key=num/1000000000;// convert to ulong
    //} else {
      key=num;
      //}
    uint_64 *ptr=(uint_64*)get_object(ht,key);
    while ( ptr!=NULL ) {
      if ( *ptr==num) break;
      ptr=(uint_64*)get_next_object(ht,key);
    }
    if ( ptr==NULL ) {
      uint_64 *e=(uint_64*)malloc(sizeof(uint_64));
      if ( e==NULL ) {
	  PRINT_ERROR("Failed to allocate memory\n");
	  exit(1);
      }
      *e=num;
      insere(ht,key,e);
    }
  }
  fclose(fd);
  fprintf(stderr,"Loading whitelist from %s...done.\n",file);
  
  hashtable_stats(ht);

  return(ht);
}

// Matrix Market format
// Header: rows columns entries
#define MM_SEP " "
void write2MM(const char* file, LABELS *rows_map, BLABELS *cols_map,  hashtable uniq_ht, UNIQ_KEYS* ukey, uint min_num_reads) {

  FILE *fd=NULL;
  if ((fd=fopen(file,"w+"))==NULL) {
    PRINT_ERROR("Failed to open file %s", file);  
    exit(1);
  }
  //
  fprintf(stderr,"Saving MM file %s...\n",file);
  write_map2fileL(file,"rows",rows_map);
  write_map2fileB(file,"cols",cols_map);

  // avoid gziping directly for now
  fprintf(fd,"%%%%MatrixMarket matrix coordinate real general\n");
  fprintf(fd,"%lu %lu ",rows_map->ctr,cols_map->ctr);
  long loc=ftell(fd);
  fprintf(fd,"%-10lu\n",0L);

  uint_64 ctr=0;
  while ( ukey!=NULL ) {
    ulong n=count_uniq_entries(uniq_ht,ukey,min_num_reads);
    if ( n>=min_num_reads) {
      fprintf(fd,"%lu%s%lu%s%lu\n",ukey->feat_id,MM_SEP,get_col_id(ukey->cell, ukey->sample,cols_map),MM_SEP,n);
    }
    ukey=ukey->next;
    ++ctr;
  }

  if ( ctr >= 9999999999 ) {
    fprintf(stderr,"ERROR: integer overflow (please contact the author of the program).\n");
    exit(1);
  }
  
  if ( ctr == 0 ) {
    fprintf(stderr,"ERROR: 0 quantified features.\n");   
    exit(1);
  }
  // finish header
  fseek(fd,loc,SEEK_SET);
  fprintf(fd,"%-10lu",ctr);
  fclose(fd);
  
  fprintf(stderr,"Saving MM file...done.\n");
}


void print_usage(int exit_status) {
    PRINT_ERROR("Usage: bam_umi_count --bam in.bam --ucounts output_filename [--min_reads 0] [--uniq_mapped|--multi_mapped]  [--dump filename] [--tag GX|TX] [--known_umi file_one_umi_per_line] [--ucounts_MM |--ucounts_tsv] [--ucounts_MM|--ucounts_tsv] [--ignore_sample]");
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
  unsigned long long num_umis_discarded=0;
  unsigned long long num_cells_discarded=0;
  UNIQ_KEYS_HT* key_list=NULL;
  LABELS* features_map=NULL;
  BLABELS* cols_map=NULL;

  char *bam_file=NULL;
  char *ucounts_file=NULL;
  char *dump_file=NULL;
  char *known_umi_file=NULL;
  char *known_cells_file=NULL;
  static int  mm_format=FALSE; // tsv by default
  static int verbose=0;  
  static int help=FALSE;
  static int ignore_sample=FALSE;
  static struct option long_options[] = {
    {"verbose", no_argument,       &verbose, TRUE},
    {"multi_mapped", no_argument,      &uniq_mapped_only, FALSE},
    {"uniq_mapped", no_argument,       &uniq_mapped_only, TRUE},
    {"ignore_sample", no_argument,       &ignore_sample, TRUE},
    {"ucounts_MM",  no_argument, &mm_format, TRUE},
    {"ucounts_tsv",  no_argument, &mm_format, FALSE},
    {"help",   no_argument, &help, TRUE},
    {"bam",  required_argument, 0, 'b'},
    {"known_umi",  required_argument, 0, 'k'},
    {"known_cells",  required_argument, 0, 'c'},
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
    
    int c = getopt_long (argc, argv, "b:u:d:t:x:c:",
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
    case 'c':
      known_cells_file=optarg;
      break;
    case 'x':
      strncpy(feat_tag,optarg,3);
      break;
    case 't':
      min_num_reads=atol(optarg);
      break;
    default:
      //print_usage(1);
      break;
    }
  }
  if (help ) {
    print_usage(0);
  }

  if ( bam_file == NULL ) print_usage(1);
  if ( ucounts_file == NULL ) print_usage(1);
  
  
  // fprintf(stderr,"HASHSIZE=%u\n",HASHSIZE);
  hashtable ht=new_hashtable(HASHSIZE);
  // uniq=feature|cell|sample 
  hashtable uniq_ht=new_hashtable(HASHSIZE);
  hashtable kumi_ht=NULL;   // UMIs white list
  hashtable kcells_ht=NULL; // cells white list

  
  // known UMIs
  if ( known_umi_file!=NULL ) {
    kumi_ht=load_whitelist(known_umi_file,1000001);
    fprintf(stderr,"UMIs whitelist %lu\n",kumi_ht->n_entries);
  }
  // known cells
  if ( known_cells_file!=NULL ) {
    kcells_ht=load_whitelist(known_cells_file,500001);
    fprintf(stderr,"Cells whitelist %lu\n",kcells_ht->n_entries);
  }
  
  // Open file and exit if error
  in = strcmp(bam_file, "-")? bam_open(bam_file, "rb") : bam_dopen(fileno(stdin), "rb"); 
  
  if (in == 0 ) {  
    PRINT_ERROR("Failed to open BAM file %s", bam_file);  
    return(PARAMS_ERROR_EXIT_STATUS);  
  }  

  FILE *ucounts_fd=NULL,*dump_fd=NULL;

  if ( dump_file!=NULL) 
    if ((dump_fd=fopen(dump_file,"w+"))==NULL) {
      PRINT_ERROR("Failed to open file %s", dump_file);  
      exit(1);
    }

  fprintf(stderr,"@min_num_reads=%u\n",min_num_reads);
  fprintf(stderr,"@uniq mapped reads=%u\n",uniq_mapped_only);
  fprintf(stderr,"@tag=%s\n",feat_tag);
  fprintf(stderr,"@unique counts file=%s\n",ucounts_file);

  features_map=init_labels(1000001);
  cols_map=init_blabels(100001);

  //
  // 
  bam1_t *aln=bam_init1();

  fprintf(stderr,"Processing %s\n",bam_file);

  // read header
  bam_header_read(in);

  // tmp buffers
  char buf[500];
  char buf2[500];
  uint8_t *nh;
  char *feat,*umi,*cell,*sample;
  // TODO: change alns to entries
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
    //if (aln->core.flag & BAM_FSECONDARY) continue; // use only primary alignments
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
      if ( !ignore_sample)
	sample=get_tag(aln,"BC");
      else
	sample=NULL;
#ifdef DEBUG      
      fprintf(stderr,"umi2-->%s %s %s %s",umi,cell,sample,feat);
#endif
      // convert the different barcodes to uint_64 ASAP
      uint_64 umi_i=char2uint_64(umi);
      
      // skip if the UMI is not valid
      if (kumi_ht!=NULL)
	if ( !valid_barcode(kumi_ht,umi_i) ) {
	  num_umis_discarded++;
	  continue;
	}

      uint_64 cell_i=char2uint_64(cell);
      if ( kcells_ht!=NULL) 
	if ( !valid_barcode(kcells_ht,cell_i) ) {
	  num_cells_discarded++;
	  continue;
	}

      uint_64 sample_i=char2uint_64(sample);
	    
      // 
      get_col_id(cell_i,sample_i, cols_map);

      //
      // feature id
      char *f=strtok(feat,",");
      char *prev_f=NULL;
      int n_feat=0;
      // number of features
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
	  int len1=strlen(f);
	  assert( len1+1 < FEAT_ID_MAX_LEN );

	  uint feat_id=label_str2id(f,features_map);
	    
	  COUNT_ENTRY* lookup=get_count_key(feat_id,umi_i,cell_i,sample_i,ht);
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
	  // 
	  uniq_ht=add_count_entry_to_uht(uniq_ht,lookup);
	  // new entry
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
  fprintf(stderr,"%lld UMIs discarded\n",num_umis_discarded);
  fprintf(stderr,"%lld cells discarded\n",num_cells_discarded);
  fprintf(stderr,"%lld features\n",label_entries(features_map));
  fprintf(stderr,"%lld cells/samples\n",blabel_entries(cols_map));
  if ( !num_tags_found ) {
    fprintf(stderr,"ERROR: no valid alignments tagged with %s were found in %s.\n",feat_tag,bam_file);
    exit(1);
  }
  hashtable_stats(uniq_ht);
  //hashtable_stats(key_list->ht);
  // uniq UMIs that overlap each gene per cell (and optionally per sample)
  // traverse the hashtable and count how many distinct UMIs are assigned to each gene (with the number of reads above minimum number of reads threshold)
  if ( ucounts_file !=NULL) {
    UNIQ_KEYS *ukey=NULL;
    if ( key_list!=NULL) ukey=key_list->keys;
    if (mm_format ) {
      write2MM(ucounts_file,features_map,cols_map,uniq_ht,ukey,min_num_reads);
    } else {
      // TSV
      if ((ucounts_fd=fopen(ucounts_file,"w+"))==NULL) {
	PRINT_ERROR("Failed to open file %s", ucounts_file);  
	exit(1);
      }      
      uint_64 ctr=0;
      short pheader=TRUE;  
      while ( ukey!=NULL ) {
	ulong n=count_uniq_entries(uniq_ht,ukey,min_num_reads);
	if ( n>=min_num_reads) {
	  print_ucount(ukey,n,"\t",ucounts_fd,features_map,pheader);
	  pheader=FALSE;
	}
	ukey=ukey->next;
	++ctr;
      }
      if ( ctr == 0 ) {
	fprintf(stderr,"ERROR: 0 quantified features.");   
	exit(1);
      }
      fclose(ucounts_fd);
    }
  }
  // dump the counts
  if ( dump_file != NULL ) {
    init_hash_traversal(ht);
    COUNT_ENTRY* e;
    short pheader=TRUE;  
    while((e=(COUNT_ENTRY*)next_hash_object(ht))!=NULL) {
      print_count_entry(e,"\t",dump_fd,features_map,pheader,min_num_reads);
      pheader=FALSE;
    }
    fclose(dump_fd);
  }
  return(0);
}
