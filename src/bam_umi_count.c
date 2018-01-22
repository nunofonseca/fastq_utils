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


//#########################################
#define uint_64 unsigned long long
#include "hash.h"
#include "fastq.h"
#include "sam_tags.h"

#define FEAT_ID_MAX_LEN 25
#define MAX_BARCODE_LEN 19

// default values
#define MAX_CELLS    1000000
#define MAX_FEATURES 100000
#define MAX_SAMPLES  1


// 
typedef struct count_ENTRY {
  uint_64 umi;
  uint sample;
  float umi_obs;
  float reads_obs;
  struct count_ENTRY* next;
} COUNT_ENTRY;

typedef struct feature_ENTRY {
  uint feat_id;
  float tot_umi_obs;
  float tot_reads_obs;
  COUNT_ENTRY** samples; //the array should be close  to the number of expected samples
} FEATURE_ENTRY;

typedef struct cell_ENTRY {
 // uint_64 cell;
  hashtable features;
  float tot_umi_obs;
  float tot_reads_obs;
} CELL;

typedef struct db {
  CELL* cells; // only one cell at a time
  uint max_features;
  uint max_cells;
  uint features_cell;
  float tot_umi_obs;
  float tot_reads_obs;
} DB;



// ---------------------------------------------
// single label => feature
// uint max 4,294,967,295
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

// barcodes as uint_64
typedef struct blabel2id {
  uint id;
  uint_64 label;
  struct blabel2id* next;
} BLABEL2ID;

typedef struct blabels {
  BLABEL2ID *first;
  BLABEL2ID *last;
  uint  ctr;
  hashtable ht;
} BLABELS;

// ------------------------------------------------


static ulong hash_str(const char *str);


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
uint blabel_entries(const BLABELS *lm) {
  return(lm->ctr);
}
// return the id or 0 if label is not in the mapping
uint_64 label2id(const uint_64 lab,BLABELS* lm ) {
  assert(lm!=NULL);
  // lookup
  ulong ikey=lab;
  BLABEL2ID *e=(BLABEL2ID*)get_object(lm->ht,ikey);
  
  // look for the match
  while (e!=NULL) {
    if ( e->label==lab ) break;
    e=e->next;
  }
  if (e==NULL) return 0;
  return(e->id);  
}
// barcode based label to id
uint_64 blabel2id(const uint_64 lab,BLABELS* lm ) {
  //
  assert(lm!=NULL);
  // lookup
  ulong ikey=lab;
  BLABEL2ID *e=(BLABEL2ID*)get_object(lm->ht,ikey);
  
  // look for the match
  while (e!=NULL) {
    if ( e->label==lab ) break;
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
  uint_642char(e->label,&buf[0]);	      
  return(&buf[0]);
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
      fprintf(fd,"%u%s%s\n",(uint)e->id,MAPSEP,e->label);
      e=e->next;
      ++ctr;
    }
  }
  fclose(fd);
}

void write_map2fileB(const char* file,const char *labname,BLABELS* map,char* suffix) {
  FILE *fd;
  char buf[300];
  sprintf(&buf[0],"%s_%s",file,labname);
  file=buf;
  char empty[1];
  if (suffix==NULL) {
    empty[0]='\0';
    suffix=&empty[0];
  }
  
  if ((fd=fopen(buf,"w+"))==NULL) {
    PRINT_ERROR("Failed to open file %s for writing", buf);  
    exit(1);
  }
  if ( map->first!=NULL) {
    uint_64 ctr=0;
    BLABEL2ID *e=map->first;    
    while ( e!=NULL ) {
      fprintf(fd,"%u%s%s%s\n",e->id,MAPSEP,blabel_id2str(e),suffix);
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
DB* new_db(uint max_cells,uint max_features,uint features_cell) {

  DB *new=(DB*)malloc(sizeof(DB));
  if (new==NULL) { return(NULL);}
  new->max_cells=max_cells+1;// index is from 1-max
  new->max_features=max_features;
  new->features_cell=features_cell;
  new->tot_umi_obs=new->tot_reads_obs=0;
  new->cells=(CELL*)malloc(new->max_cells*sizeof(CELL));
  if (new->cells==NULL) { return(NULL);}
  memset(new->cells,0L,sizeof(CELL)*new->max_cells);
  if (max_cells==0) new->max_cells=0; // sorted bam
  return(new);
}

DB* quick_reset_db(DB *db) {

  int x=0;
  db->cells[0].tot_umi_obs=0;
  db->cells[0].tot_reads_obs=0;
  if (db->cells[0].features!=NULL ) {
    // travesse the features
    init_hash_traversal(db->cells[0].features);
    FEATURE_ENTRY* e;
    while ((e=(FEATURE_ENTRY*)next_hash_object(db->cells[0].features))!=NULL) {
      e->tot_umi_obs=0;
      e->tot_reads_obs=0;
      for (x=0;x<=MAX_SAMPLES;++x) {
	if (e->samples[x]!=NULL) {
	  COUNT_ENTRY *todel=e->samples[x];	    
	  e->samples[x]=e->samples[x]->next;
	  free(todel);
	}
      }
    }    
  }
  return(db);
}

/* // check if an entry exists, returns a pointer for the the existing COUNT_ENTRY if available or creates a new one and increments the read/umis counters with the incr value */
COUNT_ENTRY* get_count_entry(const uint feat_id,const uint_64 umii,const uint cell_id,const uint sample_id,DB* db,float incr) {

  uint cell_idx=cell_id;
  if( cell_id>db->max_cells && db->max_cells>0 ) {
    PRINT_ERROR("Too many cells %u - please rerun and increase the cells using the --max_cells parameter\n",cell_id);
    exit(1);
  }
  // lookup
  if ( db->max_cells==0) cell_idx=0;
  if ( db->cells[cell_idx].features==NULL ) {
    // initialize entry
    //db->cells[cell_idx]->features=new_hashtable(db->max_features/2);
    // only a fraction of the genes are expressed
    db->cells[cell_idx].features=new_hashtable(db->features_cell);
    db->cells[cell_idx].tot_umi_obs=0;
    db->cells[cell_idx].tot_reads_obs=0;
  }

  
  // lookup for the feature
  ulong key=feat_id;
  FEATURE_ENTRY* e=(FEATURE_ENTRY*)get_object(db->cells[cell_idx].features,key);
  while (e!=NULL) {
    if ( feat_id == e->feat_id ) break;
    e=(FEATURE_ENTRY*)get_next_object(db->cells[cell_idx].features,key);
  }
  
  if ( e==NULL) {
    // new feature entry
    e=(FEATURE_ENTRY*)malloc(sizeof(FEATURE_ENTRY));
    if (e==NULL) return(NULL);
    e->tot_umi_obs=e->tot_reads_obs=0;
    e->feat_id=feat_id;
    if (insere(db->cells[cell_idx].features,feat_id,e)<0) {
      return(NULL);
    }    
    e->samples=(COUNT_ENTRY**)malloc(sizeof(COUNT_ENTRY*)*MAX_SAMPLES);
    if (e->samples==NULL) return(NULL);
    memset(e->samples,0L,sizeof(COUNT_ENTRY*)*MAX_SAMPLES);    
  }
  ulong sample_key=sample_id%MAX_SAMPLES;
  COUNT_ENTRY *ce=e->samples[sample_key];
  while(ce!=NULL) {
    if ( ce->sample==sample_id && ce->umi==umii ) break;
    ce=ce->next;
  }
  if (ce==NULL ) {
    COUNT_ENTRY *new_ce=(COUNT_ENTRY*)malloc(sizeof(COUNT_ENTRY));
    if (new_ce==NULL) return(NULL);
    new_ce->umi=umii;
    new_ce->sample=sample_id;
    new_ce->umi_obs=incr;
    new_ce->reads_obs=incr;
    new_ce->next=e->samples[sample_key];
    e->samples[sample_key]=new_ce;       
    e->tot_umi_obs+=incr;
    db->tot_umi_obs+=incr;

    ce=new_ce;
  }
  e->tot_reads_obs+=incr;
  db->tot_reads_obs+=incr;

  return(ce);
}

char EMPTY_STRING[]="";
char *get_tag(bam1_t *aln,const char tagname[2]) {

  uint8_t *s=bam_aux_get(aln,tagname);
  if (s==0) return(EMPTY_STRING);

  char *s2=bam_aux2Z(s);
  if (s2==0) return(EMPTY_STRING);

  return(s2);
}

int valid_barcode(hashtable ht,const uint_64 barcode_id) {
  if (ht==NULL) return(TRUE); // by default all BARCODEs are valid

  ulong key=barcode_id;

  uint_64 *ptr=(uint_64*)get_object(ht,key);
  while ( ptr!=NULL ) {
    if ( *ptr==barcode_id) return(TRUE);
    ptr=(uint_64*)get_next_object(ht,key);
  }
  return(FALSE);
}

// get a column id - based on the cell/sample name
//uint_64 get_col_id(const uint_64 cell, const uint_64 sample,BLABELS *map) {
//  return(blabel2id(cell,sample,map));
//}


hashtable load_whitelist(const char* file,uint hashsize,BLABELS* map) {

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
    ulong key=num;
    if ( map!=NULL ) {
      //add to the map
      num=key=blabel2id(num,map);
    }
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
      insere(ht,num,e);
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
void write2MM(const char* file, DB*db,LABELS *rows_map, BLABELS *cols_map,uint min_num_reads,uint min_num_umis,char *cell_suffix,int UMI) { 

  FILE *fd=NULL;
  if ((fd=fopen(file,"w+"))==NULL) {
    PRINT_ERROR("Failed to open file %s", file);
    exit(1);
  }
  //
  fprintf(stderr,"Saving MM file %s...\n",file);
  write_map2fileL(file,"rows",rows_map);
  write_map2fileB(file,"cols",cols_map,cell_suffix);

  // avoid gziping directly for now
  fprintf(fd,"%%%%MatrixMarket matrix coordinate real general\n");
  fprintf(fd,"%u %u ",rows_map->ctr,cols_map->ctr);
  long loc=ftell(fd);
  fprintf(fd,"%-15lu\n",0L);

  // traverse the full DB
  // todo: handle samples
  uint cell_id=0;
  uint_64 tot_ctr=0;
  uint_64 tot_cells=0;
  uint_64 tot_feat_cells=0;
  

  FEATURE_ENTRY *fe;
  while ( cell_id<=db->max_cells) {
    // traverse the feature hashtable
    if (db->cells[cell_id].features!=NULL) {
      init_hash_traversal(db->cells[cell_id].features);
      ++tot_cells;
      while ( (fe=(FEATURE_ENTRY*)next_hash_object(db->cells[cell_id].features))!=NULL) {
	//fprintf(stderr,">>%u%s%u%s%f %f %u %u\n",fe->feat_id,MM_SEP,cell_id,MM_SEP,fe->tot_umi_obs,fe->tot_reads_obs,min_num_reads,min_num_umis);
	// do not go down through the samples
	if ( fe->tot_reads_obs>=min_num_reads*1.0 &&
	     fe->tot_umi_obs>=min_num_umis*1.0 ) {	  
	  if ( UMI==TRUE && (uint)fe->tot_umi_obs>=1 ) {
	    fprintf(fd,"%u%s%u%s%u\n",fe->feat_id,MM_SEP,cell_id,MM_SEP,(uint)round(fe->tot_umi_obs));
	    tot_ctr+=(uint)fe->tot_umi_obs;
	    ++tot_feat_cells;
	  } else if ( (uint)fe->tot_reads_obs>=1)  {
	    fprintf(fd,"%u%s%u%s%u\n",fe->feat_id,MM_SEP,cell_id,MM_SEP,(uint)round(fe->tot_reads_obs));
	    tot_ctr+=(uint)fe->tot_reads_obs;
	    ++tot_feat_cells;
	  }
	}
      }
    }
    ++cell_id;    
  }  
  if ( tot_feat_cells >= 9999999999 ) {
    fprintf(stderr,"ERROR: integer overflow (please contact the author of the program).\n");
    exit(1);
  }
  
  if ( tot_feat_cells == 0 ) {
    fprintf(stderr,"ERROR: 0 quantified features.\n");
    exit(1);
  }
  // finish header
  fseek(fd,loc,SEEK_SET);
  fprintf(fd,"%-15llu",tot_feat_cells);
  fclose(fd);
  
  fprintf(stderr,"Saving MM file...done.\n");
  fprintf(stderr,"#cells/features: %llu\n",tot_feat_cells);
  fprintf(stderr,"#cells: %llu\n",tot_cells);
  fprintf(stderr,"#tot expr: %llu\n",tot_ctr);
}


void cell2MM(DB*db, FILE *counts_fd,int UMI,uint min_num_reads,uint min_num_umis,uint_64* tot_ctr,uint_64* tot_feat_cells,const uint cell_id) {

  uint cell_idx=cell_id;
  if ( db->max_cells == 0 ) cell_idx=0;
  
  if (db->cells[cell_idx].features!=NULL) {
    FEATURE_ENTRY *fe;
    init_hash_traversal(db->cells[cell_idx].features);
    while ( (fe=(FEATURE_ENTRY*)next_hash_object(db->cells[cell_idx].features))!=NULL) {
      // do not go down through the samples
      if ( fe->tot_reads_obs>=min_num_reads*1.0 &&
	   fe->tot_umi_obs>=min_num_umis*1.0 ) {	  
	if ( UMI==TRUE && (uint)fe->tot_umi_obs>=1 ) {
	  fprintf(counts_fd,"%u%s%u%s%u\n",fe->feat_id,MM_SEP,cell_id,MM_SEP,(uint)round(fe->tot_umi_obs));
	  *tot_ctr+=(uint)fe->tot_umi_obs;
	  ++*tot_feat_cells;
	} else if ( (uint)fe->tot_reads_obs>=1)  {
	  fprintf(counts_fd,"%u%s%u%s%u\n",fe->feat_id,MM_SEP,cell_id,MM_SEP,(uint)round(fe->tot_reads_obs));
	  *tot_ctr+=(uint)fe->tot_reads_obs;
	  ++*tot_feat_cells;
	}
      }
    }
  }
}

FILE* MM_header(const char* counts_file,long *header_loc) {
  FILE *fd;
  if ((fd=fopen(counts_file,"w+"))==NULL) {
    PRINT_ERROR("Failed to open file %s", counts_file);
    exit(1);
  }
  //
  fprintf(stderr,"Creating MM file %s...\n",counts_file);
  // save the header with estimates of the matrix size
  // avoid gziping directly for now
  fprintf(fd,"%%%%MatrixMarket matrix coordinate real general\n");
  *header_loc=ftell(fd);
  fprintf(fd,"%-10lu %-10lu %-15llu\n",0L,0L,0LL);
  return(fd);
}

void print_usage(int exit_status) {
    PRINT_ERROR("Usage: bam_umi_count --bam in.bam --ucounts output_filename [--min_reads 0] [--min_umis 0] [--uniq_mapped|--multi_mapped]  [--dump filename] [--tag gx|tx] [--known_umi file_one_umi_per_line] [--ucounts_MM |--ucounts_tsv] [--ucounts_MM|--ucounts_tsv] [--ignore_sample] [--cell_suffix suffix] [--max_cells number] [--max_feat number] [--feat_cell number] [--cell_tag tag] [--sorted_by_cell]");
    if ( exit_status>=0) exit(exit_status);
}


int main(int argc, char *argv[])  
{  
  bamFile in; 
  uint min_num_reads=0;
  uint min_num_umis=0;

  char feat_tag[]=GENE_ID_TAG;
  char cell_tag[]=CELL_TAG;
  unsigned long long num_alns=0;
  unsigned long long num_tags_found=0;
  unsigned long long num_umis_discarded=0;
  unsigned long long num_cells_discarded=0;

  // mappings
  LABELS* features_map=NULL;
  BLABELS* cells_map=NULL;
  BLABELS* samples_map=NULL;
  // TODO: allow these values to be passed as arguments
  ulong max_features=MAX_FEATURES;
  ulong max_cells=MAX_CELLS;
  ulong max_samples=MAX_SAMPLES;
  ulong features_cell=4000;
  ulong ncells=0;
  
  char *bam_file=NULL;
  char *ucounts_file=NULL;
  char *rcounts_file=NULL;

  char *known_umi_file=NULL;
  char *known_cells_file=NULL;
  char *cell_suffix=NULL;

  static int bam_sorted_by_cell=FALSE; // reduced memory usage
  static int uniq_mapped_only=FALSE;
  static int verbose=0;  
  static int help=FALSE;
  static int ignore_sample=FALSE;
  static struct option long_options[] = {
    {"verbose", no_argument,       &verbose, TRUE},
    {"multi_mapped", no_argument,      &uniq_mapped_only, FALSE},
    {"uniq_mapped", no_argument,       &uniq_mapped_only, TRUE},
    {"sorted_by_cell", no_argument,       &bam_sorted_by_cell, TRUE},
    {"ignore_sample", no_argument,       &ignore_sample, TRUE},
    {"help",   no_argument, &help, TRUE},
    {"bam",  required_argument, 0, 'b'},
    {"cell_suffix",  required_argument, 0, 's'},
    {"known_umi",  required_argument, 0, 'k'},
    {"known_cells",  required_argument, 0, 'c'},
    {"ucounts",  required_argument, 0, 'u'},
    {"rcounts",  required_argument, 0, 'r'},
    {"tag",  required_argument, 0, 'x'},
    {"cell_tag",  required_argument, 0, 'X'},
    {"min_reads",  required_argument, 0, 't'},
    {"min_umis",  required_argument, 0, 'U'},
    {"max_cells",  required_argument, 0, 'C'},
    {"max_feat",  required_argument, 0, 'F'},
    {"feat_cell",  required_argument, 0, 'T'},
    {0,0,0,0}
  };

  fprintf(stderr,"bam_umi_count version %s\n",VERSION);
  // process arguments
  while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    int c = getopt_long (argc, argv, "F:T:C:b:U:u:r:t:x:c:s:h",
		     long_options, &option_index);      
    if (c == -1) // no more options
      break;

    switch (c) {
    case 'h':
      help=TRUE;
      break;
    case 'b':
      bam_file=optarg;
      break;
    case 'u':
      ucounts_file=optarg;
      break;
    case 'r':
      rcounts_file=optarg;
      break;
    case 's':
      cell_suffix=optarg;
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
    case 'X':
      strncpy(cell_tag,optarg,3);
      break;
    case 't':
      min_num_reads=atol(optarg);
      break;
    case 'U':
      min_num_umis=atol(optarg);
      break;
    case 'C':
      max_cells=atol(optarg);
      break;
    case 'F':
      max_features=atol(optarg);
      break;
    case 'T':
      features_cell=atol(optarg);
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
  

    // mappings
  features_map=init_labels(max_features);
  cells_map=init_blabels(max_cells);
  if (!ignore_sample)
    samples_map=init_blabels(max_samples);

  DB *db=NULL;
  if ( bam_sorted_by_cell ) {
    db=new_db(0,max_features,features_cell);
  } else {
    db=new_db(max_cells,max_features,features_cell);
  }

  // white lists
  hashtable kumi_ht=NULL;   // UMIs white list
  hashtable kcells_ht=NULL; // cells white list
  
  // known UMIs
  if ( known_umi_file!=NULL ) {
    kumi_ht=load_whitelist(known_umi_file,1000001,NULL);
    fprintf(stderr,"UMIs whitelist %llu\n",kumi_ht->n_entries);
  }
    
  // known cells
  if ( known_cells_file!=NULL ) {
    kcells_ht=load_whitelist(known_cells_file,max_cells/2,cells_map);
    fprintf(stderr,"Cells whitelist %llu\n",kcells_ht->n_entries);
  }

  // Open file and exit if error
  in = strcmp(bam_file, "-")? bam_open(bam_file, "rb") : bam_dopen(fileno(stdin), "rb"); 
  
  if (in == 0 ) {  
    PRINT_ERROR("Failed to open BAM file %s", bam_file);  
    return(PARAMS_ERROR_EXIT_STATUS);  
  }  

  fprintf(stderr,"@min_num_reads=%u\n",min_num_reads);
  fprintf(stderr,"@min_num_umis=%u\n",min_num_umis);
  fprintf(stderr,"@uniq mapped reads=%u\n",uniq_mapped_only);
  fprintf(stderr,"@sorted bam=%u\n",bam_sorted_by_cell);
  fprintf(stderr,"@tag=%s\n",feat_tag);
  fprintf(stderr,"@unique counts file=%s\n",ucounts_file);
  if (cell_suffix!=NULL)
    fprintf(stderr,"@cell_suffix=%s\n",cell_suffix);

  long header_loc=0;
  long rheader_loc=0;
  FILE *counts_fd=NULL;
  FILE *rcounts_fd=NULL;
  //
  // 
  bam1_t *aln=bam_init1();
  // read header
  bam_header_read(in);
  
  fprintf(stderr,"Processing %s\n",bam_file);

  if ( bam_sorted_by_cell ) {
    if ( ucounts_file !=NULL) { 
      counts_fd=MM_header(ucounts_file,&header_loc);
    }
    if ( rcounts_file !=NULL) { 
      rcounts_fd=MM_header(rcounts_file,&rheader_loc);
    }
  }

  // tmp buffers
  uint8_t *nh;
  char *feat,*umi,*cell,*sample;
  // map to ids
  uint cell_id=0;
  uint prev_cell_id=0;
  uint sample_id=0;

  // 
  uint_64 tot_ctr=0;
  uint_64 tot_feat_cells=0;

  // TODO: change alns to entries
  num_alns=0;
  if ( bam_sorted_by_cell ) fprintf(stderr,"Cells processed\n");
  while(bam_read1(in,aln)>=0) { // read alignment
    if ( num_alns == ULLONG_MAX ) {
      PRINT_ERROR("counter overflow (number of alignments) - %llu\n",num_alns);
      exit(3);
    }
    ++num_alns;
    if ( ! bam_sorted_by_cell && num_alns%100000==0) { fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%llu",num_alns); fflush(stderr); }

    if (aln->core.tid < 0) continue;//ignore unaligned reads
    if (aln->core.flag & BAM_FUNMAP) continue;
    //if (aln->core.flag & BAM_FSECONDARY) continue; // use only primary alignments
    if (aln->core.flag & BAM_FPAIRED & BAM_FPROPER_PAIR & BAM_FREAD2 ) continue; // avoid double counting

    //assert(r!=NULL);
    nh = bam_aux_get(aln, "NH");
    int nh_i=1;
    // exclude multimaps?
    if (nh!=NULL) {
      nh_i=bam_aux2i(nh);
      if ( nh_i > 1 && uniq_mapped_only) continue;
    }
    //int is_paired=aln->core.flag & BAM_FPAIRED;
    feat=get_tag(aln,feat_tag);
    if (feat[0]!='\0') {
      num_tags_found++;
      umi=get_tag(aln,UMI_TAG);
      // TODO: remove support for this tag and update tests
      if ( umi[0]=='\0') // backwards compatibility
	umi=get_tag(aln,"UM");
      cell=get_tag(aln,cell_tag);
      if ( !ignore_sample)
	sample=get_tag(aln,SAMPLE_TAG);
      else
	sample=NULL;
#ifdef DEBUG      
      fprintf(stderr,"umi2-->%s %s %s %s",umi,cell,sample,feat);
#endif
      // convert the different barcodes to uint_64 
      uint_64 umi_i=char2uint_64(umi);
      
      // skip if the UMI is not valid
      if (kumi_ht!=NULL)
	if ( !valid_barcode(kumi_ht,umi_i) ) {
	  num_umis_discarded++;
	  continue;
	}
      
      uint_64 cell_i=char2uint_64(cell);

      //fprintf(stderr,"aaaa1:%s-->%llu-->%lu\n",cell,cell_i,cell_id);      
      if ( kcells_ht!=NULL) {
	cell_id=label2id(cell_i,cells_map);
	if ( !valid_barcode(kcells_ht,cell_id) ) {
	  num_cells_discarded++;
	  continue;
	}
      } else {
	cell_id=blabel2id(cell_i,cells_map);
      }
      if ( bam_sorted_by_cell ) {
	if ( prev_cell_id != cell_id ) {
	  if ( cell_id <= prev_cell_id ) {
	    fprintf(stderr,"Error: The BAM file does not seem to be sorted by CR\n");
	    exit(1);
	  }
	  
	  if ( prev_cell_id!=0 ) {
	    ++ncells;
	    fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b%-10llu",ncells);
	    cell2MM(db,counts_fd,TRUE,min_num_reads,min_num_umis,&tot_ctr,&tot_feat_cells,prev_cell_id);
	    if ( rcounts_fd!=NULL )
	      cell2MM(db,rcounts_fd,FALSE,min_num_reads,min_num_umis,&tot_ctr,&tot_feat_cells,prev_cell_id);
	    // init/reset data structures
	    db=quick_reset_db(db);
	  }
	}
	prev_cell_id=cell_id;
      }      
    
      uint_64 sample_i=char2uint_64(sample);
      // map to ids
      
      if ( ! ignore_sample )
	sample_id=blabel2id(sample_i,samples_map);
      
      //
      // feature id
      char *f=strtok(feat,",");
      char *prev_f=NULL;
      int n_feat=0;
      // number of features (count the number of , or ;
      while (f !=NULL ) {
	if (prev_f==NULL || !strcmp(f,prev_f)) 
	  ++n_feat;
	prev_f=f;
	f=strtok(NULL,",");
      }

      f=strtok(feat,",");
      prev_f=NULL;
      float incr=1.0/(n_feat*nh_i);
      while (f !=NULL ) {
	if (prev_f==NULL || !strcmp(f,prev_f)) {
	  int len1=strlen(f);
	  assert( len1+1 < FEAT_ID_MAX_LEN );
	  uint feat_id=label_str2id(f,features_map);

	  COUNT_ENTRY* lookup=get_count_entry(feat_id,umi_i,cell_id,sample_id,db,incr);
	  //fprintf(stderr,">>>>%u-->%f\n",cell_id,incr);
	  if (lookup==NULL) {
	    fprintf(stderr,"Memory allocation failed!\n");
	    exit(1);
	  }
	  // check for overflow
	  if ( lookup->reads_obs + 1 >= FLT_MAX ) {
	    PRINT_ERROR("Counter overflow (umi_obs) %.2lf\n",lookup->reads_obs);	    
	    exit(1);
	  }
	  if ( lookup->umi_obs + 1 >= FLT_MAX ) {
	    PRINT_ERROR("Counter overflow (umi_obs) %.2lf\n",lookup->umi_obs);	    
	    exit(1);
	  }
	}
	prev_f=f;
	f=strtok(NULL,",");
      }
    }
  }
  if ( bam_sorted_by_cell ) {
    // last cell
    if ( cell_id!=0 ) {
      ++ncells;
      fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b%-10llu",ncells);
      cell2MM(db,counts_fd,TRUE,min_num_reads,min_num_umis,&tot_ctr,&tot_feat_cells,cell_id);
      if ( rcounts_fd!=NULL ) 
	cell2MM(db,rcounts_fd,FALSE,min_num_reads,min_num_umis,&tot_ctr,&tot_feat_cells,cell_id);
    }
  }

  fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\n");fflush(stderr);
  bam_destroy1(aln);
  // write output
  fprintf(stderr,"Alignments processed: %llu\n",num_alns);
  fprintf(stderr,"%s encountered  %llu times\n",feat_tag,num_tags_found);
  fprintf(stderr,"%lld UMIs discarded\n",num_umis_discarded);
  fprintf(stderr,"%lld cells discarded\n",num_cells_discarded);
  fprintf(stderr,"%u features\n",label_entries(features_map));
  fprintf(stderr,"%u cells\n",blabel_entries(cells_map));

  if (samples_map!=NULL)
    fprintf(stderr,"%u samples\n",blabel_entries(samples_map));
  fprintf(stderr,"%f total reads\n",db->tot_reads_obs);
  fprintf(stderr,"%f total UMI\n",db->tot_umi_obs);
  if ( !num_tags_found ) {
    fprintf(stderr,"ERROR: no valid alignments tagged with %s were found in %s.\n",feat_tag,bam_file);
    exit(1);
  }

  if ( bam_sorted_by_cell ) {
    if (counts_fd!=NULL) {
      // finish header
      fseek(counts_fd,header_loc,SEEK_SET);
      // update the header
      fprintf(counts_fd,"%-10u %-10u %-15llu",features_map->ctr,cells_map->ctr,tot_feat_cells);
      // write the two aux files
      write_map2fileL(ucounts_file,"rows",features_map);
      write_map2fileB(ucounts_file,"cols",cells_map,cell_suffix);  
      fclose(counts_fd);
    }
    if (rcounts_fd!=NULL) {
      fseek(rcounts_fd,rheader_loc,SEEK_SET);
      // update the header
      fprintf(rcounts_fd,"%-10u %-10u %-15llu",features_map->ctr,cells_map->ctr,tot_feat_cells);
      // write the two aux files
      write_map2fileL(rcounts_file,"rows",features_map);
      write_map2fileB(rcounts_file,"cols",cells_map,cell_suffix);
    }
    exit(0);
  }

  // uniq UMIs that overlap each gene per cell (and optionally per sample)
  if ( ucounts_file !=NULL) { 
    // todo: use a cols_map to take the sample barcode into account
    write2MM(ucounts_file,db,features_map,cells_map,min_num_reads,min_num_umis,cell_suffix,TRUE); 
  }
  // dump the counts */
  if ( rcounts_file != NULL ) {
    write2MM(rcounts_file,db,features_map,cells_map,min_num_reads,min_num_umis,cell_suffix,FALSE); 
  } 

  return(0);
}
