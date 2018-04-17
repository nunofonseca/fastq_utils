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

#include "bam.h"

#include "hash.h"
#include "fastq.h"
#include "sam_tags.h" 

typedef enum  { READ1=1, READ2=2, READ3=3 , READ4=4, READ5=5 } READ_IDX;
#define INDEX1 3
#define INDEX2 4
#define INDEX3 5

typedef long FASTQ_SIZE;

struct params_s {
  char *file[INDEX3+1];
  char *outfile[2+1]; // index files are not touched
  short phred_encoding; //33/64 - map the ascii code to the 0 phred value
  FASTQ_BOOLEAN paired; // maximum read length
  FASTQ_BOOLEAN verbose;
  FASTQ_READ_OFFSET read_offset[2+1];//
  FASTQ_SIZE read_size[2+1]; //
  READ_IDX cell_read;  // cell barcode
  FASTQ_READ_OFFSET cell_offset;//
  FASTQ_SIZE cell_size;  //
  READ_IDX sample_read;  // sample barcode
  FASTQ_READ_OFFSET sample_offset;// sample barcode
  FASTQ_SIZE sample_size;  //
  READ_IDX umi_read;  // Read with the UMI (1/2)
  FASTQ_READ_OFFSET umi_offset; //UMI offset
  FASTQ_SIZE umi_size; //
  short interleaved[2]; //
  short has_interleaved_entries;
  short min_qual; // discard read if the base quality of the UMI/cell barc code is below the threshold
  short num_input_files;
};
typedef struct params_s Params;


#define PRINT_VERBOSE(p,s...) {if (p->verbose) fprintf(stderr,##s );  }

#define PRINT_SAM(s...) {fprintf(stdout,##s );}


READ_IDX read_index2read_idx(const char *s) {

  if (!strcmp(s,"read1")) return READ1;
  if (!strcmp(s,"read2")) return READ2;
  if (!strcmp(s,"index1")) return READ3;
  if (!strcmp(s,"index2")) return READ4;
  if (!strcmp(s,"index3")) return READ5;

  PRINT_ERROR("invalid file reference %s (valid values are read1,read2, index1,index2,index3)\n",s);
  exit(1);    
}

void validate_options(Params* options) {

  if ( options->file[1]==NULL ) {
    PRINT_ERROR("missing input file (-read1)");
    exit(1);
  }
  if ( options->paired && options->file[2]==NULL ) {
    PRINT_ERROR("if paired_end is used then two fastq files should be provided - missing input file (-read2)");
    exit(PARAMS_ERROR_EXIT_STATUS);
  }

  if ( options->outfile[1]==NULL ) {
    PRINT_ERROR("if single_end then -outfile1 should be provided");
    exit(PARAMS_ERROR_EXIT_STATUS);
  }
  
}


Params* Params_new(void) {
  Params* new=(Params*)malloc(sizeof(Params));
  if (new==NULL) {
    PRINT_ERROR("unable to allocate %ld bytes of memory",sizeof(Params));
    exit(SYS_INT_ERROR_EXIT_STATUS);
  }

  new->file[READ1]=NULL;
  new->file[READ2]=NULL;
  new->file[READ3]=NULL;
  new->file[READ4]=NULL;
  new->file[READ5]=NULL;
  new->outfile[READ1]=NULL;
  new->outfile[READ2]=NULL;
  new->phred_encoding=64;
  new->paired=FALSE; // maximum read length
  new->read_offset[READ1]=UNDEF;//
  new->read_size[READ1]=0; //
  new->read_offset[READ2]=UNDEF;//
  new->read_size[READ2]=0; //
  new->cell_read=UNDEF;  // cell barcode: default read1
  new->cell_offset=UNDEF;//
  new->cell_size=0;  //
  new->sample_read=UNDEF;  // sample barcode
  new->sample_offset=UNDEF;// sample barcode
  new->sample_size=0;  //
  new->umi_read=UNDEF;  // Read with the UMI (1/2)
  new->verbose=FALSE;
  new->umi_offset=UNDEF; //
  new->umi_size=0; //
  new->has_interleaved_entries=FALSE;
  new->interleaved[0]=0;
  new->interleaved[1]=0;
  new->min_qual=0; // discard read if the base quality of the UMI/cell barc code is below the threshold
  new->num_input_files=0;
  return(new);
}

void set_input_file(Params* p,char* filename,READ_IDX rdx) {
  if (rdx<READ1 || rdx>READ5 ) {
    PRINT_ERROR("internal error set_input_file");
    exit(SYS_INT_ERROR_EXIT_STATUS);
  }

  if (filename!=NULL && p->file[rdx]==NULL)
    p->num_input_files++;
  p->file[rdx]=filename;
}

void slice_read(FASTQ_ENTRY* m,const Params* p,READ_IDX cur_read) {

  //fprintf(stderr,">%d %ld %ld\n",cur_read,p->read_offset[cur_read],p->read_size[cur_read]);
  if (p->read_offset[cur_read]==UNDEF) return; // do nothing
  if (cur_read<READ3 && p->read_offset[cur_read]==0 && p->read_size[cur_read]==-1 ) return; // do not slice
  // truncate header2
  m->hdr2[1]='\n';
  m->hdr2[2]='\0';
  if (p->read_size[cur_read]==0) {
    // truncate - nothing to output
    m->seq[0]=m->qual[0]='\n';
    m->seq[1]=m->qual[1]='\0';
    return;
  }
  //
  FASTQ_READ_OFFSET offset=p->read_offset[cur_read];  
  if ( offset > 0 ) { // copy
    FASTQ_READ_OFFSET len=p->read_size[cur_read];
    FASTQ_READ_OFFSET x;
    if ( len==-1 ) len=m->read_len;
    for (x=0; x<=len;++x) {
      m->seq[x]=m->seq[x+offset];
      m->qual[x]=m->qual[x+offset];
    }
  }
  offset=p->read_size[cur_read];
  m->seq[offset]='\n';
  m->seq[offset+1]='\0';
  m->qual[offset]='\n';
  m->qual[offset+1]='\0';
}

void add_tags2readname(FASTQ_ENTRY* m,char*cell,char*umi,char*sample) {
  // add the tags right to the beginning of the read
  // STAGS_CELL=  _UMI= _SAMPLE= _ETAGS_ 31
  #define TAGS_NAMES_LEN 31
  FASTQ_READ_OFFSET  offset=TAGS_NAMES_LEN+strlen(cell)+strlen(umi)+strlen(sample);
  // go to the end
  unsigned int s=1;

  if (offset==TAGS_NAMES_LEN) return; // nothing to do

  //fprintf(stderr,"->%s\n",m->hdr1);
  while(m->hdr1[s]!='\0') ++s;
  // create a gap
  while (s>=1) {
    m->hdr1[s+offset]=m->hdr1[s];
    --s;
  }
  char c=m->hdr1[1+offset];
  sprintf(&m->hdr1[1],"STAGS_CELL=%s_UMI=%s_SAMPLE=%s_ETAGS_",cell,umi,sample);
  m->hdr1[1+offset]=c;
  // truncate the 2nd header
  m->hdr2[1]='\n';
  m->hdr2[2]='\0';
  //fprintf(stderr,"---->%s",m->hdr1);
}

short get_barcode(const FASTQ_ENTRY *m1,
		  short phred_encoding,
		  const READ_IDX  read,
		  const FASTQ_READ_OFFSET offset,
		  const FASTQ_SIZE size,
		  const short min_qual,		  
		  char* dest,
		  char* qual,
		  int get_qual) {
  dest[0]='\0';
   
  if ( read==UNDEF || offset==UNDEF || size==0 )
    return(TRUE); // Nothing TODO but all OK
  // check boundaries
  if ( offset >  m1->read_len-1 || offset+size > m1->read_len-1 ) {
    //
    fprintf(stderr,"Warning: Read too short - barcode not found\n");
    return(FALSE);
  }
  if ( min_qual > 0 ) {
    // check quality
    FASTQ_READ_OFFSET x;
    for (x=offset;x<offset+size;++x) {
      
      if (m1->qual[x]-phred_encoding < min_qual) {
	//fprintf(stderr,"Checking quality %c -> %d\n",m1->qual[x],m1->qual[x]-phred_encoding);
#ifdef DEBUG
	fprintf(stderr,"Warning: skipping read due to low quality in barcode sequence %d\n",m1->qual[x]-phred_encoding);
#endif	
	return(FALSE);	
      }
    }
  }
  //
  strncpy(dest,&m1->seq[offset],size);
  dest[size]='\0';
  if ( get_qual ) {
    strncpy(qual,&m1->qual[offset],size);
    qual[size]='\0';
  }
  return(TRUE);
}

// Extract all required information from a read based on the given configuration
int extract_info(const FASTQ_ENTRY *m1,const Params *p,
		 const READ_IDX read,// 1/2/3/4/5
		 char *cell,
		 char *umi,
		 char *sample,
		 char *cell_qual,
		 char *umi_qual,
		 char *sample_qual,
		 int keep_qual) {
  // UMI
  if ( p->umi_read==read) 
    if (!get_barcode(m1,p->phred_encoding,p->umi_read,p->umi_offset,p->umi_size,p->min_qual, umi,umi_qual,keep_qual))
      return(FALSE);    		    
  // sample
  if ( p->sample_read==read) 
    if (!get_barcode(m1,p->phred_encoding,p->sample_read,p->sample_offset,p->sample_size,p->min_qual, sample,sample_qual,keep_qual))
      return(FALSE);
  // cell
  if ( p->cell_read==read) 
    if (!get_barcode(m1,p->phred_encoding,p->cell_read,p->cell_offset,p->cell_size,p->min_qual, cell,cell_qual, keep_qual))
      return(FALSE);
    
  return(TRUE);
}

// return TRUE if at least one fd reached EOF
FASTQ_BOOLEAN fastq_files_eof(FASTQ_FILE** fdi,int n) {
  int k=1;
  while(k<=n) {
    if ( fdi[k]!=NULL )
      if (gzeof(fdi[k]->fd))
	return TRUE; // do not ensure that all files 
    ++k;
  }
  return FALSE;
}


char* format_read_name(char *s) {
  unsigned int i=0;
  while ( s[i]!='\0' ) {
    if ( s[i]==' ' ) s[i]='@';
    if ( s[i]=='\n' ) s[i]='\0';
    i++;
  }
  return(&s[1]);
}
      

void print_usage(void) {


  char msg[]="\
  --verbose    :increase level of messages printed to stderr\n\
  --brief      :decrease level of messages printed to stderr\n\
  --help       :print the usage\n\
  --read1 <filename> :fastq (optional gzipped) file name \n\
  --read2 <filename> :fastq (optional gzipped) file name \n\
  --index1 <filename> :fastq (optional gzipped) file name \n\
  --index2 <filename> :fastq (optional gzipped) file name \n\
  --index3 <filename> :fastq (optional gzipped) file name \n\
  --phred_encoding (33|64) :phred encoding used in the input files\n\
  --min_qual [0-40]        :defines the minimum quality that all bases in the UMI, CELL or Sample should have (reads that do not pass the criteria are discarded). 0 disables the filter. \n\
  --outfile1 <filename>    :file name for ouputing the reads from file1\n\
  --outfile2 <filename>    :file name for ouputing the reads from file2\n\
  --outfile3 <filename>    :file name for ouputing the reads from file3\n\
  --interleaved (read1|read2|index1|index2|index3),(read1|read2|index1|index2|index3)    :interleaved data\n\
  --umi_read (read1|read2|index1|index2|index3)       :in which input file can the UMI be found\n\
  --umi_offset integer     :offset \n\
  --umi_size               :number of bases after the offset\n\
  --cell_read (read1|read2|index1|index2|index3)      :in which input file can the cell be found\n\
  --cell_offset integer    :offset \n\
  --cell_size integer      :number of bases after the offset\n\
  --sample_read (read1|read2|index1|index2|index3)    :in which input file can the sample barcode be found\n\
  --sample_offset integer  :offset \n\
  --sample_size integer    :number of bases after the offset\n\
  --read1_offset integer   :\n\
  --read1_size integer     :\n\
  --read2_offset integer   :\n\
  --read2_size integer     :\n\
";
  fprintf(stderr,"usage: fastq_pre_barcodes --read1 fastq_file --outfile1 out_file [optional parameters]\n");
  fprintf(stderr,"%s\n",msg);
}
//
int main(int argc, char **argv ) {
  //
  Params *p=Params_new();
  int c;
  static int verbose=FALSE;
  static int paired=FALSE;
  static int help=FALSE;
  static int out_sam=FALSE;
  opterr = 0;

  fastq_print_version();
  
  static struct option long_options[] = {
    /* These options set a flag. */
    {"verbose", no_argument,       &verbose, TRUE},
    {"brief",   no_argument,       &verbose, FALSE},
    {"paired_end",   no_argument,       &paired, TRUE},
    {"single_end",   no_argument,       &paired, FALSE},
    {"sam",   no_argument,       &out_sam, TRUE},
    {"fastq",   no_argument,       &out_sam, FALSE},
    {"help",   no_argument, &help, TRUE},
    /* These options don’t set a flag.
       We distinguish them by their indices. */
    {"umi_read",     required_argument,       0, 'a'},
    {"umi_offset",  required_argument,       0, 'b'},
    {"umi_size",  required_argument, 0, 'c'},
    {"read1_offset",  required_argument, 0, 'd'},
    {"read1_size",    required_argument, 0, 'e'},
    {"read2_offset",  required_argument, 0, 'f'},
    {"read2_size",    required_argument, 0, 'g'},
    {"min_qual",     required_argument,       0, 'h'},
    {"cell_read",     required_argument,       0, 'i'},
    {"cell_offset",  required_argument,       0, 'j'},
    {"cell_size",  required_argument, 0, 'k'},
    {"read1",  required_argument, 0, 'l'},
    {"read2",  required_argument, 0, 'm'},
    {"index1",  required_argument, 0, 't'},
    {"index2",  required_argument, 0, 'v'},
    {"index3",  required_argument, 0, 'u'},
    {"outfile1",  required_argument, 0, 'n'},
    {"outfile2",  required_argument, 0, 'o'},
    {"interleaved",  required_argument, 0, 'z'},
    {"sample_read",     required_argument,       0, 'p'},
    {"sample_offset",  required_argument,       0, 'q'},
    {"sample_size",  required_argument, 0, 'r'},
    {"phred_encoding",  required_argument, 0, 's'},
    {0,0,0,0}
  };

  // process arguments
  while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    c = getopt_long (argc, argv, "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:z:",
		     long_options, &option_index);
    int xx=0;
    char *token;
    if (c == -1) // no more options
      break;
      
    switch (c) {

    case 'z':
      token = strtok(optarg,",");
      while( token != NULL ) {
	p->interleaved[xx]=read_index2read_idx(token);
	++xx;
	token = strtok(NULL, ",");
      }
      if ( xx!=2 ) {
	PRINT_ERROR("two file references should be passed to --interleaved");
      }
      p->has_interleaved_entries=TRUE;
      break;
    case 'a':
      p->umi_read=read_index2read_idx(optarg);
      break;
      
    case 'b':
      p->umi_offset=atol(optarg);
      break;
      
    case 'c':
      p->umi_size=atol(optarg);
      break;
	  
    case 'd':
      p->read_offset[READ1]=atol(optarg);
      break;

    case 'e':
      p->read_size[READ1]=atol(optarg);
      break;

    case 'f':
      p->read_offset[READ2]=atol(optarg);
      break;

    case 'g':
      p->read_size[READ2]=atol(optarg);
      break;
      
    case 'h':
      p->min_qual=atoi(optarg);
      break;

    case 'i':
      p->cell_read=read_index2read_idx(optarg);
      break;
      
    case 'j':
      p->cell_offset=atol(optarg);
      break;
      
    case 'k':
      p->cell_size=atol(optarg);
      break;

    case 'l':
      set_input_file(p,optarg,READ1);
      break;

    case 'm':
      set_input_file(p,optarg,READ2);
      break;

    case 't':
      set_input_file(p,optarg,READ3);
      break;

    case 'v':
      set_input_file(p,optarg,READ4);
      break;

    case 'u':
      set_input_file(p,optarg,READ5);
      break;

    case 'n':
      p->outfile[READ1]=optarg;
      break;

    case 'o': // only necessary for PE data
      p->outfile[READ2]=optarg;
      break;

      //    case 'u': 
      //      p->outfile[READ3]=optarg;
      //      break;

    case 'p':
      p->sample_read=read_index2read_idx(optarg);
      break;
      
    case 'q':
      p->sample_offset=atol(optarg);
      break;
      
    case 'r':
      p->sample_size=atol(optarg);
      break;

    case 's':
      p->phred_encoding=atoi(optarg);
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
  p->paired=paired;
  p->verbose=verbose;
  // 
  PRINT_INFO("Validating options...");
  validate_options(p);
  PRINT_INFO("Options OK.");
  // Assumptions:
  // 1) the fastq files have been validated and
  // 2) reads have the same order in the multiple files
  unsigned long discarded_reads=0;
  unsigned long processed_reads=0;

  // place holders to keep the info
  char sample[MAX_BARCODE_LENGTH];
  char umi[MAX_BARCODE_LENGTH];
  char cell[MAX_BARCODE_LENGTH];
  char sample_qual[MAX_BARCODE_LENGTH];
  char umi_qual[MAX_BARCODE_LENGTH];
  char cell_qual[MAX_BARCODE_LENGTH];

  PRINT_INFO("input files %d",p->num_input_files);
  // 
  short x;
  FASTQ_BOOLEAN skip;
  FASTQ_FILE* fdi[INDEX3+1]={NULL,NULL,NULL,NULL,NULL};
  FASTQ_ENTRY* m[INDEX3+1]={NULL,NULL,NULL,NULL,NULL};

  FASTQ_FILE* fdw[INDEX3+1]={NULL,NULL,NULL,NULL,NULL}; // out files

  char rnames[INDEX3+1][MAX_LABEL_LENGTH];
  
  for (x=READ1;x<=INDEX3;++x) {
    if ( p->file[x] != NULL ) {
      PRINT_VERBOSE(p,"Opening %s",p->file[x]);
      fdi[x]=fastq_new(p->file[x],FALSE,"r");
      fdi[x]->is_pe=TRUE;
    }
  }
  PRINT_VERBOSE(p,"done\n");
  
  // initialize the objects to hold an entry
  for (x=READ1;x<=INDEX3;++x) m[x]=fastq_new_entry();

  // output files
  if ( ! out_sam ) 
    for (x=READ1;x<=READ2;++x) {
      if ( p->outfile[x]!=NULL)
	fdw[x]=fastq_new(p->outfile[x],FALSE,"w4"); 
    }
  else {
    PRINT_SAM("@HD\tVN:1.0 SO:unknown\n");    
    PRINT_SAM("@PG\tID:1 PN:fastq_pre_barcodes CL:%s",argv[0]);
    for(c=1;c<argc-1;c++) {
      PRINT_SAM(" %s",argv[c]);
    }
    PRINT_SAM("\n");
  }

  //
  while(!fastq_files_eof(fdi,INDEX3) ) {
    skip=FALSE;

    for (x=READ1;x<=INDEX3;++x) // read one entry from each input  file
      if ( p->file[x] != NULL )
	if (fastq_read_next_entry(fdi[x],m[x])==0) goto end_loop; // EOF
    // handle interleaved
    if (p->has_interleaved_entries) {
      // jump to next read
      if (fastq_read_next_entry(fdi[p->interleaved[1]],m[p->interleaved[1]])==0) goto end_loop; // EOF
    }
    // check if the read names match
    if (p->num_input_files>1) {
      unsigned long len;
      for (x=READ1;x<=INDEX3;++x)
	if ( p->file[x] != NULL )
	  fastq_get_readname(fdi[x],m[x],&rnames[x][0],&len,TRUE);

      if ( p->file[READ2] != NULL )
	if (strcmp(rnames[READ1],rnames[READ2]) ) {
	  PRINT_ERROR("Readnames do not match across files (read #%ld)",processed_reads+1);
	  exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
	}

      if ( p->file[READ3] != NULL )
	if (strcmp(rnames[READ1],rnames[READ3]) ) {
	  PRINT_ERROR("Readnames do not match across files (read #%ld)",processed_reads+1);
	  exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
	}

      if ( p->file[READ4] != NULL )
	if (strcmp(rnames[READ1],rnames[READ4]) ) {
	  PRINT_ERROR("Readnames do not match across files (read #%ld)",processed_reads+1);
	  exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
	}

      if ( p->file[READ5] != NULL )
	if (strcmp(rnames[READ1],rnames[READ5]) ) {
	  PRINT_ERROR("Readnames do not match across files (read #%ld)",processed_reads+1);
	  exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
	}     
    }
    //
    ++processed_reads;
    // initialize
    cell[0]=sample[0]=umi[0]='\0';
    cell_qual[0]=sample_qual[0]=umi_qual[0]='\0';
    for (x=READ1;x<=INDEX3;++x) // extract the tags from each file
      if ( p->file[x]!=NULL) {
	if (extract_info(m[x],p,x,&cell[0],&umi[0],&sample[0],&cell_qual[0],&umi_qual[0],&sample_qual[0],out_sam)==FALSE ) {
	  // failed due to bad quality
	  // other errors would have resulted in the exit of the program
	  PRINT_VERBOSE(p,"Discarded %s %s %s <- %s\n",cell,umi,sample,m[x]->hdr1);
	  ++discarded_reads;
	  skip=TRUE;
	  break;
	}
      }
  
    if(skip) continue;
    // output the read(s) with the modified read name
    PRINT_VERBOSE(p,"_STAGS_CELL=%s_UMI=%s_SAMPLE=%s_ETAGS\n",cell,umi,sample);
    
    if ( out_sam ) {
      int se=(p->outfile[READ2]!=NULL?FALSE: TRUE);
      // prepare the "alignment"
      // slice_read(m[x],p,x);
      int flag=0;
      flag=BAM_FUNMAP; 
      if ( !se ) flag=BAM_FUNMAP|BAM_FMUNMAP|BAM_FPAIRED|BAM_FREAD1; // mate is unmapped
      
      // readname, quality and seq
      slice_read(m[READ1],p,READ1);
      // qname flag RNAME@* POS@0 MAPQ@255 CIGAR@* RNEXT(next mate*) PNEXT@0 TLEN@0 SEQ  QUAL Z:QX Z:OQ  ...      
      // rname = unsigned int
      unsigned int len=strlen(m[READ1]->seq);
      m[READ1]->seq[len-1]='\0';
      m[READ1]->qual[strlen(m[READ1]->qual)-1]='\0';
      char *rn1=format_read_name(m[READ1]->hdr1);
      PRINT_SAM("%ld\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%u",
	processed_reads,flag,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,len-1);


      PRINT_SAM("\t%s\t%s\t%s:Z:%s",m[READ1]->seq,&m[READ1]->qual[0],ORIG_RN_TAG,rn1);
      PRINT_SAM("\t%s:Z:%s",ORIG_QUAL_TAG,m[READ1]->qual);
      //
      if ( umi[0]!='\0'  ) 
        PRINT_SAM("\t%s:Z:%s\t%s:Z:%s",UMI_TAG,umi,UMI_QUAL_TAG,umi_qual);
      if ( cell[0]!='\0'  ) 
        PRINT_SAM("\t%s:Z:%s\t%s:Z:%s",CELL_TAG,cell,CELL_QUAL_TAG,cell_qual);
      if ( sample[0]!='\0'  ) 
        PRINT_SAM("\t%s:Z:%s\t%s:Z:%s",SAMPLE_TAG,sample,SAMPLE_QUAL_TAG,sample_qual);
    
      PRINT_SAM("\n");
    
      if ( !se ) {      // 2nd mate
        flag=BAM_FUNMAP|BAM_FMUNMAP|BAM_FPAIRED|BAM_FREAD2; // mate is unmapped
         slice_read(m[READ2],p,READ2);
	 // qname flag RNAME@* POS@0 MAPQ@255 CIGAR@* RNEXT(next mate*) PNEXT@0 TLEN@0 SEQ  QUAL Z:QX Z:OQ  ...      
	 // rname = unsigned int
	 len=strlen(m[READ2]->seq);
	 m[READ2]->seq[len-1]='\0';
	 m[READ2]->qual[strlen(m[READ2]->qual)-1]='\0';
	 char *rn=format_read_name(m[READ2]->hdr1);
	 PRINT_SAM("%ld\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%u",processed_reads,flag,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,len);
	 PRINT_SAM("\t%s\t%s\t%s:Z:%s",m[READ2]->seq,m[READ2]->qual,ORIG_RN_TAG,rn);
	 PRINT_SAM("\t%s:Z:%s",ORIG_QUAL_TAG,m[READ2]->qual);
	 //
	 if ( umi[0]!='\0'  ) 
	   PRINT_SAM("\t%s:Z:%s\t%s:Z:%s",UMI_TAG,umi,UMI_QUAL_TAG,umi_qual);
	 if ( cell[0]!='\0'   ) 
	   PRINT_SAM(" %s:Z:%s\t%s:Z:%s",CELL_TAG,cell,CELL_QUAL_TAG,cell_qual);
	 if ( sample[0]!='\0'  ) 
	   PRINT_SAM("\t%s:Z:%s\t%s:Z:%s",SAMPLE_TAG,sample,SAMPLE_QUAL_TAG,sample_qual);
	 
	 PRINT_SAM("\n");
      }
    } else {
      // fastq
      for (x=1;x<=READ2;++x)
	if (fdw[x]!=NULL) {
	  add_tags2readname(m[x],cell,umi,sample);
	  slice_read(m[x],p,x);	
	  fastq_write_entry(fdw[x],m[x]);
	}
    }
    PRINT_READS_PROCESSED(fdi[READ1]->cline/4,100000);
    // handle interleaved
    if (p->has_interleaved_entries) {
      // jump to next read
      if (fastq_read_next_entry(fdi[p->interleaved[0]],m[p->interleaved[0]])==0) goto end_loop; // EOF
    }

  }
 end_loop: 
  //   extract the info, change read name, trim the read, write
  PRINT_INFO("Reads processed: %ld",processed_reads);
  PRINT_INFO("Reads discarded: %ld",discarded_reads);
  if ( !out_sam )
    for (x=READ1;x<=READ2;++x)
      if (fdw[x]!=NULL) 
	fastq_destroy(fdw[x]);
  
  
  exit(0);
}

