/*
 * =========================================================
 * Copyright 2012-2018,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
 * You should have received a copy of the GNU General Public License
 * if not, see <http://www.gnu.org/licenses/>.
 *
 *
 * =========================================================
 */
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "fastq.h"
#include "hash.h"
#include "range_list.h"


int main(int argc, char **argv ) {

  fastq_print_version();
  
  //if (argc!=3) {
  //  fprintf(stderr,"Usage: fastq_tests\n");
  //  exit(PARAMS_ERROR_EXIT_STATUS);
  // }
  void *p;
  int v1=1,v2=2;
  hashtable ht=new_hashtable(1000);
  assert(ht!=NULL);
  hashtable_stats(ht);
  reset_hashtable(ht);
  //
  p=get_object(ht,100);
  assert(p==NULL);
  insere(ht,100,&v1);
  p=get_object(ht,100);
  assert(p!=NULL);
  insere(ht,100,&v1);
  insere(ht,110,&v2);
  hashtable_stats(ht);
  //
  p=get_object(ht,100);
  assert(p!=NULL);
  p=delete(ht,100,&v1);
  assert(p!=NULL);
  p=delete(ht,100,&v1);
  assert(p!=NULL);
  // 
  init_hash_traversal(ht);
  p=next_hash_object(ht);
  assert(p!=NULL);
  int *i=(int*)p;
  printf(">%d\n",*i);
  p=next_hash_object(ht);
  assert(p==NULL);
  // stats
  hashtable_stats(ht);
  insere(ht,110,&v2);
  //
  free_hashtable(ht);
  exit(0);
}

