/*
# =========================================================
# Copyright 2012-2018,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
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
#ifndef HASH
#define HASH
#include <stdlib.h>
#if defined (__cplusplus) || (defined (__STDC__) && __STDC__)
#define __ptr_t         void *
#else /* Not C++ or ANSI C.  */
#define __ptr_t         char *
#endif /* C++ or ANSI C.  */           

#ifndef unsignedlong
#define unsignedlong unsigned long long
#endif

#ifndef NULL
#define NULL    0
#endif  


struct bucket {
 struct bucket *next;
 unsignedlong value;      /* Value >=0 used as key in the hashing*/ 
 __ptr_t  obj;     /* pointer to a object*/
};
typedef struct bucket  hashnode;


struct hashtable_s {
  hashnode **buckets; //
  hashnode **buckets_last; // pointer to the last element of each bucket
  unsignedlong size;         // number of buckets
  unsignedlong last_bucket; // used in searchs/ hash traversals
  unsignedlong n_entries; // number of entries in the hashtable
  hashnode* last_node;
};

//typedef hashnode **hashtable;
typedef struct hashtable_s* hashtable;

/* functions */
hashtable new_hashtable(unsignedlong hashsize);
__ptr_t get_next_object(hashtable,unsignedlong);
__ptr_t delete(hashtable,unsignedlong,__ptr_t);
//__ptr_t replace_object(hashtable,unsignedlong,__ptr_t);
__ptr_t get_object(hashtable,unsignedlong);
int insere(hashtable,unsignedlong,__ptr_t);
void free_hashtable(hashtable);

void init_hash_traversal(hashtable table);
__ptr_t next_hash_object(hashtable table);
__ptr_t next_hashnode(hashtable table);
void hashtable_stats(hashtable table);
#endif
