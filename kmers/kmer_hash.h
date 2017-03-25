#ifndef KMER_HASH_H
#define KMER_HASH_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include "contig_generation.h"

/* Creates a hash table and (pre)allocates memory for the memory heap */
// hash_table_t* create_hash_table(int64_t nEntries, memory_heap_t *memory_heap)
// {
//    hash_table_t *result;
//    int64_t n_buckets = nEntries * LOAD_FACTOR;

//    result = (hash_table_t*) malloc(sizeof(hash_table_t));
//    result->size = n_buckets;
//    result->table = (bucket_t*) calloc(n_buckets , sizeof(bucket_t));
   
//    if (result->table == NULL) {
//       fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %lld buckets of %lu bytes\n", n_buckets, sizeof(bucket_t));
//       exit(1);
//    }
   
//    memory_heap->heap = (kmer_t *) malloc(nEntries * sizeof(kmer_t));
//    if (memory_heap->heap == NULL) {
//       fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
//       exit(1);
//    }
//    memory_heap->posInHeap = 0;
   
//    return result;
// }

/* Auxiliary function for computing hash values */
int64_t hashseq(int64_t  hashtable_size, char *seq, int size)
{
   unsigned long hashval;
   hashval = 5381;
   for(int i = 0; i < size; i++) {
      hashval = seq[i] +  (hashval << 5) + hashval;
   }
   
   return hashval % hashtable_size;
}

/* Returns the hash value of a kmer */
int64_t hashkmer(int64_t  hashtable_size, char *seq)
{
   return hashseq(hashtable_size, seq, KMER_PACKED_LENGTH);
}

/* Looks up a kmer in the hash table and returns a pointer to that entry */
shared kmer_t* lookup_kmer(int display, shared hash_table_t *hashtable, const unsigned char *kmer)
{
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(hashtable->size, (char*) packedKmer);
   bucket_t cur_bucket;
   shared kmer_t *result;
   
   cur_bucket = hashtable->table[hashval];
   result = cur_bucket.head;

   // if (display) {
   //    fprintf(stderr, "packedKmer: %s, localKmer: %s\n", packedKmer, localKmer);
   // }
   
   for (; result!=NULL; ) {

      char localKmer[KMER_PACKED_LENGTH];
      upc_memget(localKmer, result->kmer, KMER_PACKED_LENGTH * sizeof(unsigned char));

      if ( memcmp(packedKmer, localKmer, KMER_PACKED_LENGTH * sizeof(unsigned char)) == 0 ) {
         // fprintf(stderr, "Right character: %c\n", result->r_ext);
         // fprintf(stderr, "Left character: %c\n", result->l_ext);
         return result;
      }
      result = result->next;
   }
   return NULL;
}

/* Adds a kmer and its extensions in the hash table (note that a memory heap should be preallocated. ) */
int add_kmer(int64_t desired_position, shared hash_table_t *hashtable, shared memory_heap_t *memory_heap, shared [LINE_SIZE] const unsigned char *buffer, char left_ext, char right_ext)
{
   char bufferKmer[KMER_LENGTH];
   upc_memget(bufferKmer, buffer, KMER_LENGTH * sizeof(unsigned char));

   /* Pack a k-mer sequence appropriately */
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(bufferKmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(hashtable->size, (char*) packedKmer);
   int64_t pos = desired_position;
   
   /* Add the contents to the appropriate kmer struct in the heap */
   upc_memput((memory_heap->heap[pos]).kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(unsigned char));
   (memory_heap->heap[pos]).l_ext = left_ext;
   (memory_heap->heap[pos]).r_ext = right_ext;
   
   // fprintf(stderr, "HI\n");
   // // upc_lock_t* lock = hashtable->locks[hashval];
   // fprintf(stderr, "BYE\n");
   // // fprintf(stderr, "upc lock attempt: %d\n", upc_lock_attempt(hashtable->locks[hashval]));
   // fprintf(stderr, "hashval: %ld\n", hashval);
   // fprintf(stderr, "size: %ld\n", hashtable->size);
   
   upc_lock(hashtable->table[hashval].lock);
   // fprintf(stderr, "PASTLOCK");
   /* Fix the next pointer to point to the appropriate kmer struct */
   (memory_heap->heap[pos]).next = hashtable->table[hashval].head;
   /* Fix the head pointer of the appropriate bucket to point to the current kmer */
   hashtable->table[hashval].head = &(memory_heap->heap[pos]);
   upc_unlock(hashtable->table[hashval].lock);
   // fprintf(stderr, "HELLO");
   
   return 0;
}

/* Adds a k-mer in the start list from desired position of the memory heap (the k-mer was "just added" in the memory heap at position posInHeap - 1) */
void addKmerToStartList(int64_t desired_position, shared memory_heap_t *memory_heap, start_kmer_t **startKmersList)
{
   start_kmer_t *new_entry;
   shared kmer_t *ptrToKmer;
   
   int64_t prevPosInHeap = desired_position;
   ptrToKmer = &(memory_heap->heap[prevPosInHeap]);
   new_entry = (start_kmer_t*) malloc(sizeof(start_kmer_t));
   new_entry->next = (*startKmersList);
   new_entry->kmerPtr = ptrToKmer;
   (*startKmersList) = new_entry;
}

/* Deallocation functions */
int dealloc_heap(memory_heap_t *memory_heap)
{
   upc_free(memory_heap->heap);
   return 0;
}

int dealloc_hashtable(hash_table_t *hashtable)
{
   upc_free(hashtable->table);
   return 0;
}


#endif // KMER_HASH_H
