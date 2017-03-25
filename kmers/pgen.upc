#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <upc_io.h>

#include "packingDNAseq.h"
#include "kmer_hash.h"

shared int64_t posInContig, contigID = 0, totBases = 0, nKmers, cur_chars_read, total_chars_to_read;
shared [LINE_SIZE] unsigned char *working_buffer;
shared hash_table_t* hashtable;
shared memory_heap_t memory_heap;

int main(int argc, char *argv[]){

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;
	char cur_contig[MAXIMUM_CONTIG_SIZE], unpackedKmer[KMER_LENGTH+1], left_ext, right_ext, *input_UFX_name;
	unpackedKmer[KMER_LENGTH] = '\0';
	start_kmer_t *startKmersList = NULL, *curStartNode;
	FILE *parallelOutputFile;

    /* Read the input file name */
    input_UFX_name = argv[1];

	/** Read input **/
	upc_barrier;
	inputTime -= gettime();
	/* Initialize lookup table that will be used for the DNA packing routines */
    init_LookupTable(); 

    /* Extract the number of k-mers in the input file */
    nKmers = getNumKmersInUFX(input_UFX_name);

    /* Create a hash table */
    int64_t n_buckets = nKmers * LOAD_FACTOR;

    hashtable = (shared hash_table_t*) upc_all_alloc(1, sizeof(hash_table_t));
    hashtable->size = n_buckets;
    hashtable->table = (shared [1] bucket_t*) upc_all_alloc(n_buckets, sizeof(bucket_t));
    upc_barrier;

    int64_t i;
    for(i=0; i<n_buckets; i++) {
    	hashtable->table[i].lock = upc_global_lock_alloc();
      if (i == n_buckets-1) {
        fprintf(stderr, "HIII\n");
        fprintf(stderr, "upc lock attempt for %ld: %d\n", i, upc_lock_attempt(hashtable->table[i].lock));
        upc_unlock(hashtable->table[i].lock);
        fprintf(stderr, "unlock success\n");
        fprintf(stderr, "upc lock attempt for %ld: %d\n", i, upc_lock_attempt(hashtable->table[i].lock));
        upc_unlock(hashtable->table[i].lock);
        fprintf(stderr, "unlock success 2\n");
      }
    }
    i = n_buckets-1;
    fprintf(stderr, "hi");
    fprintf(stderr, "upc lock attempt for %ld: %d\n", i, upc_lock_attempt(hashtable->table[i].lock));
    upc_unlock(hashtable->table[i].lock);
    fprintf(stderr, "hiiii");
    upc_barrier;
    i = 3787674;
    fprintf(stderr, "upc lock attempt for %ld: %d\n", i, upc_lock_attempt(hashtable->table[i].lock));
    upc_unlock(hashtable->table[i].lock);

    if (hashtable->table == NULL) {
       fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %lld buckets of %lu bytes\n", n_buckets, sizeof(bucket_t));
       exit(1);
    }
    
    memory_heap.heap = (shared kmer_t *) upc_all_alloc(nKmers, sizeof(kmer_t));
    if (memory_heap.heap == NULL) {
       fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
       exit(1);
    }
    // done creating hash table

    i = 3787674;
    fprintf(stderr, "upc lock attempt for %ld: %d\n", i, upc_lock_attempt(hashtable->table[i].lock));
    upc_unlock(hashtable->table[i].lock);

    total_chars_to_read = nKmers * LINE_SIZE;
    upc_file_t *inputFile = upc_all_fopen(input_UFX_name, UPC_RDONLY|UPC_COMMON_FP, 0, NULL);
    working_buffer = (shared [LINE_SIZE] unsigned char*) upc_all_alloc(nKmers, LINE_SIZE * sizeof(unsigned char));
    cur_chars_read = upc_all_fread_shared(inputFile, (shared void*) working_buffer, LINE_SIZE, sizeof(unsigned char), total_chars_to_read, UPC_IN_ALLSYNC|UPC_OUT_ALLSYNC);
    upc_all_fclose(inputFile);

	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();

	int64_t ptr;
    upc_forall(ptr = 0; ptr < cur_chars_read; ptr += LINE_SIZE; &working_buffer[ptr]) {
    	left_ext = (char) working_buffer[ptr+KMER_LENGTH+1];
        right_ext = (char) working_buffer[ptr+KMER_LENGTH+2];  
        int64_t ctr = ptr / LINE_SIZE;

        /* Add k-mer to hash table */
        add_kmer(ctr, hashtable, &memory_heap, &working_buffer[ptr], left_ext, right_ext);  

        /* Create also a list with the "start" kmers: nodes with F as left (backward) extension */
        /* Each thread has its own linked list of start kmers */
        if (left_ext == 'F') {
           addKmerToStartList(ctr, &memory_heap, &startKmersList);
        } 
    }
	upc_barrier;
	constrTime += gettime();

	/** Graph traversal **/
	traversalTime -= gettime();
	char *parallelOutputFileName;
	shared kmer_t *cur_kmer_ptr;
	sprintf(parallelOutputFileName, "thread%d.out", MYTHREAD);
	parallelOutputFile = fopen(parallelOutputFileName, "w"); 

    /* Pick start nodes from the startKmersList */
    curStartNode = startKmersList; 

    while (curStartNode != NULL ) {
       /* Need to unpack the seed first */
       cur_kmer_ptr = curStartNode->kmerPtr;
       unsigned char* localKmer;
       upc_memget(localKmer, cur_kmer_ptr, KMER_PACKED_LENGTH * sizeof(unsigned char));
       unpackSequence((unsigned char*) localKmer,  (unsigned char*) unpackedKmer, KMER_LENGTH);
       /* Initialize current contig with the seed content */
       memcpy(cur_contig, unpackedKmer, KMER_LENGTH * sizeof(char));
       posInContig = KMER_LENGTH;
       right_ext = cur_kmer_ptr->r_ext; 

       /* Keep adding bases while not finding a terminal node */
       while (right_ext != 'F') {
          cur_contig[posInContig] = right_ext;
          posInContig++;
          /* At position cur_contig[posInContig-KMER_LENGTH] starts the last k-mer in the current contig */
          cur_kmer_ptr = lookup_kmer(hashtable, (const unsigned char *) &cur_contig[posInContig-KMER_LENGTH]);
          right_ext = cur_kmer_ptr->r_ext;
       } 

       /* Print the contig since we have found the corresponding terminal node */
       cur_contig[posInContig] = '\0';
       fprintf(parallelOutputFile,"%s\n", cur_contig);
       contigID++;
       totBases += strlen(cur_contig);
       /* Move to the next start node in the list */
       curStartNode = curStartNode->next;
    } 

    fclose(parallelOutputFile);
	upc_barrier;
	traversalTime += gettime();

	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
