#include <stdio.h>
//#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include "packingDNAseq.h"
#include "kmer_hash.h"

shared hash_table_t * hashtable;
shared memory_heap_t memory_heap;
upc_lock_t* shared *lock_array;


/* Creates a hash table and (pre)allocates memory for the memory heap */


int main(int argc, char **argv) {
    time_t start, end;
    double constrTime, traversalTime;
    char cur_contig[MAXIMUM_CONTIG_SIZE], unpackedKmer[KMER_LENGTH + 1], left_ext, right_ext, *input_UFX_name;
    int64_t posInContig, contigID = 0, totBases = 0, ptr = 0, nKmers, cur_chars_read, total_chars_to_read;
    unpackedKmer[KMER_LENGTH] = '\0';
    kmer_t *cur_kmer_ptr;
    start_kmer_t *startKmersList = NULL, *curStartNode;
    unsigned char *working_buffer;
    FILE *inputFile, *serialOutputFile;
    /* Read the input file name */
    input_UFX_name = argv[1];
    printf("Start reading!\n");
    /* ============== GRAPH CONSTRUCTION ============== */

    start = clock();
    /* Initialize lookup table that will be used for the DNA packing routines */
    init_LookupTable();
    printf("Start 00000\n");
    /* Extract the number of k-mers in the input file */
    upc_lock_t * file_lock =upc_all_lock_alloc();
    upc_lock(file_lock);
    printf("end 00000\n");
    nKmers = getNumKmersInUFX(input_UFX_name);
    printf("start 11111!\n");
    upc_unlock(file_lock);
    printf("Finish 11111!\n");

    /* Read the kmers from the input file and store them in the working_buffer */

    const int Kmers_per_thread = (nKmers + THREADS - 1) / THREADS;
    // read the Kumer info into local threads
    const int Kmer_local = (MYTHREAD == THREADS - 1) ? nKmers - (THREADS - 1) * Kmers_per_thread : Kmers_per_thread;

    printf("Finish 2222!\n");
    /* Read the kmers from the input file and store them in the working_buffer */
    total_chars_to_read = nKmers * LINE_SIZE;
    unsigned char * full_buffer = (unsigned char *) malloc(total_chars_to_read * sizeof(unsigned char));
    upc_lock(file_lock);
    inputFile = fopen(input_UFX_name, "r");
    cur_chars_read = fread(full_buffer, sizeof(unsigned char), total_chars_to_read, inputFile);
    fclose(inputFile);
    upc_unlock(file_lock);
    printf("start 3333!\n");
    working_buffer = (unsigned char *) malloc(Kmer_local*LINE_SIZE * sizeof(unsigned char));
    printf("kerms local: %d\n",Kmer_local);
    memcpy(working_buffer,&full_buffer[MYTHREAD*Kmers_per_thread*LINE_SIZE],Kmer_local*LINE_SIZE);
    printf("Finish 3333!\n");
    upc_barrier;

    /* Create a hash table */
    //local nkerms
    hashtable = create_hash_table(nKmers,&memory_heap);
    printf("after create, hashtable->size is %d\n",hashtable->size);
    upc_barrier;
    /* Process the working_buffer and store the k-mers in the hash table */
    /* Expected format: KMER LR ,i.e. first k characters that represent the kmer, then a tab and then two chatacers, one for the left (backward) extension and one for the right (forward) extension */
    printf("testing!!!!!\n");
    int64_t table_size = hashtable->size;
    lock_array =upc_all_alloc(table_size,sizeof(upc_lock_t*));
    //for (int i = 0; i < table_size; i++) {
    //    //printf("%d\n",i);
    //    lock_array[i] = upc_all_lock_alloc();
    //}
    int i;
    upc_forall(i=0;i<table_size;i++;i){
        lock_array[i]= upc_all_lock_alloc();
    }
    upc_barrier;
    printf("Start 4444!\n");

    printf("In main, hashtable->size is %d\n",hashtable->size);

    while (ptr < Kmer_local*LINE_SIZE) {
        /* working_buffer[ptr] is the start of the current k-mer                */
        /* so current left extension is at working_buffer[ptr+KMER_LENGTH+1]    */
        /* and current right extension is at working_buffer[ptr+KMER_LENGTH+2]  */

        left_ext = (char) working_buffer[ptr + KMER_LENGTH + 1];
        right_ext = (char) working_buffer[ptr + KMER_LENGTH + 2];
        printf("Testing\n");
        printf("ptr is %d, limit is %d\n",ptr,Kmer_local*LINE_SIZE);

        /* Pack a k-mer sequence appropriately */
        char packedKmer[KMER_PACKED_LENGTH];
        printf("T1\n");
        packSequence(&working_buffer[ptr], (unsigned char *) packedKmer, KMER_LENGTH); // &working_buffer[ptr] = kmer
        printf("%d, %s\n",hashtable->size,packedKmer);
        int64_t hashval = hashkmer(hashtable->size, (char *) packedKmer);
        printf("Start add_kmer\n");
        upc_lock(lock_array[hashval]);
        /* Add k-mer to hash table */
        add_kmer(hashtable, &memory_heap, &working_buffer[ptr], left_ext, right_ext, hashval);
        upc_unlock(lock_array[hashval]);
        printf("End add_kmer\n");
        /* Create also a list with the "start" kmers: nodes with F as left (backward) extension */
        if (left_ext == 'F') {
            addKmerToStartList(&memory_heap, &startKmersList);
        }

        /* Move to the next k-mer in the input working_buffer */
        ptr += LINE_SIZE;
    }
    printf("Finish 55555!\n");

    upc_barrier;
    end = clock();
    constrTime = 1.0 * (end - start) / CLOCKS_PER_SEC;

    /* ============== GRAPH TRAVERSAL ============== */

    start = clock();
    char filename[20];
    sprintf(filename,"serial%d.out",MYTHREAD);
    serialOutputFile = fopen(filename, "w");

    /* Pick start nodes from the startKmersList */
    curStartNode = startKmersList;

    while (curStartNode != NULL) {
        /* Need to unpack the seed first */
        cur_kmer_ptr = curStartNode->kmerPtr;
        unpackSequence((unsigned char *) cur_kmer_ptr->kmer, (unsigned char *) unpackedKmer, KMER_LENGTH);
        /* Initialize current contig with the seed content */
        memcpy(cur_contig, unpackedKmer, KMER_LENGTH * sizeof(char));
        posInContig = KMER_LENGTH;
        right_ext = cur_kmer_ptr->r_ext;

        /* Keep adding bases while not finding a terminal node */
        while (right_ext != 'F') {
            cur_contig[posInContig] = right_ext;
            posInContig++;
            /* At position cur_contig[posInContig-KMER_LENGTH] starts the last k-mer in the current contig */
            cur_kmer_ptr = lookup_kmer(hashtable, (const unsigned char *) &cur_contig[posInContig - KMER_LENGTH]);
            right_ext = cur_kmer_ptr->r_ext;
        }

        /* Print the contig since we have found the corresponding terminal node */
        cur_contig[posInContig] = '\0';
        fprintf(serialOutputFile, "%s\n", cur_contig);
        contigID++;
        totBases += strlen(cur_contig);
        /* Move to the next start node in the list */
        curStartNode = curStartNode->next;
    }

    fclose(serialOutputFile);
    end = clock();

    /* Print timing and output info */
    printf("Generated %lld contigs with %lld total bases\n", contigID, totBases);
    traversalTime = 1.0 * (end - start) / CLOCKS_PER_SEC;
    printf("Total execution time: %f seconds (%f graph construction / %f graph traversal)\n",
           constrTime + traversalTime, constrTime, traversalTime);

    return 0;
}
