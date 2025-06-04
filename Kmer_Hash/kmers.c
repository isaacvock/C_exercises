#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>


/*
Strategy: Iterate over each kmer in a given read and
hash it to find which entry of the array needs to be
incremented. 
*/

static const char ALPHABET[4] = {'A', 'C', 'G', 'T'};

typedef struct {
    size_t count;
    char *kmer;
} kmer_ht;

void exit_nomem(void){
    fprintf(stderr, "out of memory\n");
    exit(1);
}


/* 
FUNCTIONS I WILL FIGURE OUT HOW TO IMPLEMENT:
*/

// Create kmer hash table
kmer_ht* ht_create(void);

// Get item with a given key (NUL-terminated) from
// hash table. Return the value (which was set with ht_set),
// or NULL if key not found
void* ht_get(kmer_ht* table, const char* key);

// Free memory allocated for hash table, including allocated keys
void ht_destroy(kmer_ht* table);

// Set item with given key to value. If not present in table,
// key is copied to newly allocated memory. Return address
// of copied key, or NULL if out of memory.
const char* ht_set(kmer_ht* table, const char* key, void* value);

// Hash table iterator: create with ht_iterator(), iterate with ht_next()
typedef struct{
    const char* key; // current key
    void* value; // current value

    // Don't use these fields directly
    kmer_ht* _table; // ref to hash table being iterated
    size_t _index; // current index into ht._entries
} kmer_hti;


// Return new hash table iterator (for use with ht_next)
kmer_hti ht_iterator(kmer_ht* table);

// Move iterator to next item in hash table, update iterator's key
// and value to current item, and return true. If there are
// no more items, return false. Don't call ht_set during iteration
bool ht_next(kmer_hti* it);

// Return number of items in hash table.
size_t ht_length(kmer_ht* table);


int main(int argc, char *argv[]){


    /* Make sure correct number of arguments */
    if (argc != 3){
        fprintf(stderr,
        "Usage: %s <read sequence> <k> \n"
        " <read sequence> : a sequence of As, Gs, Ts, and Cs\n"
        " <k> : k-mer length to search for\n",
        argv[0]
        );
        return EXIT_FAILURE;
    }

    /* Parse input sequence */

    const char *seq = argv[1];
    const char *seqcheck = argv[1];

    size_t seqlen = strlen(seq);

    /* Check validity of input sequence */

    while(*seqcheck != '\0'){

        for(size_t i = 0; i < 4; i++){

            if(ALPHABET[i] == *seqcheck){
                seqcheck++;
                break;
            }else if(i == 3){

                fprintf(stderr,
                "Error: <read sequence> must be a sequence of As, Gs, Ts, and Cs\n");
                return EXIT_FAILURE;

            }

        }

    }
    
    /* Parse k */

    char *endptr = NULL;
    unsigned long long tmp = strtoull(argv[2], &endptr, 10);

    errno = 0;

    if (errno == ERANGE ||           /* out of range for unsigned long long */
        *endptr != '\0' ||           /* trailing junk like "123abc"         */
        tmp == 0) {                  /* 0 is not a sensible k     */
        fprintf(stderr,
                "Error: <k> must be a positive integer between 1 and %zu\n",
                SIZE_MAX);
        return EXIT_FAILURE;
    }

    size_t k = (size_t)tmp;

    /* Count all kmers */

    size_t total_kmers = pow(4, k);


    kmer_ht* kcounts = ht_create();
    if(kcounts == NULL){
        exit_nomem()
    }


    char *kmer = malloc(k + 1);

    for(size_t i = 0; i < (seqlen - k); i++){

        strncpy(kmer, seq, k);
        seq++;

        void* value = ht_get(kcounts, kmer);
        if(value != NULL){

            int* pcount = (int*)value;
            (*pcount)++;
            continue;

        }

        // kmer not found, allocate space for new int and set to 1
        int* pcount = malloc(sizeof(int));
        if (pcount == NULL) {
            exit_nomem();
        }

        *pcount = 1;

        if(ht_set(kcounts, kmer, pcount) == NULL){
            exit_nomem();
        }



    }


    /*
    Print out kmers and frequencies, freeing values as we go.
    */

    kmer_hti it = ht_iterator(kcounts);
    while(ht_next(&it)){
        printf("%s %d\n", it.key, *(int*)it.value);
        free(it.value);
    }

    ht_destroy(kcounts);

}




// Create kmer hash table
kmer_ht* ht_create(void){

}

// Get item with a given key (NUL-terminated) from
// hash table. Return the value (which was set with ht_set),
// or NULL if key not found
void* ht_get(kmer_ht* table, const char* key){

}

// Free memory allocated for hash table, including allocated keys
void ht_destroy(kmer_ht* table){

}

// Set item with given key to value. If not present in table,
// key is copied to newly allocated memory. Return address
// of copied key, or NULL if out of memory.
const char* ht_set(kmer_ht* table, const char* key, void* value){

}


// Return new hash table iterator (for use with ht_next)
kmer_hti ht_iterator(kmer_ht* table){

}

// Move iterator to next item in hash table, update iterator's key
// and value to current item, and return true. If there are
// no more items, return false. Don't call ht_set during iteration
bool ht_next(kmer_hti* it){

}

// Return number of items in hash table.
size_t ht_length(kmer_ht* table){

}

