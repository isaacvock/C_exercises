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

#define FNV_OFFSET 14695981039346656037UL
#define FNV_PRIME 1099511628211UL

#define INITIAL_CAPACITY 65536

// Return 64-bit FNV-1a hash for key (NUL-terminated). See description:
// https://en.wikipedia.org/wiki/Fowler–Noll–Vo_hash_function
static uint64_t hash_key(const char* key) {
    uint64_t hash = FNV_OFFSET;
    for (const char* p = key; *p; p++) {
        hash ^= (uint64_t)(unsigned char)(*p);
        hash *= FNV_PRIME;
    }
    return hash;
}

static const char ALPHABET[4] = {'A', 'C', 'G', 'T'};

typedef struct{
    const char* key;
    void* value;
} ht_entry;

typedef struct {
    size_t length;
    size_t capacity;
    ht_entry* entries; // Hash slots
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
        exit_nomem();
    }


    char *kmer = malloc(k + 1);

    for(size_t i = 0; i < (seqlen - k); i++){

        strncpy(kmer, seq, k);
        seq++;

        void* value = ht_get(kcounts, kmer);
        if(value != NULL){

            // Increment value
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

    // Allocate space for hash table struct.
    kmer_ht *new_ht; 
    new_ht = malloc(sizeof(kmer_ht));

    if(new_ht == NULL){
        return NULL;
    }

    new_ht->length = 0;
    new_ht->capacity = INITIAL_CAPACITY;

    // Allocate (zero'd) space for entry buckets
    new_ht->entries = calloc(new_ht->capacity, sizeof(ht_entry));
    if(new_ht == NULL){
        free(new_ht);
        return NULL;
    }

    return new_ht;

}

// Get item with a given key (NUL-terminated) from
// hash table. Return the value (which was set with ht_set),
// or NULL if key not found
void* ht_get(kmer_ht* table, const char* key){

    uint64_t hash = hash_key(key);

    // Modulo for nerds (if capacity is power of 2)
    size_t index = (size_t)(hash & (uint64_t)(table->capacity - 1));


    // Loop until we find an empty entry, because that
    // will let us know that the key does not exist since
    // we use linear probing when inserting new objects 
    // into the hash table
    while(table->entries[index].key != NULL){

        if(strcmp(key, table->entries[index].key) == 0){

            // Found key, return value.
            return table->entries[index].value;

        }

        index++;
        if(index >= table->capacity){
            //Go back to beginning
            index = 0;
        }

    }

    return NULL;



}

// Free memory allocated for hash table, including allocated keys
void ht_destroy(kmer_ht* table){

    // Free allocated keys
    for(size_t i = 0; i < table->capacity; i++){
        // Need to cast because it is declared as a const char*
        // rather than just a char*
        free((void*)table->entries[i].key);
    }

    // Free entries array and table itself
    free(table->entries);
    free(table);

}


// Better ht_set implementation
const char* ht_set(kmer_ht* table, const char* key, void* value){


    // First check if length is 1/2 of capacity
    // If so, expand the hash table capacity
    if(table->length + 1 > table->capacity / 2){

        ht_expand(table);

    }


    // Then, hash and check where the key needs to be inserted.
    // If that space is taken, perform linear probing until
    // an open space is found

    uint64_t hash = hash_key(key);
    size_t index = (size_t)(hash & (uint64_t)(table->capacity - 1));

    while(table->entries[index].key != NULL){
        index++;

        if(index >= table->capacity){
            index = 0;
        }
    }

    strcpy(table->entries[index].key, key);
    table->entries[index].value = 1;


}


void ht_expand(kmer_ht* table){

    uint64_t hash;
    size_t index;

    // Double capacity size; inefficient but makes hashing
    // a lot more efficient because of the trick we can employ
    size_t old_capacity = table->capacity;
    table->capacity = table->capacity * 2;
    ht_entry* new_entries = calloc(table->capacity, sizeof(ht_entry));

    // Rehash everything
    for(size_t i = 0; i<old_capacity; i++){

        hash = hash_key(table->entries[i].key);
        index = (size_t)(hash & (uint64_t)(table->capacity - 1));

        while(new_entries[index].key != NULL){
            index++;

            if(index >= table->capacity){
                index = 0;
            }
        }

        strcpy(new_entries[index].key, table->entries[i].key);
        new_entries[index].value = table->entries[i].value;
        free((void*)table->entries[i].key);
        free(table->entries);

    }

    

}



// Set item with given key to value. If not present in table,
// key is copied to newly allocated memory. Return address
// of copied key, or NULL if out of memory.
const char* ht_set(kmer_ht* table, const char* key, void* value){

    uint64_t hash = hash_key(key);

    // Modulo for nerds
    size_t index = (size_t)(hash & (uint64_t)(table->capacity - 1));

    if(table->entries[index].key == '\0'){

        // Not present; add it
        table->entries[index].key = malloc(sizeof(key));
        strcpy(table->entries[index].key, key);
        table->entries[index].value = value;

        return &table->entries[index].key;

    }else if(strcmp(table->entries[index].key, key) == 0){

        // Present; return key
        return &table->entries[index].key;

    }else{

        // Find an empty element of entries array
        size_t init_index = index;
        index++;
        while(table->entries[index].key != '\0'){

            index++;

            if(index >= table->capacity){
                index = 0;
            }

            // If we have looped through everything, we know the
            // array is full and needs to be resized
            if(index == init_index){
                
                size_t init_capacity = table->capacity;

                // If capacity is not a power of 2, our fast modulo alternative will break
                table->capacity = table->capacity * 2;
                void *tmp = realloc(table->entries, table->capacity * sizeof(ht_entry));

                if(tmp == NULL){
                    return NULL;
                }

                table->entries = tmp;

                strcpy(table->entries[init_capacity].key, key);

                return &table->entries[init_capacity].key;

            }

        }

        return &table->entries[index].key;




    }

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

