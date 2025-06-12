// Simple hash table implemented in C.

#include "ht.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdbool.h>
#include <stddef.h>
#include <errno.h>
#include <stdio.h>
#include <math.h>

// Hash table entry (slot may be filled or empty).
typedef struct {
    const char* key;  // key is NULL if this slot is empty
    void* value;
} ht_entry;

void exit_nomem(void){
    fprintf(stderr, "out of memory\n");
    exit(1);
}


// Hash table structure: create with ht_create, free with ht_destroy.
struct ht {
    ht_entry* entries;  // hash slots
    size_t capacity;    // size of _entries array
    size_t length;      // number of items in hash table
};

#define FNV_OFFSET 14695981039346656037UL
#define FNV_PRIME 1099511628211UL
#define INITIAL_CAPACITY 65536

static const char ALPHABET[5] = {'A', 'C', 'G', 'T', 'N'};

int main(int argc, char *argv[]){


    


    /* Make sure correct number of arguments */
    if (argc != 3){
        fprintf(stderr,
        "Usage: %s <read sequence> <k> \n"
        " <fastq> : Path to a .fastq file\n"
        " <k> : k-mer length to search for\n",
        argv[0]
        );
        return EXIT_FAILURE;
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


    /* Open file */

    const char *path = argv[1];

    FILE *file = fopen(path, "r");

    if(!file){
        perror("fopen");
        fprintf(stderr, "x Failed to open '%s'\n", path);
        return EXIT_FAILURE;
    }

    /* LOOP OVER FILE COUNT KMERS IN ALL READS */

    char *line = NULL;
    size_t n = 0;
    size_t lineno = 0;

    // Data structures used throughout
    ht* kcounts = ht_create();
    if(kcounts == NULL){
        exit_nomem();
    }
    char *kmer = malloc(k + 1);


    while(getline(&line, &n, file) != -1){

        if(lineno == 29100){

            printf("What's wrong?");
        }

        if(lineno % 4 != 1){

            lineno++;
            continue;
            

        }



        lineno++;


        /* Parse input sequence */
        line[strcspn(line, "\r\n")] = '\0';
        const char *seq = line;

        const char *seqcheck = line;

        size_t seqlen = strlen(seq);

        /* Check validity of input sequence */

        while(*seqcheck != '\0'){

            for(size_t i = 0; i < 5; i++){

                if(ALPHABET[i] == *seqcheck){
                    seqcheck++;
                    break;
                }else if(i == 4){

                    fprintf(stderr,
                    "Error: <read sequence> must be a sequence of As, Gs, Ts, Cs, and Ns\n");
                    return EXIT_FAILURE;

                }

            }

        }


        /* Count all kmers */

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
            

    }


    /* Free stuff */
    free(line);
    free(kmer);

    /*
    Print out kmers and frequencies, freeing values as we go.
    */

    hti it = ht_iterator(kcounts);
    while(ht_next(&it)){
        printf("%s %d\n", it.key, *(int*)it.value);
        free(it.value);
    }

    ht_destroy(kcounts);

}





ht* ht_create(void) {
    // Allocate space for hash table struct.
    ht* table = malloc(sizeof(ht));
    if (table == NULL) {
        return NULL;
    }
    table->length = 0;
    table->capacity = INITIAL_CAPACITY;

    // Allocate (zero'd) space for entry buckets.
    table->entries = calloc(table->capacity, sizeof(ht_entry));
    if (table->entries == NULL) {
        free(table); // error, free table before we return!
        return NULL;
    }
    return table;
}

void ht_destroy(ht* table) {
    // First free allocated keys.
    for (size_t i = 0; i < table->capacity; i++) {
        free((void*)table->entries[i].key);
    }

    // Then free entries array and table itself.
    free(table->entries);
    free(table);
}

#define FNV_OFFSET 14695981039346656037UL
#define FNV_PRIME 1099511628211UL

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

void* ht_get(ht* table, const char* key) {
    // AND hash with capacity-1 to ensure it's within entries array.
    uint64_t hash = hash_key(key);
    size_t index = (size_t)(hash & (uint64_t)(table->capacity - 1));

    // Loop till we find an empty entry.
    while (table->entries[index].key != NULL) {
        if (strcmp(key, table->entries[index].key) == 0) {
            // Found key, return value.
            return table->entries[index].value;
        }
        // Key wasn't in this slot, move to next (linear probing).
        index++;
        if (index >= table->capacity) {
            // At end of entries array, wrap around.
            index = 0;
        }
    }
    return NULL;
}

// Internal function to set an entry (without expanding table).
static const char* ht_set_entry(ht_entry* entries, size_t capacity,
        const char* key, void* value, size_t* plength) {
    // AND hash with capacity-1 to ensure it's within entries array.
    uint64_t hash = hash_key(key);
    size_t index = (size_t)(hash & (uint64_t)(capacity - 1));

    // Loop till we find an empty entry.
    while (entries[index].key != NULL) {
        if (strcmp(key, entries[index].key) == 0) {
            // Found key (it already exists), update value.
            entries[index].value = value;
            return entries[index].key;
        }
        // Key wasn't in this slot, move to next (linear probing).
        index++;
        if (index >= capacity) {
            // At end of entries array, wrap around.
            index = 0;
        }
    }

    // Didn't find key, allocate+copy if needed, then insert it.
    if (plength != NULL) {
        key = strdup(key);
        if (key == NULL) {
            return NULL;
        }
        (*plength)++;
    }
    entries[index].key = (char*)key;
    entries[index].value = value;
    return key;
}

// Expand hash table to twice its current size. Return true on success,
// false if out of memory.
static bool ht_expand(ht* table) {
    // Allocate new entries array.
    size_t new_capacity = table->capacity * 2;
    if (new_capacity < table->capacity) {
        return false;  // overflow (capacity would be too big)
    }
    ht_entry* new_entries = calloc(new_capacity, sizeof(ht_entry));
    if (new_entries == NULL) {
        return false;
    }

    // Iterate entries, move all non-empty ones to new table's entries.
    for (size_t i = 0; i < table->capacity; i++) {
        ht_entry entry = table->entries[i];
        if (entry.key != NULL) {
            ht_set_entry(new_entries, new_capacity, entry.key,
                         entry.value, NULL);
        }
    }

    // Free old entries array and update this table's details.
    free(table->entries);
    table->entries = new_entries;
    table->capacity = new_capacity;
    return true;
}

const char* ht_set(ht* table, const char* key, void* value) {
    assert(value != NULL);
    if (value == NULL) {
        return NULL;
    }

    // If length will exceed half of current capacity, expand it.
    if (table->length >= table->capacity / 2) {
        if (!ht_expand(table)) {
            return NULL;
        }
    }

    // Set entry and update length.
    return ht_set_entry(table->entries, table->capacity, key, value,
                        &table->length);
}

size_t ht_length(ht* table) {
    return table->length;
}

hti ht_iterator(ht* table) {
    hti it;
    it._table = table;
    it._index = 0;
    return it;
}

bool ht_next(hti* it) {
    // Loop till we've hit end of entries array.
    ht* table = it->_table;
    while (it->_index < table->capacity) {
        size_t i = it->_index;
        it->_index++;
        if (table->entries[i].key != NULL) {
            // Found next non-empty item, update iterator key and value.
            ht_entry entry = table->entries[i];
            it->key = entry.key;
            it->value = entry.value;
            return true;
        }
    }
    return false;
}