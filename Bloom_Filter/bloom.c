#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <errno.h>

/*
This is my attempt to implement a Bloom filter. Will use it to
determine whether or not there are kmers in a fastq file that
are not present in a genome.
A bloom filter uses several hash functions to identify elements
that are definitely not in a given set (though it is possible
for there to be false positive hits).
*/

// Hash table entry (slot may be filled or empty).
typedef struct{
    const char*  key; // key is NULL if this slot is empty
    void* value;
} ht_entry;

void exit_nomem(void){
    fprintf(stderr, "out of memory\n");
    exit(1);
}

// Hash table structure: create with ht_create, free with ht_destroy.
struct ht{
    ht_entry* entries; // hash slots
    size_t capacity; // size of _entries array
    size_t length; // number of items in hash table
};


#define FNV_OFFSET 14695981039346656037UL
#define FNV_PRIME 1099511628211UL
#define INITIAL_CAPACITY 65536

static const char ALPHABET[5] = {'A', 'C', 'G', 'T', 'N'};

int main(int argc, char* argv[]){

    /* Need two arguments */
    if(argc != 4){
        fprintf(stderr,
        "Usage: %s <FASTA file> <fastq file> <k> \n"
        " <FASTA file> : path to a .fasta or .fa file\n"
        " <fastq file> : path to a .fastq or .fq file"
        " <k> : length of k-mer to look for",
        argv[0]
        );
        return EXIT_FAILURE;
    }

    /* Open FASTA file */

    const char *fasta_path = argv[1];

    FILE *fasta_file = fopen(fasta_path, "r");

    if(!fasta_file){
        perror("fopen");
        fprintf(stderr, "x Failed to open '%s'\n", fasta_path);
        return EXIT_FAILURE;
    }

    /* Open FASTQ file */

    const char *fastq_path = argv[2];

    FILE *fastq_file = fopen(fastq_path, "r");

    if(!fastq_file){
        perror("fopen");
        fprintf(stderr, "x Failed to open '%s'\n", fastq_path);
        return EXIT_FAILURE;
    }

    /* Parse k */

    char *endptr = NULL;
    unsigned long long tmp = strtoull(argv[3], &endptr, 10);

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




}