#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <stdint.h>
#include <math.h>


static const char ALPHABET[4] = {'A', 'C', 'G', 'T'};

typedef struct {
    size_t count;
    char *kmer;
} kmer_tab;


size_t count_kmers(const char *seq, char *query);

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

    kmer_tab *kmer_cnt;
    kmer_cnt = malloc(total_kmers * sizeof(kmer_tab));


    for(size_t i = 0; i < total_kmers; i++){


        /* Get the kmer of interest */
        char *s = malloc(k + 1);

        size_t tmp = i;
        for(int pos = k -1; pos >=0; --pos){

            /* Cutsey trick to output in lexicographical order */
            s[pos] = ALPHABET[tmp & 3];
            tmp >>= 2;

        }
        s[k] = '\0';
        
        kmer_cnt[i].kmer = malloc(k + 1);
        strcpy(kmer_cnt[i].kmer, s);


        /* Count instances of kmers */
        size_t kcount = count_kmers(seq, s);
        kmer_cnt[i].count = kcount;
        
        printf("%s appears %zu times\n", s, kcount);

    }


    /* Free memory */
    for(size_t i = 0; i < total_kmers; i++){

        free(kmer_cnt[i].kmer);

    }
    free(kmer_cnt);

}

size_t count_kmers(const char *seq, char *query){

    int count = 0;
    const char *tmp = seq;

    while(tmp = strstr(tmp, query)){

        count++;
        
        /* Move to next character of string */
        tmp++;

    }

    return count;
    
}

