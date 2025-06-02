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

    /* Parse input sequence */

    const char *seq = argv[1];
    
    /* Parse k */

    char *endptr = NULL;
    unsigned long long tmp = strtoull(argv[2], &endptr, 10);
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