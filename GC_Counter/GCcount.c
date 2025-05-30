#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_NAME_LEN 30
#define MAX_SEQ_LEN 200

typedef struct {
    size_t len;
    char name[MAX_NAME_LEN + 1];
    char seq[MAX_SEQ_LEN + 1];
} FastaSeq;


int main(int argc, char *argv[]){


    FastaSeq seq = {
        strlen(argv[2]),
        "",
        "",
    };

    strncpy(seq.name, argv[1], MAX_NAME_LEN - 1);
    strncpy(seq.seq, argv[2], MAX_SEQ_LEN - 1);


    
    double GCcount = 0;

    for(int i = 0; i < seq.len; i++){

        if(seq.seq[i] == 'G' | seq.seq[i] == 'C'){
            GCcount = GCcount + 1;
        }

    }


    printf("GC content is %f \n", GCcount / seq.len);



}