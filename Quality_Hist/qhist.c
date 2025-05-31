#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_READ_LEN 300
#define MAX_QNAME_LEN 100 

typedef struct{
    size_t read_len;
    char seq[MAX_READ_LEN + 1];
    char qualities[MAX_READ_LEN + 1];
    char name[MAX_QNAME_LEN + 1];
} fastq_entry;


int main(int argc, char *argv[]){

    /* Open file */

    const char *path = argv[1];

    FILE *file = fopen(path, "r");

    if (!file){
        perror("fopen");
        fprintf(stderr, "x Failed to open '%s'\n", path);
        return EXIT_FAILURE;
    }


    /* How many reads to read? */

    char *endptr = NULL;
    unsigned long long tmp = strtoull(argv[2], &endptr, 10);
    size_t num_reads = (size_t)tmp;

    /* Parse fastq entries */
    fastq_entry *fastqs;

    fastqs = malloc(num_reads * sizeof(fastq_entry));

    char line[1024];
    size_t lineno = 0;
    size_t entry_cnt = 0;

    while(fgets(line, sizeof line, file) && lineno < num_reads){
        
        
        if(lineno % 4 == 0){

            strncpy(fastqs[entry_cnt].name, line, MAX_QNAME_LEN);

        }else if(lineno % 4 == 1){

            strncpy(fastqs[entry_cnt].seq, line, MAX_READ_LEN);
            fastqs[entry_cnt].read_len = MAX_READ_LEN;

        }else if(lineno % 4 == 2){
            
            lineno++;
            continue;


        }else if(lineno % 4 == 3){

            strncpy(fastqs[entry_cnt].qualities, line, MAX_READ_LEN);
            entry_cnt++;

        }
        
        lineno++;

    }

    /* Print some stats */
    size_t readlen;
    float *GCconts;
    GCconts = malloc(lineno * sizeof(float));


    for(int i = 0; i<lineno; i++){

        float GCcount = 0;

        readlen = (float)fastqs[i].read_len;

        for(int j = 0; j < readlen; j++){

            if(fastqs[i].seq[j] == 'G' | fastqs[i].seq[j] == 'C'){
                GCcount = GCcount + 1;
            }

        }

        GCconts[i] = GCcount / readlen;

    }



    float avg_GC = 0;

    for(int i = 0; i <lineno; i++){

        avg_GC = avg_GC + (GCconts[i] / (lineno + 0.0));

    }


    printf("Average GC fraction is %.4f\n", avg_GC);
}