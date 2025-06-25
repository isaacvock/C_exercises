#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <stdint.h>
#include <omp.h>
#include <time.h>


typedef struct{
    size_t read_len;
    char *seq;
    char *qualities;
    char *name;
} fastq_entry;


int main(int argc, char *argv[]){

    clock_t start_time = clock();

    /* Make sure correct number of arguments */
    if (argc != 3){
        fprintf(stderr,
        "Usage: %s <fastq file> <num_reads> \n"
        " <fastq file> : path to a .fastq or .fq file\n"
        " <num_reads> : positive integer (how many reads to parse)\n",
        argv[0]
        );
        return EXIT_FAILURE;
    }

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

    errno = 0;

    if (errno == ERANGE ||           /* out of range for unsigned long long */
        *endptr != '\0' ||           /* trailing junk like "123abc"         */
        tmp == 0) {                  /* 0 is not a sensible read count      */
        fprintf(stderr,
                "Error: <num_reads> must be a positive integer between 1 and %zu\n",
                SIZE_MAX);
        return EXIT_FAILURE;
    }


    size_t num_reads = (size_t)tmp;

    /* Parse fastq entries */

    fastq_entry *fastqs;

    fastqs = malloc(num_reads * sizeof(fastq_entry));

    char *line = NULL;
    size_t n = 0;
    size_t lineno = 0;
    size_t entry_cnt = 0;

    while(getline(&line, &n, file) && lineno < (num_reads * 4)){
        
        line[strcspn(line, "\r\n")] = '\0';
        
        size_t mem_to_allocate = strlen(line) + 1;

        if(lineno % 4 == 0){

            fastqs[entry_cnt].name = malloc(mem_to_allocate);
            strcpy(fastqs[entry_cnt].name, line);

        }else if(lineno % 4 == 1){

            fastqs[entry_cnt].seq = malloc(mem_to_allocate);
            strcpy(fastqs[entry_cnt].seq, line);
            fastqs[entry_cnt].read_len = mem_to_allocate - 1;

        }else if(lineno % 4 == 3){

            fastqs[entry_cnt].qualities = malloc(mem_to_allocate);
            strcpy(fastqs[entry_cnt].qualities, line);
            entry_cnt++;

        }
        
        lineno++;
        n = 0;

    }

    /* Free stuff */
    free(line);

    /* Print some stats */
    size_t readlen;
    float *GCconts;
    GCconts = malloc(num_reads * sizeof(float));


    for(size_t i = 0; i<num_reads; i++){

        size_t GCcount = 0;

        readlen = fastqs[i].read_len;

        for(size_t j = 0; j < readlen; j++){

            if(toupper(fastqs[i].seq[j]) == 'G' || toupper(fastqs[i].seq[j]) == 'C'){
                GCcount = GCcount + 1;
            }

        }

        GCconts[i] = (double)GCcount / (double)readlen;

        free(fastqs[i].seq);
        free(fastqs[i].name);
        free(fastqs[i].qualities);

    }

    free(fastqs);

    double avg_GC = 0;

    for(size_t i = 0; i <num_reads; i++){

        avg_GC = avg_GC + (GCconts[i] / ((double)num_reads));

    }

    clock_t end_time = clock();

    printf("Average GC fraction is %.4f\n", avg_GC);

    double run_time = (double)(end_time - start_time)/CLOCKS_PER_SEC;
    printf("Run time is: %f\n", run_time);
}
