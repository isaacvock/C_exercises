#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <stdint.h>
#include <omp.h>
#include <time.h>


int main(int argc, char *argv[]){

    double start_time = omp_get_wtime();

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

    char *line = NULL;
    size_t n = 0;
    size_t lineno = 0;
    size_t entry_cnt = 0;

    double avg_GC = 0.0;
    size_t GCcount = 0;
    size_t readlen;

    char **seq = malloc(num_reads * sizeof *seq);
    size_t *len = malloc(num_reads * sizeof *len);

    ssize_t l;

    for(size_t i = 0; i < num_reads; i++){

        getline(&line, &n, file); // Read name

        l = getline(&line, &n, file); // Sequence

        line[strcspn(line, "\r\n")] = '\0';

        seq[i] = strdup(line);
        len[i] = l - 1;

        entry_cnt++;
        lineno++;

        getline(&line, &n, file); // '+'
        getline(&line, &n, file); // Qualities

    }

    /* There were less than the requested number of reads*/
    if(entry_cnt < num_reads){
        num_reads = entry_cnt;
    }

    /* GC counting */
    #pragma omp parallel for reduction(+:avg_GC)
    for(size_t i = 0; i < num_reads; i++){

        GCcount = 0;
        readlen = len[i];
        for(size_t j = 0; j < readlen; j++){

            if(seq[i][j] == 'G' || seq[i][j] == 'C' || seq[i][j] == 'g' || seq[i][j] == 'c'){
                GCcount++;
            }

        }

        avg_GC += ((double)GCcount / (double)readlen)/((double)num_reads);

    }

    double end_time = omp_get_wtime();

    printf("Average GC fraction is %.4f\n", avg_GC);

    printf("Run time is: %.3f s\n", end_time - start_time);
}
