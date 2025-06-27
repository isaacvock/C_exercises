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

    double avg_GC;
    double GCcont_track = 0;
    size_t GCcount = 0;
    size_t readlen;
    double total_GC = 0;

    #pragma omp parallel for    \
     reduction(+:total_GC)
    for(int i = 0; i < (num_reads*4); i++){

        if(getline(&line, &n, file)){

            line[strcspn(line, "\r\n")] = '\0';

            readlen = strlen(line);
            
            if(lineno % 4 == 1){

                for(size_t j = 0; j < readlen; j++){

                    char current_char = toupper(line[j]);
                    if(current_char == 'g' || current_char == 'c' || current_char == 'G' || current_char == 'C'){
                        GCcount = GCcount + 1;
                    }

                }

                GCcont_track = ((double)GCcount / (double)readlen);
                total_GC += GCcont_track;

            }

            entry_cnt++;
            lineno++;

        }

    }

    /* There were less than the requested number of reads*/
    if(entry_cnt < num_reads){
        num_reads = entry_cnt;
    }

    avg_GC = total_GC / num_reads;

    double end_time = omp_get_wtime();

    printf("Average GC fraction is %.4f\n", avg_GC);

    printf("Run time is: %.3f s\n", end_time - start_time);
}
