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

#define MATCH_SCORE 2
#define MISMATCH_SCORE -1
#define GAP_SCORE -2

static const char ALPHABET[4] = {'A', 'C', 'G', 'T'};

typedef struct{
    char* align1;
    char* align2;
} alignment;


void exit_nomem(void){
    fprintf(stderr, "out of memory\n");
    exit(1);
}

int main(int argc, char *argv){


    /* INPUT PARSING */

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

    /* Parse 1st input sequence */

    const char *seq1 = argv[1];
    const char *seqcheck = argv[1];

    size_t seqlen = strlen(seq1);

    /* Check validity of 1st input sequence */

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

    /* Parse 2nd input sequence */

    const char *seq2 = argv[2];
    const char *seqcheck2 = argv[2];

    size_t seqlen = strlen(seq2);

    /* Check validity of input sequence */

    while(*seqcheck2 != '\0'){

        for(size_t i = 0; i < 4; i++){

            if(ALPHABET[i] == *seqcheck2){
                seqcheck2++;
                break;
            }else if(i == 3){

                fprintf(stderr,
                "Error: <read sequence> must be a sequence of As, Gs, Ts, and Cs\n");
                return EXIT_FAILURE;

            }

        }

    }

    /* DP */

    alignment* align = get_alignment(seq1, seq2);

}


/* Run N-W algorithm*/
alignment* get_alignment(char* s1, char* s2){

    /* Copy strings to iterate along */

    /* Dimensions of DP array */
    size_t N = strlen(s1); // Rows
    size_t M = strlen(s2); // Columns

    /* DP array */
    int* DP_array;
    DP_array = malloc(((N+1)*(M+1))*sizeof(int));

    if(DP_array == NULL){
        exit_nomem();
    }

    // Lots of variables used throughout
    size_t element = 0;
    size_t last_col = 0;
    size_t last_row = 0;
    size_t diagnoal_above_element;
    size_t above_element;
    size_t beside_element;
    int score_diag;
    int score_abov;
    int score_besi;
    int double_include_score;

    for(size_t i = 0; i<=N; i++){

        for(size_t j = 0; j<=M; j++){

            if(i == 0 && j == 0){
                // Initial score = 0

                DP_array[0] = 0;

            } else if(i == 0){
                // Only gap penalities at this point

                last_col = element-1;
                DP_array[element] = DP_array[last_col] + GAP_SCORE;

            } else if(j == 0){
                // Only gap penalties at this point

                last_row = element-(M+1);
                DP_array[element] = DP_array[last_row] + GAP_SCORE;

            } else{

                // Three possible paths:
                diagnoal_above_element = element - (M+2);
                above_element = element - (M+2);
                beside_element = element - 1;

                // Is inclusion of both a match?
                if(strcmp(*s1, *s2) == 0){

                    double_include_score = MATCH_SCORE;

                }else{

                    double_include_score = MISMATCH_SCORE;

                }

                // Calc scores for each
                score_diag = DP_array[diagnoal_above_element] + double_include_score;
                score_abov = DP_array[above_element] + GAP_SCORE;
                score_besi = DP_array[beside_element] + GAP_SCORE;

                DP_array[element] = find_max3(score_diag, score_abov, score_besi);



            }


            element++;

        }

    }


}

int find_max3(int a, int b, int c){

    if(a >= b){

        if(b >= c){

            return b;

        }else{
            return c;
        }

    }else if(a >= c){

        return a;

    }else{

        return c;
    }

}