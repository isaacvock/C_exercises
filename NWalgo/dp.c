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

int* get_score_matrix(const char* s1, const char* s2);
int find_max3(int a, int b, int c);
size_t which_max3(int a, int b, int c);
alignment* get_best_alignment(int *S, const char* s1, const char* s2);

int main(int argc, char *argv[]){

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

    // Best scores
    int* scores = get_score_matrix(seq1, seq2);

    // Best alignment

    alignment* alignment = get_best_alignment(scores, seq1, seq2);

    // Print alignment 
    for(const char *p = alignment->align1; *p; ++p){
        putchar(*p);
    }
    putchar('\n');
    for(const char *p = alignment->align2; *p; ++p){
        putchar(*p);
    }
    putchar('\n');


    // Free memory
    free(alignment->align1);
    free(alignment->align2);
    free(scores);

}


/* Walk through score matrix to get alignment */
alignment* get_best_alignment(int *S, const char* s1, const char* s2){

    /* Dimensions of DP array */
    size_t N = strlen(s1); // Rows
    size_t M = strlen(s2); // Columns

    size_t final_index = ((N+1)*(M+1)) - 1;
    size_t current_index = final_index;

    /* Alignment object */
    alignment* alignment = malloc(sizeof(alignment));
    alignment->align1 = malloc(N+M+1); // Longest possible alignment string
    alignment->align2 = malloc(N+M+1);

    /* Variables used throughout */

    // Squares we can travel to
    size_t diagnoal_above_element;
    size_t above_element;
    size_t beside_element;
    
    // Scores in squares we can travel to
    int score_diag;
    int score_abov;
    int score_besi;

    // Where we are in the matrix
    size_t current_row = N;
    size_t current_col = M;

    size_t which_score; // Which square do we go to next (highest score)?
    size_t len = 0; // Current alignment length

    while(current_index != 0){

        if(current_row == 0){

            // Only going to the left at this point
            current_index = current_index - 1;
            alignment->align1[len] = '-';
            alignment->align2[len] = s2[current_col - 1];
            current_col--;


        } else if(current_col == 0){

            // Only going up at this point
            current_index = current_index - (M+2);
            alignment->align1[len] = s1[current_row - 1];
            alignment->align2[len] = '-';
            current_row--;


        } else{

            // Three possible paths:
            diagnoal_above_element = current_index - (M+2);
            above_element = current_index - (M+1);
            beside_element = current_index - 1;


            // Calc scores for each
            score_diag = S[diagnoal_above_element];
            score_abov = S[above_element];
            score_besi = S[beside_element];

            which_score = which_max3(score_diag, score_abov, score_besi);

            if(which_score == 0){

                // Go to diagonal square
                // Include next nt of both strings; seq1[current_row - 1], seq2[current_col - 1]
                alignment->align1[len] = s1[current_row - 1];
                alignment->align2[len] = s2[current_col - 1];
                current_row--;
                current_col--;
                current_index = diagnoal_above_element;


            }else if(which_score == 1){
                // Go to square above
                // Include next nt of seq1, and gap
                alignment->align1[len] = s1[current_row - 1];
                alignment->align2[len] = '-';
                current_row--;
                current_index = above_element;

            }else{

                // Go to square to the left
                // Include next nt of seq2, and gap
                alignment->align1[len] = '-';
                alignment->align2[len] = s2[current_col - 1];
                current_col--;
                current_index = beside_element;

            }


        }

        len++;


    }

    // NULL-terminate
    alignment->align1[len] = '\0';
    alignment->align2[len] = '\0';

    return alignment;


}

/* Run N-W algorithm*/
int* get_score_matrix(const char* s1, const char* s2){

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
                above_element = element - (M+1);
                beside_element = element - 1;

                // Is inclusion of both a match?
                if(s1[i-1] == s2[j-1]){

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

    return DP_array;

}

int find_max3(int a, int b, int c){

    if(a >= b){

        if(a >= c){

            return a;

        }else{
            return c;
        }

    }else if(b >= c){

        return b;

    }else{

        return c;
    }

}


size_t which_max3(int a, int b, int c){

    if(a >= b){

        if(a >= c){

            return 0;

        }else{
            return 2;
        }

    }else if(b >= c){

        return 1;

    }else{

        return 2;
    }

}