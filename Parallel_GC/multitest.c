#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>


int main(int argc, char* argv[]){

    clock_t start_time = clock();

    int array[1000];

    #pragma omp parallel for
    for(int i = 0; i<1000; i++){

        array[i] = rand();

    }

    clock_t end_time = clock();

    double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("Time spent running: %f\n", time_spent);

}