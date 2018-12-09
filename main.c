#include <stdio.h>
#include <math.h>
#include <tgmath.h>
#include <stdbool.h>
#include <stdlib.h>
#include <omp.h>

unsigned long long product = 71;
//double primes[8] = {2,3,5,7,11,13,17,19};

unsigned long long * generatePrimes(long double max) {
    double start;
    double end;
    unsigned long long cap = ((unsigned long long) max) + 1;
    printf("Cap of primes to generate: %lld \n", cap);
    unsigned long long count = 0;
    printf("Initializing sieve array... ");
    start = omp_get_wtime();
    int * array = (int *)malloc(cap*sizeof(int));
    end = omp_get_wtime();
    printf("Done (%lf s)\n", end-start);

//mark 0 and 1 as not prime
    array[0] = 1;
    array[1] = 1;

    printf("Sieving primes... ");
    start = omp_get_wtime();
    //for every array index that is still marked 1 for prime
    for (unsigned long long i = 2; (i < cap) && (array[i] == 0); i++) {
        //multiply the index number by consecutive integers and mark each product as not prime
        for (int j = 2; (i * j) < cap; j++) {
            array[i * j] = 1;
        }
    }
    end = omp_get_wtime();
    printf("Done (%lf s)\n", end-start);
    printf("Counting primes... ");
    start = omp_get_wtime();
    //count primes for output array size
    for (unsigned long long i = 2; i < cap; i++) {
        if (array[i] == 0) {
            //printf("%d\n", i);
            count++;
        }
    }
    end = omp_get_wtime();
    printf("Done (%lf s)\n", end-start);
    printf("Creating array of %lld primes... ", count);
    start = omp_get_wtime();
    unsigned long long *output = (unsigned long long *)malloc(count*sizeof(unsigned long long));
    count = 0;
    //add primes to output array
    for (unsigned long long i = 2; i < cap; i++) {
        if (array[i] == 0) {
            output[count] = i;
            count++;
        }
    }
    end = omp_get_wtime();
    printf("Done (%lf s)\n", end-start);
    printf("Free sieve array and return list of primes array\n");
    free(array);
    return output;
}

int main(int argc, char** argv) {
    double start = omp_get_wtime();
    //get product from terminal arguments
    if (argc > 1) {
        char * pEnd;
        product = strtoull(argv[1],&pEnd,10);
    } else {
        printf("No args found; using default value.\n");
    }
    printf("----------------------------------------------------------\n");
    printf("User Input = %lld \n", product);
    long double maxLower = sqrtl(product);
    printf("Maximum Lower Prime = %Lf \n", maxLower);
    printf("----------------------------------------------------------\n");
    //generate a list of primes up to the max lower prime
    unsigned long long * primeList = generatePrimes(maxLower);
    //check each prime against product
    printf("Attempting to factor... ");
    double fStart = omp_get_wtime();
    for (unsigned long long i = 0; primeList[i] != 0; i++) {
        if ( fmodl(product, primeList[i]) == 0.0 ) {
            double end = omp_get_wtime();
            printf("Done (%lf s)\n", end-fStart);
            printf(">>>>> Primes found: %lld, %lld <<<<<\n", primeList[i], (product/primeList[i]));
            printf("total run time: %lf", end-start);
            return 0;
        }
    }
    printf("Prime factors not found.\n");
    return 1;
}