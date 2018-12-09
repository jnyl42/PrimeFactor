/* Factor Semiprimes
 * By Jason Nyland and Anthony Volkov
 * December 9, 2018
 *
 * This program uses OpenMP to implement multi-threading for various aspects involving
 * factoring the product of two prime numbers.
 *
 * We first generate a list of primes using the Sieve of Eratosthenes algorithm, and then
 * check each result by searching for a prime where 'input % prime == 0'.
 *
 */

#include <stdio.h>
#include <math.h>
#include <tgmath.h>
#include <stdbool.h>
#include <stdlib.h>
#include <omp.h>


unsigned long long product;  //user input goes here
unsigned long long primeArraySize;  //size of array returned by generatePrimes (needed later)

//
// generates prime numbers up to max; returns a pointer to an array of prime numbers
//
unsigned long long * generatePrimes(long double max) {

    double start, end;  //timekeeping values
    unsigned long long sieveSize = ((unsigned long long) max) + 1;

    printf("Initializing sieve array... ");
    start = omp_get_wtime();
    bool * sieveArray = (bool *)malloc(sieveSize*sizeof(bool));
    end = omp_get_wtime();
    printf("Done (%lf s)\n", end-start);

    printf("[OMP] Setting initial values... ");
    start = omp_get_wtime();
    //mark 0 and 1 as not prime
    sieveArray[0] = false;
    sieveArray[1] = false;
    #pragma omp parallel  //using OMP to fill the array for speed (not really effective)
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();

        for (int i = tid+2; i < sieveSize; i+=nthreads) {
            sieveArray[i] = true;
        }
    }
    end = omp_get_wtime();
    printf("Done (%lf s)\n", end-start);

    printf("[OMP] Sieving primes... ");
    start = omp_get_wtime();
    //for every index < sqrt(sieveSize)
    for (unsigned long long i = 2; (i < sqrtl(sieveSize)); i++) {
        //if value is prime [true]
        if (sieveArray[i] == true) {
            //use omp to mark factors as non-prime [false] for speed (not really effective)
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                int nthreads = omp_get_num_threads();
                for (unsigned long long j = (unsigned long long) tid + 2; (i * j) < sieveSize; j+=nthreads) {
                    sieveArray[i * j] = false;
                }
            }
        }
    }
    end = omp_get_wtime();
    printf("Done (%lf s)\n", end-start);

    printf("Counting primes... ");
    start = omp_get_wtime();
    //count primes for output array size
    unsigned long long count = 0;
    for (unsigned long long i = 2; i < sieveSize; i++) {
        if (sieveArray[i] == true) {
            count++;
        }
    }
    end = omp_get_wtime();
    printf("Done (%lf s)\n", end-start);

    printf("Creating array of %lld primes... ", count);
    start = omp_get_wtime();
    unsigned long long *output = (unsigned long long *)malloc(count*sizeof(unsigned long long));
    primeArraySize = count;
    count = 0;
    //add primes to output array
    for (unsigned long long i = 2; i < sieveSize; i++) {
        if (sieveArray[i] == true) {
            output[count] = i;
            count++;
        }
    }
    end = omp_get_wtime();
    printf("Done (%lf s)\n", end-start);

    free(sieveArray);
    return output;
}

int main(int argc, char** argv) {
    double start = omp_get_wtime();
    //get product from terminal arguments
    if (argc > 1) {
        char * pEnd;
        product = strtoull(argv[1],&pEnd,10);
    } else {
        printf("No args found; exiting.\n");
        return 1;
    }
    printf("----------------------------------------------------------\n");
    printf("  Number to Factor:   %lld \n", product);
    printf("----------------------------------------------------------\n");

    //generate a list of primes up to the max lower prime
    unsigned long long * primeList = generatePrimes(sqrtl(product));

    //use OMP to check each prime against product for speed (very effective)
    printf("[OMP] Attempting to factor... ");
    double fStart = omp_get_wtime();
    bool found = false;
    unsigned long long i;
    #pragma omp parallel shared(found) private (i)
    {
        int nthreads = omp_get_num_threads();
        int tid = omp_get_thread_num();

        for (i = (unsigned long long) tid; i < primeArraySize && !found; i+=nthreads) {
            if ( fmodl(product, primeList[i]) == 0.0 ) {
                double end = omp_get_wtime();
                printf("Done (%lf s)\n", end-fStart);
                printf("----------------------------------------------------------\n");
                printf("  Primes found:   %lld, %lld \n", primeList[i], (product/primeList[i]));
                printf("----------------------------------------------------------\n");
                printf("total run time: %lf", end-start);
                found = true;
            }
        }
    }

    if (found) {
        free(primeList);
        return 0;
    }
    else {
        printf("Prime factors not found.\n");
        free(primeList);
        return 1;
    }
}