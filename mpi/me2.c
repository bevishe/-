#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a, b)  ((a)<(b)?(a):(b))

int main(int argc, char *argv[]) {
    unsigned long int count;        /* Local prime count */
    double elapsed_time; /* Parallel execution time */
    unsigned long int first;        /* Index of first multiple */
    unsigned long int global_count = 0; /* Global prime count */
    unsigned long long int high_value;   /* Highest value on this proc */
    unsigned long int i;
    int id;           /* Process ID number */
    unsigned long int index;        /* Index of current prime */
    unsigned long long int low_value;    /* Lowest value on this proc */
    char *marked;       /* Portion of 2,...,'n' */
    unsigned long long int n;            /* Sieving from 2, ..., 'n' */
    int p;            /* Number of processes */
    unsigned long int proc0_size;   /* Size of proc 0's subarray */
    unsigned long int prime;        /* Current prime */
    unsigned long int size;         /* Elements in 'marked' */
    int mod;
    MPI_Init(&argc, &argv);
	
	char* primes;
	int prime_multiple;
	

    /* Start the timer */

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    //make sure have the number of n
    if (argc != 2) {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    //know the number we have to mark
    n = atoll(argv[1]);
    int m = n;
    n = (n%2==0)?(n/2-1):((n-1)/2);

    //给你不同进程id来分配不同的low_value和high_value
    if(n%p != 0){
	mod = n%p;
	if(id<=mod - 1){
	    low_value=3+2*id*(n/p+1);
	    high_value=3+2*id*(n/p+1)+n/p*2;
	}else{
	    low_value=3+2*(n/p+1)*mod+2*(id-mod)*(n/p);
	    high_value=3+2*(n/p+1)*mod+2*(id-mod)*(n/p)+(n/p-1)*2;
	}
    }else{
        low_value=3+2*id*(n/p);
	high_value=3+2*id*(n/p)+2*(n/p-1);
    }
    // the size of verey process's list we have to malloc
    size = (high_value-low_value)/2 + 1;

    proc0_size = (high_value-low_value)/2;
	
    //判断根号n是否在零号进程标注的marked中
    if ((2 + proc0_size) < (int) sqrt((double) n)){
        if (!id) printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }
	
    /* Allocate this process's share of the array. */
    marked = (char *) malloc(size);

    if (marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }

    // compute primes from 2 to sqrt(n);
    int sqrt_n = sqrt(m);
    primes = (char*)malloc(sqrt_n+1);
	int prime_i = 2;
	int prime_first;
	int prime_first_;
	do{
		prime_first = prime_i*prime_i;
		for(prime_first_ = prime_first;prime_first_<=sqrt_n;prime_first_ = prime_first_ + prime_i){
			primes[prime_first_] = 1;
		}
		while(primes[++prime_i]==1);
	}while(prime_i*prime_i<=sqrt_n);


    for (i = 0; i < size; i++) marked[i] = 0;
    if (!id) index = 0;
    prime = 3;
	for (prime = 3; prime <= sqrt_n; prime++){
		//求每个进程当前prime对应的first（求的first其实是对应的index）
		if (primes[prime] == 1)
            continue;
        if (prime * prime > low_value)
            first = (prime * prime - low_value)/2;
        else {
            if (!(low_value % prime)) first = 0;
            else first = ((low_value%prime)%2==0)?((low_value+2*prime-low_value%prime-low_value)/2):((low_value+prime-low_value%prime-low_value)/2);
        }
        for (i = first; i < size; i = i + prime  ) marked[i] = 1;
		
		//给所有的primes数组进行marked标记
	}
    count = 0;
    for (i = 0; i < size; i++)
        if (!marked[i]) count++;
    if (p > 1)
        MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);

    /* Stop the timer */

    elapsed_time += MPI_Wtime();


    

    // 统计中没有算上2，需要在global——count中加上1
    if (!id) {
        printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count+1, elapsed_time, p);
    } 
    MPI_Finalize();
    return 0;

}
