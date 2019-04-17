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

    for (i = 0; i < size; i++) marked[i] = 0;
    if (!id) index = 0;
    prime = 3;
    do {
		//求每个进程当前prime对应的first（求的first其实是对应的index）
        if (prime * prime > low_value){
		    first = (prime * prime - low_value)/2;
		}
        else {
            if (low_value % prime == 0) first = 0;
			//当不可以出尽时候，然后算出来的first的真是数偶数和奇数的情况
            else first = ((low_value%prime)%2==0)?((low_value+2*prime-low_value%prime-low_value)/2):((low_value+prime-low_value%prime-low_value)/2);
        }
		
        for (i = first; i < size; i = i + prime) marked[i] = 1;
		
		printf("prime:%ld,id:%ld,first:%ld",prime,id,first);
		
		//计算下一个prime
        if (!id) {
            while (marked[++index]);
            prime =2*(index+1) + 1;
        }
        if (p > 1) MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } while (prime * prime <= m);
    count = 0;
    for (i = 0; i < size; i++)
        if (!marked[i]) count++;
    if (p > 1)
        MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);

    /* Stop the timer */

    elapsed_time += MPI_Wtime();


    /* Print the results */
    printf("id:%d\n",id);
    printf("low_value:%lld,high_value:%lld\n",low_value,high_value);
    printf("my id is %d,the number of su shu is %ld\n",id,count);

    // 统计中没有算上2，需要在global——count中加上1
    if (!id) {
        printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count+1, elapsed_time, p);
    }
    MPI_Finalize();
    return 0;

}
