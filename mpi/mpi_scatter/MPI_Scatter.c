
#include <stdlib.h>
#include "mpi.h"
#include <stdio.h>
#include <string.h>

#define N 8

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Read_vec
 *  Description:  
 * =====================================================================================
 */
	void
Read_vec ( int *local_a,
		   int local_n,
		   int n,
		   int myrank,
		   MPI_Comm comm)
{
	int *a = NULL;
	int i = 0;
	if (myrank == 0){
		a = malloc(n * sizeof(int));
		printf("Enter the number\n");
		for (i=0;i<n;i++)
			scanf("%d", &a[i]);

		MPI_Scatter(a,  local_n, MPI_INT, local_a, local_n, MPI_INT, 0, comm);
		free(a);
	}
	else{	
		MPI_Scatter(a,  local_n, MPI_INT, local_a, local_n, MPI_INT, 0, comm);

		/*-----------------------------------------------------------------------------
		 *  发现1：这里发现如果对于进程号不为0的不做处理的话那么其他进程在等待输入
		 *  的时候就已经开始执行导致无法接收到正确数据。
		 *  发现2：在等待输入过程中，CPU的使用率很高，这条语句是不是用于接收?
		 *-----------------------------------------------------------------------------*/
	}
}		/* -----  end of function Read_vec  ----- */
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
	int
main ( int argc, char *argv[] )
{
	MPI_Init(NULL, NULL);
	int myrank = 0;
	int size = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int *local_a = malloc( (N/size)* sizeof(int));
	int local_n = N/size;
	int j = 0;			
	Read_vec(local_a, local_n, N, myrank, MPI_COMM_WORLD);
	for(j=0;j<local_n;j++)
		printf("i am %d, i receive %d\n", myrank, local_a[j]);
	free(local_a);
	MPI_Finalize();
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
