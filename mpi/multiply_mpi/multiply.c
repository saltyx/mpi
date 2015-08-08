/*
 * =====================================================================================
 *
 *       Filename:  multiply.c
 *
 *    Description:  矩阵相乘
 *
 *        Version:  1.0
 *        Created:  2015年08月08日 14时54分19秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  shiyan (), shiyan233@hotmail.com
 *   Organization:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mpi.h"

/*-----------------------------------------------------------------------------
 *  n*m矩阵,n行m列
 *-----------------------------------------------------------------------------*/
#define N 4
#define M 2
void Transpose(double *b,double *tmp ,int n, int row, int col);
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Read_vec
 *  Description:  
 * =====================================================================================
 */
	void
Read_vec (double *local_a,
		 double *local_b,
		 int local_n,
		 int n,
		 int myrank,
		 MPI_Comm comm)
{
	double *a = NULL;
	double *b = NULL;
	double *tmp = NULL;
	if (myrank == 0)
	{
		a = malloc( n*sizeof(double) );
		int i = 0;
		printf("a is \n");
		for (i=0; i<n;i++)
		{
			a[i] = i;
			printf("%lf ", a[i]);
		}
		printf("\nb is \n");	
		b = malloc( n*sizeof(double) );
		i = 0;
		for (i=0; i<n;i++)
		{
			b[i] = i ;
			printf("%lf ", b[i]);
		}
		tmp = malloc( n*sizeof(double));
		Transpose(b, tmp, n, N, M);//将b转置之后放入tmp中
		for(i=0;i<n;i++)
			printf("\ntmp is %lf\n", tmp[i]);
		MPI_Scatter(tmp, local_n, MPI_DOUBLE, local_b, local_n, MPI_DOUBLE, 0, comm);
		MPI_Scatter(a, local_n, MPI_DOUBLE, local_a, local_n, MPI_DOUBLE, 0, comm);
		free(a);
		free(b);
		free(tmp);
	}else{
		MPI_Scatter(tmp, local_n, MPI_DOUBLE, local_b, local_n, MPI_DOUBLE, 0, comm);
		MPI_Scatter(a, local_n, MPI_DOUBLE, local_a, local_n, MPI_DOUBLE, 0, comm);	
	}	
}		/* -----  end of function Read_vec  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Gather_vec
 *  Description:  
 * =====================================================================================
 */
	void
Gather_vec (
		double *local_c,
		int local_n,
		int n,
		int myrank,
		MPI_Comm comm
		)
{
	double *t = NULL;
	int i = 0;
	if (myrank == 0){
		t = malloc(n*sizeof(double));
		MPI_Gather(local_c, local_n, MPI_DOUBLE, t, local_n, MPI_DOUBLE, 0, comm);
		printf("result:\n");
		for(i=0;i<n;i++)
			printf("%lf ", t[i]);
		putchar('\n');
		free(t);
	}else{
		MPI_Gather(local_c, local_n, MPI_DOUBLE, t, local_n, MPI_DOUBLE, 0, comm);
	}
}		/* -----  end of function Gather_vec  ----- */
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  transpose
 *  Description:  空间O(n)的转置
 * =====================================================================================
 */
void Transpose(double *b,double *tmp ,int n, int row, int col)
{
	int k = 0;
	for (k=0;k<n;k++)
	{
		tmp[k] = b[ (k%row)*col + k/row ];
	}
}		/* -----  end of function transpose  ----- */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Multiply
 *  Description:  
 * =====================================================================================
 */
	void
Multiply ()
{
	MPI_Init(NULL,NULL);

	int myrank;
	int size;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int local_n = (M*N)/size;
	double *local_a = malloc(local_n*sizeof(double));
	double *local_b = malloc(local_n*sizeof(double));
	double *local_c = malloc(local_n*sizeof(double));
	Read_vec(local_a, local_b, local_n, (M*N), myrank, MPI_COMM_WORLD);
	int i = 0;
	for (i=0;i<local_n;i++)
		local_c[i] = 0.0;
	for(i =0 ;i<local_n; i++){
		local_c[i] += local_a[i]*local_b[i];
	}
	Gather_vec(local_c, local_n, M*N, myrank, MPI_COMM_WORLD);	
	MPI_Finalize();
}		/* -----  end of function Multiply  ----- */
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
	int
main ( int argc, char *argv[] )
{
	Multiply();
	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
