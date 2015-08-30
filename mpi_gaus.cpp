#include"mpi_gaus.h"
#include<stdio.h>

//#define FULL_DEBUG
#ifdef FULL_DEBUG
#define DEBUG
#endif


void makeInstance(float **A,float* B,int j,int size)
{
		int k=0;
		float keyItem=A[j][j];

		for( k=0;k<size;k++)
		{ 
				A[j][k]/=keyItem;
		}
		B[j]/=keyItem;

}


void gaus(float ** A,float* B,float* X,int size)
{
		// Forward move
		// Choose special item ( diagonal item -x)
		// x 0 0 0
		// 0 x 0 0
		// 0 0 x 0
		// 0 0 0 x
		int i=0,j=0,k=0;
		float specialItem=0;
		float koeff=0;
		for( j = 0; j < size; j++ )
		{
				specialItem=A[j][j];
				//making specialItem equal to 1
				makeInstance(A,B,j,size);





				// Going downwards from specialItem 
				for( i=j+1;i<size;i++)
				{
						koeff=A[i][j];

						// Substracting all items in the row


						for( k=0;k<size;k++)
						{
								A[i][k]=A[i][k]-(koeff*A[j][k]);
						}
						B[i]-=koeff*B[j];


				}

		}
		//Backward move
		//1 0 0 0 | b1
		//0 1 0 0 | b2
		//0 0 1 0 | b3
		//0 0 0 1 | b4

		
				for( i=size-1;i>0; --i)
				{
						for( j=0;j<i;++j)
						{
								B[j]-=A[j][i]*B[i];
								A[j][i]=0;

						}
				}
				for( i=0;i<size;++i)
				{
						X[i]=B[i];
				}
				

		/*
				X[size-1]=B[size-1];
				float tmp=0;
				for(int i=size-2;i>=0;i--)
				{

						for(int j=size-1;j>i;j--)
						{
								tmp=tmp+A[i][j]*X[j];
						}
						X[i]=B[i]-tmp;
						tmp=0;
				}
				*/
}


void printMatrix(float** matrix,int size)
{

		int i=0,j=0;
	for( i=0;i<size;i++)
	{
			for( j=0;j<size;j++)
			{
					printf("%f\t",matrix[i][j]);
			//		std::cout<<matrix[i][j]<<"\t";
			}
			puts("\n");
	}	
}

void printEqualsVector(float* vec,int size)
{
		int i=0;
		for( i=0;i<size;i++)
					printf("%f\t",vec[i]);
				//std::cout<<vec[i]<<"\t";
}

void printResult(float* vec,int size)
{
		int i=0;
		for( i=0;i<size;i++)
		{
				printf("X = %d %f",i,vec[i]);
				//std::cout<<"X"<<i<<" = "<<vec[i]<<std::endl;
		}
}
