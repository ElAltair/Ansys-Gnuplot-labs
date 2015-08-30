#include<iostream>
#include<iomanip>
#include<fstream>
#include"mpi_gaus.h"
#include<string>
#include<sstream>
#include<unistd.h>
#include<math.h>
#include<vector>

//#define DEBUG


float deltaX = 0.0;
float deltaT = 0.0;
float deltaY = 0.0;

struct Node
{
		int id;
		float x;
		float y;
		Node(float _id=0,float _x=0,float _y=0):id(_id),x(_x),y(_y){};

}** matrix;

float **TemperatureMatrix;
float StartTemp=20;

struct TempConditions
{
		float topTemp;
		float leftTemp;
		float rightTemp;
		float bottomTemp;
		bool secondCondTopTemp;
		bool secondCondLeftTemp;
		bool secondCondRightTemp;
		bool secondCondBottomTemp;
		TempConditions(float top=0,float left=0,float right=0,float bottom=0,bool btop=false,bool bleft=false, bool bright=false, bool bbottom=false)
				:topTemp(top),leftTemp(left),rightTemp(right),bottomTemp(bottom),secondCondTopTemp(btop),secondCondLeftTemp(bleft),
				secondCondRightTemp(bright),secondCondBottomTemp(bbottom){};
};

TempConditions* startCond;
	//	float** A;                        // Big matrix of all equations for each node
//		float* B;                         // Vector of values
//		float* Temperature;               // Vector of unknown values

void createMatrixes(float** &A,float* &B,float* &Temperature,int nodeWidth,int nodeHeight)
{
		int allNodes=nodeWidth*nodeHeight;
		A= new float*[allNodes];

		B=new float[allNodes];
		Temperature=new float[allNodes];

		matrix= new Node*[nodeHeight];
		for(int i=0;i<nodeHeight;++i)
				matrix[i]=new Node[nodeWidth];


		for(int i=0;i<allNodes;++i)
				A[i]=new float[allNodes];

		TemperatureMatrix =  new float*[nodeHeight];
		for(int i=0;i<nodeWidth;++i)
		{
				TemperatureMatrix[i] = new float [nodeWidth];
		}
		for(int i=0;i<nodeHeight;++i)
		{
				for(int j=0;j<nodeWidth;++j)
				{
						TemperatureMatrix[i][j]=StartTemp;
				}
		}


}
/*
void makeInstance(float **A,float* B,int j,int nodeWidth)
{
		float keyItem=A[j][j];

		for(int k=0;k<nodeWidth;k++)
		{ 
				A[j][k]/=keyItem;
		}
		B[j]/=keyItem;

}
*/

void forwardGaus(float ** A,float* B,int nodeWidth,int nodeHeight)
{
		// Forward move
		// Choose special item ( diagonal item -x)
		// x 0 0 0
		// 0 x 0 0
		// 0 0 x 0
		// 0 0 0 x
		float specialItem=0;
		float koeff=0;
		for(int j = 0; j < nodeHeight; j++ )
		{
				specialItem=A[j][j];
				//making specialItem equal to 1
				makeInstance(A,B,j,nodeWidth);





				// Going downwards from specialItem 
				for(int i=j+1;i<nodeHeight;i++)
				{
						koeff=A[i][j];

						for(int k=0;k<nodeWidth;k++)
						{
								A[i][k]=A[i][k]-(koeff*A[j][k]);
						}
						B[i]-=koeff*B[j];


				}

		}
}

void backwardGaus(float** A,float* B,float* X,int size)
{
		//Backward move
		//1 0 0 0 | b1
		//0 1 0 0 | b2
		//0 0 1 0 | b3
		//0 0 0 1 | b4


		for(int i=size-1;i>0; --i)
		{
				for(int j=0;j<i;++j)
				{
						B[j]-=A[j][i]*B[i];
						A[j][i]=0;

				}
		}
		for(int i=0;i<size;++i)
		{
				X[i]=B[i];
		}



}

void printToFile(std::string route,float** mat,float* vec,int nodeWidth,int nodeHeight,float count)
{
		std::ofstream file;
		std::stringstream convert;
		convert<<route<<count<<".data";
		std::string temp= convert.str();

		const char* name = temp.c_str();

		file.open(name);
		int allNodes=nodeHeight*nodeWidth;

		for(int i=0;i<nodeHeight;++i)
		{
				for(int j=0;j<nodeWidth;++j)
				{
						file<<matrix[i][j].x<<"\t"<<nodeHeight - matrix[i][j].y<<"\t"<<TemperatureMatrix[i][j]<<"\n";
				}
				file<<std::endl;
		}
		file.close();
}

void createGnuplotFile(int count,float stopTime=0.1)
{
		std::ofstream file;
		file.open("gnuplot.gp");
		file<<"set view map"<<std::endl;
		for(int i=0;i<count;++i)
		{
				file<<"splot \"out/node"<<i<<".data\""<<" with pm3d"<<std::endl;
				file<<"pause "<<stopTime<<std::endl;
		}
		file<<"pause -1";
}


void printTempMatrix(int nodeWidth,int nodeHeight)
{
		for(int i=0;i<nodeHeight;++i)
		{
				for(int j=0;j<nodeWidth;++j)
				{
						std::cout<<std::setw(8)<<TemperatureMatrix[i][j]<<" ";
				}
				std::cout<<std::endl;
		}

}



void printNodeMatrix(int nodeWidth,int nodeHeight)
{

		for(int i=0;i<nodeHeight;++i)
		{
				for(int j=0;j<nodeWidth;++j)
						std::cout<<matrix[i][j].id<<"\t";

				std::cout<<std::endl;
		}
}

/*
   void printMatrix(int** matrix,int width,int height)
   {
   for(int i=0;i<height;++i)
   nodeWidth
   for(int j=0;j<width;++j)
   std::cout<<matrix[i][j]<<"\t";

   std::cout<<std::endl;
   }


   }
   */
void printVector(float* vec,int height)
{
		for(int i=0;i<height;++i)
				std::cout<<vec[i]<<std::endl;
}

void printFormatMatrixes(float** A,float *B,int width,int height)
{
		std::streamsize ss = std::cout.precision();
		std::cout<<" ";
		for(int i=0;i<width;++i)
				std::cout<<std::setw(3)<<i<<" ";
		std::cout<<std::endl;
		double floorNum=0.0;
		for(int i=0;i<height;++i)
		{
				std::cout<<"|";
				for(int j=0;j<width;++j)
				{

						/*
						   if((floorNum=floor(A[i][j]))==A[i][j])
						   std::cout<<std::setprecision(2)<<std::setw(3)<<A[i][j]<<" ";
						   else
						   {
						   std::cout<<std::setiosflags(std::ios::fixed)<<std::setprecision(2)<<std::setw(3)<<A[i][j]<<" ";
						   std::cout<<resetiosflags(std::ios::fixed);
						   }
						   */
						std::cout<<std::setw(3)<<A[i][j]<<" ";

				}
				std::cout<<" "<<"|"<<" "<<"|"<<" "<<"T"<<std::setw(3)<<i<<" "<<"|"<<" "<<"|"<<" "<<B[i]<<" "<<"|";
				std::cout<<std::endl;
		}
		std::cout.precision(ss);
}


void fillBlockNodeMatrix(int nodeWidth,int nodeHeight,int blockWidth,int blockCount)
{
		int boundaryI=nodeWidth+1;
		int boundaryId=nodeHeight*blockWidth*blockCount;
		int innaryId=0,prevI=0;
		int i=0,j=0,k=0;

		for(i=0;i<blockCount;++i)
		{
				prevI=i*(blockWidth+1);
				for(j=0;i<nodeHeight;++i)
				{
						for(k=0;k<blockWidth;++k)
						{
								matrix[j][k+prevI].id=innaryId++;
						}
						boundaryI=k;

				}
				for(k=0;k<nodeHeight;++k)
				{
						matrix[k][boundaryI]=boundaryId++;
				}
		}



}

void fillNodeMatrix(int nodeWidth,int nodeHeight)
{

		int count=0;
		int allNodes=nodeHeight*nodeWidth;
		int i=0,j=0;
		while(count<allNodes)
		{
				if(j!=0 && count%nodeWidth==0)
				{
						i++;
						j=0;
						continue;
				}
				matrix[i][j].id=count;
				matrix[i][j].x=j*deltaX;
				matrix[i][j].y=i*deltaY;
				j++;
				count++;
		}

}

void fillEquationResult(float* B,int nodeWidth,int nodeHeight)
{
		for(int i=0;i<nodeHeight;++i)
				for(int j=0;j<nodeWidth;++j)
				{
						if(matrix[i][j].y==0)
						{
								if(startCond->secondCondTopTemp)
								{
										B[matrix[i][j].id]=startCond->topTemp*deltaY;
								}
								else
								{
										B[matrix[i][j].id]=startCond->topTemp;
								}

						}
						else if(matrix[i][j].x==0)
						{
								if(startCond->secondCondLeftTemp)
								{
										B[matrix[i][j].id]=startCond->leftTemp*deltaX;
								}
								else
										B[matrix[i][j].id]=startCond->leftTemp;
						}

						else if(matrix[i][j].x==deltaX*(nodeWidth-1))
						{
								if(startCond->secondCondRightTemp)
								{
										B[matrix[i][j].id]=startCond->rightTemp*deltaX;
								}
								else
								{
										B[matrix[i][j].id]=startCond->rightTemp;
								}

						}
						else if(matrix[i][j].y==deltaY*(nodeHeight-1))
						{
								if(startCond->secondCondBottomTemp)
								{
										B[matrix[i][j].id]=startCond->bottomTemp*deltaY;
								}
								else
								{
										B[matrix[i][j].id]=startCond->bottomTemp;
								}
						}
						else
						{
								B[matrix[i][j].id]=TemperatureMatrix[i][j]*deltaX*deltaX*deltaY*deltaY;
						}
				}

}


void fillEquationMatrix(float** A,int nodeWidth,int nodeHeight)
{

		int allNodes= nodeWidth*nodeHeight;
		int count = 0;
		int i=0,j=0,k=0;
		int equalPos=0;
		if(!startCond)
		{
				startCond = new TempConditions(200,0,200,0,0,1,0,1);
		}




		while(count<allNodes)
		{
				equalPos=matrix[i][j].id;
				if(j!=0 && j%nodeWidth==0)
				{
						i++;
						j=0;
						continue;

				}
				if(matrix[i][j].y==0)
				{
						if(startCond->secondCondTopTemp)
						{
								//	B[count]=startCond->topTemp*deltaY;
								A[count][equalPos]=1;
								A[count][equalPos+nodeHeight]=-1;
						}
						else
						{
								//   B[count]=startCond->topTemp;
								A[count][equalPos]=1;
						}
						count++;
						j++;
						equalPos++;
						continue;
				}
				else if(matrix[i][j].x==0 )
				{
						if(startCond->secondCondLeftTemp)
						{
								//	B[count]=startCond->leftTemp*deltaX;
								A[count][equalPos]=1;
								A[count][equalPos+1]=-1;
						}
						else
						{
								A[count][equalPos]=1;
								//	B[count]=startCond->leftTemp;
						}
						count++;
						j++;
						equalPos++;
						continue;

				}
				else if(matrix[i][j].x==deltaX*(nodeWidth-1))
				{
						if(startCond->secondCondRightTemp)
						{
								//	B[count]=startCond->rightTemp*deltaX;
								A[count][equalPos]=1;
								A[count][equalPos-1]=-1;
						}
						else
						{
								//   B[count]=startCond->rightTemp;
								A[count][equalPos]=1;
						}
						count++;
						j++;
						equalPos++;
						continue;
				}
				else if(matrix[i][j].y==deltaY*(nodeHeight-1))
				{
						if(startCond->secondCondBottomTemp)
						{
								//	B[count]=startCond->bottomTemp*deltaY;
								A[count][equalPos]=1;
								A[count][equalPos-nodeHeight]=-1;
						}
						else
						{
								A[count][equalPos]=1;
								//	B[count]=startCond->bottomTemp;
						}
						count++;
						j++;
						equalPos++;
						continue;
				}

				else
				{
						A[count][equalPos]=(2*deltaX*deltaX*deltaT+2*deltaY*deltaY*deltaT)+(deltaX*deltaX*deltaY*deltaY);
						A[count][equalPos-1]=-deltaY*deltaY*deltaT;
						A[count][equalPos+1]=-deltaY*deltaY*deltaT;
						A[count][equalPos-nodeWidth]=-deltaX*deltaX*deltaT;
						A[count][equalPos+nodeWidth]=-deltaX*deltaX*deltaT;
						//	B[count]=TemperatureMatrix[i][j]*deltaX*deltaX*deltaY*deltaY; // TODO
						j++;
						equalPos++;
						count++;
						continue;
				}
				i++;
				j=0;
		}
}

void copyMatrixA ( float ** &from, float ** & to,int nodeWidth, int nodeHeight)
{
		for( int i = 0; i < nodeHeight*nodeWidth; ++i)
		{
				for( int j = 0; j < nodeWidth*nodeHeight; ++j)
						to[i][j] = from[i][j];
		}
}



int main(int argc,char** argv)
{

		if(argc!=7)
		{
				std::cout<<"Usage: /progname plateWidth plateHeight nodeX nodeY deltaT fullTime\n";
				return -1;
		}



		float plateWidth = atof(argv[1]);
		float plateHeight = atof(argv[2]);
		int nodeWidth = atoi(argv[3]);
		int nodeHeight = atoi(argv[4]);
		deltaT = atof(argv[5]);
		float fullTime = atof(argv[6]);

		deltaX = plateWidth / nodeWidth;
		deltaY = plateHeight / nodeHeight;
		std::cout<<" deltaX = "<< deltaX << "\ndeltaY = " << deltaY << "\ndeltaT = "<<deltaT << std::endl;
		float ** A;
		float ** tempA;
		float* B;
		float *Temperature;
		
		int allNodes = nodeHeight * nodeWidth;

		tempA = new float*[allNodes];
		for(int i = 0; i < allNodes; ++i)
		{
				tempA[i] = new float[allNodes];
		}

		
		/* Creation dynamic elements */
		createMatrixes(A,B,Temperature,nodeWidth,nodeHeight);




		/* end dynamic creatations */

		float time = 0.0;
		int count=0;
		fillNodeMatrix(nodeWidth,nodeHeight);  // filling the node matrix
		fillEquationMatrix(A,nodeWidth,nodeHeight);
		//forwardGaus(A,B,nodeWidth,nodeHeight);
		//
#ifdef DEBUG
		printFormatMatrixes(A,B,allNodes,allNodes);
#endif
		while(time<fullTime)
		{

				//printNodeMatrix(nodeHeight,nodeWidth);
				//printTempMatrix(nodeHeight,nodeWidth);

				copyMatrixA(A,tempA,nodeWidth,nodeHeight);
#ifdef DEBUG
				printFormatMatrixes(tempA,B,allNodes,allNodes);
#endif
				//fillEquationMatrix(A,nodeWidth,nodeHeight);
				fillEquationResult(B,nodeWidth,nodeHeight);
				std::cout<<"Time = "<<time<<std::endl;
				//printVector(B,allNodes);
				//printFormatMatrixes(A,B,allNodes,allNodes);
				std::cout<<std::endl;
				//forwardGaus(tempA,B,allNodes,allNodes);
				gaus(tempA,B,Temperature,allNodes);
				//backwardGaus(A,B,Temperature,allNodes);
				//printFormatMatrixes(A,B,allNodes,allNodes);
				printToFile("out/node",tempA,B,nodeWidth,nodeHeight,count);



				for(int i=0;i<nodeHeight;++i)
						for(int j=0;j<nodeWidth;++j)
						{

								TemperatureMatrix[i][j]=Temperature[matrix[i][j].id];
						}

#ifdef DEBUG
				printTempMatrix(nodeWidth,nodeHeight);
#endif
				time+=deltaT;
				count++;
		}
		std::cout<<"Count = "<<count<<std::endl;
		createGnuplotFile(count);



		if(fork()==0)
		{
				system("gnuplot gnuplot.gp");
				return 0;
		}
		else
		{
				return 0;
		}

}

