#include <math.h>
#include <time.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <mpi.h>
using namespace std;

int nrank, nsize;

void bucketSort(float* arr, int& n, int bucket_size) {
  // Crear buckets      
  vector<float> bucket[bucket_size];
  int i, j, bucket_index, index = 0, size_b[32]={}; 

  // asignar elementos a los buckets
  for (i = (nrank * (n/nsize)); i < ((n/nsize) * (nrank + 1)); ++i) {
    bucket_index = bucket_size*arr[i]/1000;
    bucket[bucket_index].push_back(arr[i]);
  }

  // ordenar buckets
  sort(bucket[nrank].begin(), bucket[nrank].end());

  // Obtener valores previos al bucket
  for (i = 0; i < nsize; i++)
  {
    size_b[i]= 0;
    for (j = 0; j < i; j++)
    {
      size_b[i] += bucket[j].size();
    }
  }

  index = size_b[nrank];

  // Concatenar buckets en arr[]
  for (j = 0; j < bucket[nrank].size(); j++)
  {
    arr[index] = bucket[nrank][j];

    if (nrank != 0)
    {
      MPI_Send(&arr[index], 1, MPI_FLOAT, 0, nrank, MPI_COMM_WORLD);
    }
    if (nrank == 0) 
    {
      for (int z = 1; z < nsize; z++)
      {
        MPI_Recv(&arr[size_b[z]], 1, MPI_FLOAT, z, z, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
    index++;
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

double t0=0.0, tf=0.0;

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &nrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nsize);
  
  int i, n = pow(2, 19);     
  float* randArray;
  
  srand((int)time(0));
    
  randArray = new float[n];

  for(int i = 0; i < n; ++i)
  {
    randArray[i]=(float)rand()/(float)(RAND_MAX/999.0);
  }
  
  /*
  for (i = 0; i < n; ++i) 
    	printf("%1.2f \n", randArray[i]);
	printf("\n");
  */

  t0 = MPI_Wtime(); //Tiempo Inicial

  // ordenar array en buckets
  bucketSort(randArray, n, nsize);
  
  /*
  printf("ORDENADO\n");
  for (i = 0; i < n; ++i) 
    printf("%1.2f \n", randArray[i]);
  */

  delete[] randArray;
  
  if (nrank==0) 
  {
    tf=MPI_Wtime(); //Tiempo Final
    printf("Tiempo: %1.6f s\n", (tf-t0));
  }
  MPI_Finalize();
}
