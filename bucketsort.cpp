#include <math.h>
#include <time.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

using namespace std;

void bucketSort(float* arr, int& n, int bucket_size) {
  // Crear buckets      
  std::vector<float> bucket[bucket_size];
  int i, j, bucket_index, index = 0;

  // asignar elementos a los buckets
  for (i = 0; i < n; ++i) {
    	bucket_index = bucket_size*arr[i]/1000;
      	bucket[bucket_index].push_back(arr[i]);
  }

  // ordenar buckets
  for (j = 0; j < bucket_size; ++j)
    sort(bucket[j].begin(), bucket[j].end());

  // Concatenar buckets en arr[]
  for (i = 0; i < bucket_size; i++) {
    for (j = 0; j < bucket[i].size(); j++)
      arr[index++] = bucket[i][j];
  }
}

int main(int argc, char *argv[]) {
  int i, n = pow(2, 19);     
  float* randArray;
  
  srand((int)time(0));
    
  randArray = new float[n];

  for(int i = 0; i < n; ++i)
    randArray[i]=(float)rand()/(float)(RAND_MAX/999.0);

  // ordenar array en buckets
  bucketSort(randArray, n, 8);
  
  delete[] randArray;
}
