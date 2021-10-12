#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "../Include/umfpack.h"

using namespace std;

double f(double x, double y)
{
  return -12.0 * (pow(x - y, 2) + 1e-5) * (5.0 * pow(x - y, 2) + 1e-5);
}

double u_answer(double x, double y)
{
  return pow((x - y) * (x - y) + 1e-5, 3);
}


int main()
{
  
  int n = 20;        // i = 0..49, j = 0..49
  int bn = n - 1;
  double h = 1.0 / n;
  double sh = h * h;
  int N = n * n;  
  int aN = N + 1;
  
  int dx = 1;
  int dy = 10;
  int* G_D_top = new int[bn];
  int* G_D_right = new int[bn];
  
  
  for (int a = 0; a < n; a++)
  {
    G_D_top[a] = 100;
  }
  for (int a = 0; a < n; a++)
  {
    G_D_right[a] = 100;
  }
  

  int* row_index = new int[aN];
  int NNZ = 5 * (n - 2)*(n - 2);   //without borders
  NNZ += 4 * (n - 2);              //right border, Dirichlet
  NNZ += 4 * (n - 2);              //top border, Dirichlet
  NNZ += 3;                        //top right corner, Dirichlet
  
  NNZ += 4 * (n - 2);              
  NNZ += 4 * (n - 2);              
  NNZ += 3;
  NNZ += 3;                         
  NNZ += 3;                         
                       
  int* column_index = new int[NNZ];
  double* values = new double[NNZ];
  double* b = new double[N];
  row_index[0] = 0;
  int cur = 0;
  
  for (int k = 0; k < N; k++)
  {
    int i =  k  / n;    
    int j = k % n;
    b[k] = cos(M_PI * i * h) * cos(M_PI * j * h);
  }
  
  for (int k = 0; k < N; k++)
  {
    int i = k / n;
    int j = k % n;
    cout << "i= " << i << "j= " << j << "\n";
    
    if ( (i < bn) && (j < bn) && (i > 0) && (j > 0) )
    {
      column_index[cur] = (i - 1) * n + j;
      values[cur] = - dy / sh;
      cur++;

      values[cur] = - dy / sh;
      column_index[cur] = i * n + (j - 1);
      cur++;

      values[cur] = 2 * (dx + dy) / sh;
      column_index[cur] = i * n + j;
      cur++;
      
      values[cur] = - dx / sh;
      column_index[cur] = i * n + (j + 1);
      cur++;
            
      values[cur] = - dx / sh;
      column_index[cur] = (i + 1) * n + j;
      cur++;
      
      row_index[k + 1] = row_index[k] + 5;
    }
    if ( (i == 0) && (j >= 1) && (j < bn) )  // left border
    {
      values[cur] = - dy / sh;
      column_index[cur] = i * n + (j - 1);
      cur++;
      
      values[cur] = 2 * (dx + dy) / sh;
      column_index[cur] = i * n + j;
      cur++;
      
      
      values[cur] = - dx / sh;
      column_index[cur] = i * n + (j + 1);
      cur++;
            
      values[cur] = - dx / sh;
      column_index[cur] = (i + 1) * n + j;
      cur++;
      
      row_index[k + 1] = row_index[k] + 4;      
    }
    if ( (j == bn) && (i >= 1) && (i < bn) )  // right border
    {      
      
      column_index[cur] = (i - 1) * n + j;
      values[cur] = - dy / sh;
      cur++;

      values[cur] = - dy / sh;
      column_index[cur] = i * n + (j - 1);
      cur++;

      values[cur] = 2 * (dx + dy) / sh;
      column_index[cur] = i * n + j;
      cur++;
            
      values[cur] = - dx / sh;
      column_index[cur] = (i + 1) * n + j;
      cur++;
      
      row_index[k + 1] = row_index[k] + 4;
    }
    if ( (i == bn) && (j >= 1) && (j < bn) )  // bottom border
    {
      column_index[cur] = (i - 1) * n + j;
      values[cur] = - dy / sh;
      cur++;

      values[cur] = - dy / sh;
      column_index[cur] = i * n + (j - 1);
      cur++;
      
      values[cur] = 2 * (dx + dy) / sh;
      column_index[cur] = i * n + j;
      cur++;
      
      values[cur] = - dx / sh;
      column_index[cur] = i * n + (j + 1);
      cur++;
    
      row_index[k + 1] = row_index[k] + 4;
    }
    if ( (j == 0) && (i >= 1) && (i < bn) )  // top border
    {
      column_index[cur] = (i - 1) * n + j;
      values[cur] = - dy / sh;
      cur++;

      values[cur] = 2 * (dx + dy) / sh;
      column_index[cur] = i * n + j;
      cur++;
      
      values[cur] = - dx / sh;
      column_index[cur] = i * n + (j + 1);
      cur++;
            
      values[cur] = - dx / sh;
      column_index[cur] = (i + 1) * n + j;
      cur++;
      
      row_index[k + 1] = row_index[k] + 4;

    }
    if ( (i == bn) && (j == bn) )
    {
      column_index[cur] = (i - 1) * n + j;
      values[cur] = - dy / sh;
      cur++;

      values[cur] = - dy / sh;
      column_index[cur] = i * n + (j - 1);
      cur++;
      

      values[cur] = 2 * (dx + dy) / sh;
      column_index[cur] = i * n + j;
      cur++;
      
      row_index[k + 1] = row_index[k] + 3;
    }
    if ( (i == 0) && (j == 0) )
    {
      values[cur] = 2 * (dx + dy) / sh;
      column_index[cur] = i * n + j;
      cur++;
      
      
      values[cur] = - dx / sh;
      column_index[cur] = i * n + (j + 1);
      cur++;
            
      values[cur] = - dx / sh;
      column_index[cur] = (i + 1) * n + j;
      cur++;
      
      row_index[k + 1] = row_index[k] + 3;
    } 
    if ( (i == 0) && (j == bn) )
    {
      values[cur] = - dy / sh;
      column_index[cur] = i * n + (j - 1);
      cur++;
      
      values[cur] = 2 * (dx + dy) / sh;
      column_index[cur] = i * n + j;
      cur++;
            
      values[cur] = - dx / sh;
      column_index[cur] = (i + 1) * n + j;
      cur++;
      
      row_index[k + 1] = row_index[k] + 3;
    }
    if ( (i == bn) && (j == 0) )
    {
      column_index[cur] = (i - 1) * n + j;
      values[cur] = - dy / sh;
      cur++;

      values[cur] = 2 * (dx + dy) / sh;
      column_index[cur] = i * n + j;
      cur++;
      
      values[cur] = - dx / sh;
      column_index[cur] = i * n + (j + 1);
      cur++;
      
      row_index[k + 1] = row_index[k] + 3;            
    }
    
  }
  for (int r = 0; r < NNZ; r++)
    cout << " " << column_index[r] << " ";  
  
  
  /*for (int r = 0; r <= N; r++)
    cout << r << " " << row_index[r] << "\n ";  
  */


  // Writing to file
  ofstream file("CSR.txt");
  file << N << " " << NNZ << endl;
  for (int k = 0; k < NNZ; ++k)
  {
    file << values[k] << " ";
  }
  file << endl;
  for (int k = 0; k < NNZ; ++k)
  {
    file << column_index[k] << " ";
  }
  file << endl;
  for (int i = 0; i < aN; ++i)
  {
    file << row_index[i] << " ";
  }
  file << endl;
  for (int i = 0; i < N; ++i)
  {
    file << b[i] << " ";
  }
  file << endl;
  file.close();


  //Factorizing  
  int start_fact = clock();
  double *u = new double[N];
  for (int i = 0; i < N; ++i)
  {
    u[i] = 0;
  }
  void *Symbolic, *Numeric ;
  umfpack_di_symbolic (N, N, row_index, column_index, values, &Symbolic, nullptr, nullptr);
  umfpack_di_numeric (row_index, column_index, values, Symbolic, &Numeric, nullptr, nullptr);
  double time_fact = (double) (clock() - start_fact) / CLOCKS_PER_SEC;
  cout << "time of factorization: ms " << time_fact * 1e3 << endl;

  // Solving
  int start_solve = clock();
  umfpack_di_solve (UMFPACK_A, row_index, column_index, values, u, b, Numeric, nullptr, nullptr);
  double time_solve = (double )(clock() - start_solve) / CLOCKS_PER_SEC;

  umfpack_di_free_symbolic (&Symbolic);
  umfpack_di_free_numeric (&Numeric);  
  cout << "time of solving: ms " << time_solve * 1e3 << "\n\n";

  //Checking
  for (int h = 0; h < n; h++)
  {
    uint32_t hght = h * n;
    cout << "[";
    for (int w = 0; w < n - 1; w++)
    {
      cout << u[hght + w] << ", ";
    }
    cout << u[hght + n - 1] << " ";
    cout<< "], \n";
  }

  delete[] b;
  delete[] u;
  delete[] values;
  delete[] column_index;
  delete[] row_index;
  
}
