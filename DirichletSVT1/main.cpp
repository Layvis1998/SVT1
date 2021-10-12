#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "../Include/umfpack.h"
using namespace std;

int n = 30;        // i = 0..29, j = 0..29
int bn = n - 1;
double h = 1.0 / n;
double sh = h * h;
int N = n * n;  
int aN = N + 1;

int* row_index = new int[aN];
int NNZ = 5 * (n - 2)*(n - 2) + 4 *(4 * (n - 2)) + 4 * (3);
                       
int* column_index = new int[NNZ];
double* values = new double[NNZ];
double* b = new double[N];
int cur = 0;
  
int dx = 1;
int dy = 10;
int* G_D_top = new int[bn];
int* G_D_right = new int[bn];
    
double f(double x, double y)
{
  return -12.0 * (pow(x - y, 2) + 1e-5) * (5.0 * pow(x - y, 2) + 1e-5);
}

double u_answer(double x, double y)
{
  return pow((x - y) * (x - y) + 1e-5, 3);
}

void No_borders(int i, int j, int k)
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

void Top_Dirichlet(int i, int j, int k)
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

void Right_Dirichlet(int i, int j, int k)
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

void Left_Dirichlet(int i, int j, int k)
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

void Bottom_Dirichlet(int i, int j, int k)
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

void BR_Dirichlet(int i, int j, int k)
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

void TL_Dirichlet(int i, int j, int k)
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
void BL_Dirichlet(int i, int j, int k)
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

void TR_Dirichlet(int i, int j, int k)
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

int main()
{
  for (int a = 0; a < n; a++)
  {
    G_D_top[a] = 100;
  }
  for (int a = 0; a < n; a++)
  {
    G_D_right[a] = 100;
  }
  
  for (int k = 0; k < N; k++)
  {
    int i =  k  / n;    
    int j = k % n;
    b[k] = cos(M_PI * i * h) * cos(M_PI * j * h);
  }

  row_index[0] = 0;  
  for (int k = 0; k < N; k++)
  {
    int i = k / n;
    int j = k % n;
    
    if ( (i < bn) && (j < bn) && (i > 0) && (j > 0) ) 
      No_borders(i, j, k);
    
    if ( (i == 0) && (j >= 1) && (j < bn) )
      Top_Dirichlet(i, j, k);
    
    if ( (j == bn) && (i >= 1) && (i < bn) )
      Right_Dirichlet(i, j, k);   
    
    if ( (j == 0) && (i >= 1) && (i < bn) ) 
      Left_Dirichlet(i, j, k);
    
    if ( (i == bn) && (j >= 1) && (j < bn) )
      Bottom_Dirichlet(i, j, k);

    if ( (i == bn) && (j == bn) )  // bottom right corner
      BR_Dirichlet(i, j, k);
        
    if ( (i == 0) && (j == 0) )  //top left corner      
      TL_Dirichlet(i, j, k);
     
    if ( (i == 0) && (j == bn) ) //bottom left corner
      BL_Dirichlet(i, j, k);
    
    if ( (i == bn) && (j == 0) )  // bottom left corner
      TR_Dirichlet(i, j, k);
  }
    
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
  cout << "Solution:\n ";
  for (int h = 0; h < n; h++)
  {
    uint32_t hght = h * n;
    cout << "[";
    for (int w = 0; w < n - 1; w++)
    {
      cout << u[hght + w] << ", ";
    }
    cout << u[hght + n - 1] << " ";
    cout<< "] \n";
  }

  delete[] b;
  delete[] u;
  delete[] values;
  delete[] column_index;
  delete[] row_index;
}
