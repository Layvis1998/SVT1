#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "Include/umfpack.h"
using namespace std;

int Neumann_type = 1;   //type of Neumann borders

int n = 10;        
int bn = n - 1;
int bbn = bn - 1;
int bbbn = bbn - 1;
double h = 1.0 / n;
double sh = h * h;
int N = n * n;  
int aN = N + 1;
int b_N = bn * bn;

int NNZ;
int* row_index;                       
int* column_index;
double* values;
int cur = 0;
  
int dx = 1;
int dy = 1;
    
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
  values[cur] = - dx / sh;
  cur++;

  values[cur] = - dy / sh;
  column_index[cur] = i * n + (j - 1);
  cur++;

  values[cur] = 2 * (dx + dy) / sh;
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy / sh;
  column_index[cur] = i * n + (j + 1);
  cur++;
            
  values[cur] = - dx / sh;
  column_index[cur] = (i + 1) * n + j;
  cur++;
      
  row_index[k + 1] = row_index[k] + 5;
}


void TR_Dirichlet(int i, int j, int k)
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

void Top_Dirichlet(int i, int j, int k)
{

  values[cur] = - dx / sh;
  column_index[cur] = i * n + (j - 1);
  cur++;
      
  values[cur] = 2 * (dx + dy) / sh;
  column_index[cur] = i * n + j;
  cur++;
      
  values[cur] = - dy / sh;
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
  values[cur] = - dx / sh;
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



void Left_Neumann1(int i, int j, int k)
{
  column_index[cur] = i * n;
  values[cur] = - dy / sh;
  cur++;

  column_index[cur] = i * n + 1;
  values[cur] =  dy / sh;
  cur++;
  
  row_index[k + 1] = row_index[k] + 2; 
}

void Bottom_Neumann1(int i, int j, int k)
{
  column_index[cur] = bbn * n + j;
  values[cur] = dx / sh;
  cur++;

  column_index[cur] = bn * n + j;
  values[cur] =  - dx / sh;
  cur++;
  
  row_index[k + 1] = row_index[k] + 2; 
}


void BL_Neumann1(int i, int j, int k)
{ 
  
  column_index[cur] = bn * n;
  values[cur] = - (dx + dy) / sh / sqrt(2);
  cur++;
  
  row_index[k + 1] = row_index[k] + 1;
}

void Known_Values(int i, int j, int k)
{
  values[cur] = 1;
  column_index[cur] = i * n + j;
  cur++;
  row_index[k + 1] = row_index[k] + 1; 
}


int main()
{
  double* G_Dirichlet_top = new double[bbn];
  double* G_Dirichlet_right = new double[bbn];
  double* G_Neumann_left = new double[bbn];
  double* G_Neumann_bottom = new double[bbn];
  double G_BL = 2000;
  double G_TR = 0;
  double G_TL = 2000;
  double G_BR = 0;
  double* b = new double[N];   
  
  row_index = new int[aN]; 

  for (int k = 0; k < N; k++)
    b[k] = 0;
  
  for (int a = 0; a < bbn; a++)
    G_Dirichlet_top[a] = 0;
    
  for (int a = 0; a < bbn; a++)
    G_Dirichlet_right[a] = 0;

  for (int a = 0; a < bbn; a++)
    G_Neumann_left[a] = 2000;
  
  for (int a = 0; a < bbn; a++)
    G_Neumann_bottom[a] = 0;
  
  for (int a = 1; a <= bbn; a++)
    b[a * n] = G_Neumann_left[a - 1] / h;
  
  for (int a = 1; a <= bbn; a++)
    b[bn*n + a] = G_Neumann_bottom[a - 1] / h; 

  for (int a = 1; a <= bbn; a++)                       
    b[a * n + bn] += G_Dirichlet_right[a - 1] * dy / sh;
    
  for (int a = 1; a <= bbn; a++)
    b[a] += G_Dirichlet_top[a - 1] * dx / sh;

  b[0] = G_TL;
  b[bn] = G_TR;
  b[n * n - 1] = G_BR;
  b[bn * n] = G_BL / h;

  NNZ = 5 * bbn*bbn + 2*bbn + 2*bbn + 4*bbn + 4*bbn + 3 + (bn + bn + 1);  
  column_index = new int[NNZ];
  values = new double[NNZ];  
  row_index[0] = 0;  
  

  for (int k = 0; k < N; k++)
  {
    int i = k / n;
    int j = k % n;
        
    if ( (i == 0) && (j >= 0) && (j <= bn) )
      Known_Values(i, j, k);
    
    if ( (j == bn) && (i >= 1) && (i <= bn) )     
      Known_Values(i, j, k);
    
    if ( (i > 1) && (i < bn) && (j > 0) && (j < bbn) ) 
      No_borders(i, j, k);
    
    if ( (i == 1) && (j >= 1) && (j < bbn) )
      Top_Dirichlet(i, j, k);
    
    if ( (j == bbn) && (i > 1) && (i <= bbn) )
      Right_Dirichlet(i, j, k);   
    
    if ( (i == 1) && (j == bbn) )
      TR_Dirichlet(i, j, k);

    if ( (j == 0) && (i >= 1) && (i <= bbn) )
      Left_Neumann1(i, j, k);
       
    if ( (i == bn) && (j >= 1) && (j <= bbn) )
      Bottom_Neumann1(i, j, k);
      
    if ( (j == 0) && (i == bn) )
      BL_Neumann1(i, j, k);
  }
  

  // Writing to file
  
  ofstream file("CSR2.txt");
  file << N << " " << NNZ << endl;
  

  for (int k = 0; k < NNZ; ++k)
    file << values[k] << " ";    
  file << endl;
  
  for (int k = 0; k < NNZ; ++k)
    file << column_index[k] << " ";    
  file << endl;
  
  for (int i = 0; i < aN; ++i)
    file << row_index[i] << " ";
  file << endl;
  
  for (int i = 0; i < N; ++i)
    file << b[i] << " ";    
  file << endl;
  
  file.close();


  //Factorizing  
  int start_fact = clock();
  double *u = new double[N];
  for (int i = 0; i < N; ++i)
    u[i] = 0;
    
  void *Symbolic, *Numeric ;
  umfpack_di_symbolic (N, N, row_index, column_index, values, &Symbolic, nullptr, nullptr);
  umfpack_di_numeric (row_index, column_index, values, Symbolic, &Numeric, nullptr, nullptr);
  double time_fact = (double) (clock() - start_fact) / CLOCKS_PER_SEC;
  cout << "time of factorization: ms " << time_fact * 1e3 << endl;

  //Solving
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
      cout << u[hght + w] << ", ";
    
    cout << u[hght + n - 1] << " ";
    cout<< "] \n";
  }

}
