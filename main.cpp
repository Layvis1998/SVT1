#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "Include/umfpack.h"
using namespace std;

int Neumann_type = 2;   //type of Neumann borders

int n = 30;        
int bn = n - 1;
int bbn = bn - 1;
double h = 1.0 / n;
double sh = h * h;
int N = n * n;  
int aN = N + 1;

int NNZ;
int* row_index;                       
int* column_index;
double* values;
double* b;
int cur = 0;
  
int dx = 1;
int dy = 10;
double* G_Dirichlet_top;
double* G_Dirichlet_right;
double* G_Neumann_left;
double* G_Neumann_bottom;
    
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

void BL_Dirichlet(int i, int j, int k)
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

void Left_Neumann1(int i, int j, int k)
{
  column_index[cur] = 0 + (i - 1) * n;
  values[cur] = dx / h;
  cur++;

  column_index[cur] = 0 + i * n;
  values[cur] = - dx / h;
  cur++;
  
  row_index[k + 1] = row_index[k] + 2; 
}

void Bottom_Neumann1(int i, int j, int k)
{
  column_index[cur] = bn * n + (j - 1);
  values[cur] = dy / h;
  cur++;

  column_index[cur] = bn * n + j;
  values[cur] = - dy / h;
  cur++;
  
  row_index[k + 1] = row_index[k] + 2; 
}

void Bottom_Neumann2(int i, int j, int k)
{  
  column_index[cur] = bbn * n;
  values[cur] = dy / sh;
  cur++;
  
  column_index[cur] = bn * n - 1;
  values[cur] = dx / sh / 2;
  cur++;
  
  column_index[cur] = bn * n;
  values[cur] = (dy + dx) / sh;
  cur++;
  
  column_index[cur] = bn * n + 1;
  values[cur] = dx / sh / 2;
  cur++;
  
  row_index[k + 1] = row_index[k] + 4;
}

void BL_Neumann2(int i, int j, int k)
{ 
  column_index[cur] = bbn * n;
  values[cur] = -dy / sh;
  cur++;
  
  column_index[cur] = bn * n - 1;
  values[cur] = -dx / sh;
  cur++;
  
  column_index[cur] = bn * n;
  values[cur] = (dy + dx) / sh;
  cur++;
    
  row_index[k + 1] = row_index[k] + 3;
}

void Left_Neumann2(int i, int j, int k)
{
  column_index[cur] = (i - 1) * n;
  values[cur] = - dy / sh / 2;
  cur++;

  column_index[cur] = i * n;
  values[cur] = (dx + dy) / sh;
  cur++;

  column_index[cur] = i * n + 1;
  values[cur] = -dx / sh;
  cur++;
  
  column_index[cur] = (i + 1) * n;
  values[cur] = - dy / sh / 2;
  cur++;
  
  row_index[k + 1] = row_index[k] + 4; 
}

int main()
{
  G_Dirichlet_top = new double[n];
  G_Dirichlet_right = new double[bn];
  G_Neumann_left = new double[bn];
  G_Neumann_bottom = new double[bbn];
  row_index = new int[aN];
  b = new double[N];    

  for (int k = 0; k < N; k++)
  {
    int i =  k  / n;    
    int j = k % n;
    b[k] = cos(M_PI * i * h) * cos(M_PI * j * h);
  }
  
  for (int a = 0; a < bn; a++)
    G_Dirichlet_top[a] = 100;
    
  for (int a = 0; a < bbn; a++)
    G_Dirichlet_right[a] = 100;

  for (int a = 0; a < bn; a++)
    G_Neumann_left[a] = 1;
  
  for (int a = 0; a < bbn; a++)
    G_Neumann_bottom[a] = 1.2;
  
  if (Neumann_type == 1)
  {
    NNZ = 5 * bbn*bbn + 2*bn + 2*bbn + 4*bbn + 4*bbn + 3 + 3 + 3;  
    for (int a = 0; a < n; a++)
      b[a * n] = G_Neumann_left[a];
      
    for (int a = 1; a < bn; a++)
      b[bn*n + a] = G_Neumann_bottom[a - 1]; 
  }
  if (Neumann_type == 2)
  {
    NNZ = 5 * bbn*bbn + 4*bn + 4*bbn + 4*bbn + 4*bbn + 3 + 3 + 3;  
    for (int a = 0; a < n; a++)
      b[a * n] = b[a * n] / 2 + G_Neumann_left[a] / h; 
    
    for (int a = 1; a < bn; a++)      
      b[bn*n + a] = b[bn*n + a]/2 + G_Neumann_bottom[a - 1] / h; 
  }
                     
  column_index = new int[NNZ];
  values = new double[NNZ];
  
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
    
    if ( (j == 0) && (i >= 1) && (i <= bn) )
    {
      if (Neumann_type == 1)
        Left_Neumann1(i, j, k);
      else if (i < bn)
        Left_Neumann2(i, j, k);
    }
    
    if ( (j == 0) && (i == bn) && (Neumann_type != 1))
      BL_Neumann2(i, j, k);  
      
    if ( (i == bn) && (j >= 1) && (j < bn) )
    {
      if (Neumann_type == 1)
        Bottom_Neumann1(i, j, k);
      else
        Bottom_Neumann2(i, j, k);
    }
    if ( (i == bn) && (j == bn) )  // bottom right corner
      BR_Dirichlet(i, j, k);
        
    if ( (i == 0) && (j == 0) )  //top left corner      
      TL_Dirichlet(i, j, k);
    if ( (i == 0) && (j == bn) )  // bottom left corner
      TR_Dirichlet(i, j, k);
  }
  
  // Writing to file
  ofstream file("CSR.txt");
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
      cout << u[hght + w] << ", ";
    
    cout << u[hght + n - 1] << " ";
    cout<< "] \n";
  }


}
