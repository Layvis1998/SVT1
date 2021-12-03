#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "Include/umfpack.h"
using namespace std;

int Neumann_type = 1;   //type of Neumann borders

int n = 25;        
int bn = n - 1;
int bbn = bn - 1;
int bbbn = bbn - 1;
double h = 1.0 / bn;
double sh = h * h;
int N = n * n;  
int aN = N + 1;
int b_N = bn * bn;
  
double dx = 1;
double dy = 10;

int NNZ;
int* row_index;                       
int* column_index;
double* values;
int cur = 0;

void No_borders(int i, int j, int k)
{
  values[cur] = - dx / sh;
  column_index[cur] = (i - 1) * n + j;
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
      
  row_index[k + 1] = row_index[k] + 4; 
}

void Right_Dirichlet(int i, int j, int k)
{
  values[cur] = - dx / sh;
  column_index[cur] = (i - 1) * n + j;
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
  values[cur] = - dy / sh;
  column_index[cur] = i * n + j;
  cur++;

  values[cur] =  dy / sh;
  column_index[cur] = i * n + (j + 1);
  cur++;
  
  row_index[k + 1] = row_index[k] + 2; 
}

void Bottom_Neumann1(int i, int j, int k)
{
  values[cur] = dx / sh;
  column_index[cur] = (i - 1) * n + j;
  cur++;

  values[cur] =  - dx / sh;
  column_index[cur] = i * n + j;
  cur++;
  
  row_index[k + 1] = row_index[k] + 2; 
}


void BL_Neumann1(int i, int j, int k)
{ 

  values[cur] = (dx + dy) / sh  * (1 / sqrt(2));
  column_index[cur] = bbn * n + 1;
  cur++;

  values[cur] = - (dx + dy) / sh * (1 / sqrt(2));
  column_index[cur] = bn * n;
  cur++;

  row_index[k + 1] = row_index[k] + 2;
}


int main()
{
  double* b = new double[N];   
  double* G_Dirichlet_top = new double[bbn];
  double* G_Dirichlet_right = new double[bbn];
  double* G_Neumann_left = new double[bbn];
  double* G_Neumann_bottom = new double[bbn];

  double G_BL =  (- (1/sqrt(2)) * dx * sin(M_PI * bn * h) * cos(0) * M_PI
    +(1/sqrt(2)) * dy * sin(0) * cos(M_PI * bn * h)  * M_PI);
  double G_TL = cos(0) * cos(0);   
  double G_TR = cos(0) * cos(M_PI * bn * h);
  double G_BR = cos(M_PI * bn * h) * cos(M_PI * bn* h);


  row_index = new int[aN]; 

  for (int k = 0; k < N; k++)
  {  
    int i = k / n;
    int j = k % n;
    b[k] =  (dx + dy) * cos(M_PI * i * h) * cos(M_PI * j * h) * M_PI * M_PI;
  }


  for (int a = 0; a < bbn; a++)
    G_Dirichlet_top[a] = cos(0) * cos(M_PI * (a + 1) * h);
    
  for (int a = 0; a < bbn; a++)
    G_Dirichlet_right[a] = cos(M_PI * (a + 1) * h ) * cos(M_PI * bn * h);

  for (int a = 0; a < bbn; a++)
    G_Neumann_left[a] = dy * sin(0) * cos(M_PI * (a + 1) * h) * M_PI;
  
  for (int a = 0; a < bbn; a++)
    G_Neumann_bottom[a] = -dx * sin(M_PI * bn * h) * cos(M_PI * (a + 1) * h) * M_PI;
  
  for (int a = 1; a <= bbn; a++)
    b[a * n] = G_Neumann_left[a - 1] / h;
  
  for (int a = 1; a <= bbn; a++)
    b[bn*n + a] = G_Neumann_bottom[a - 1] / h; 
  

  for (int a = 1; a <= bbn; a++)                       
    b[a * n + bn] = G_Dirichlet_right[a - 1];
    
  for (int a = 1; a <= bbn; a++)
    b[a] = G_Dirichlet_top[a - 1];


  b[0] = G_TL;
  b[bn] = G_TR;
  b[n * n - 1] = G_BR;
  b[bn * n] = G_BL / h;
   
  for (int a = 1; a <= bbn; a++)
    b[n + a] += (b[a] / sh * dx);

  for (int a = 1; a <= bbn; a++)                       
    b[a * n + bbn] += b[a * n + bn] / sh * dy;
    

  NNZ = (5 * bbn*bbn + 2*bbn + 2*bbn + 4*bbn + 4*bbn + 3) + (bn + bn + 2);               
  column_index = new int[NNZ];                                                
  values = new double[NNZ];  
  row_index[0] = 0;  
  

  for (int k = 0; k < N; k++)
  {
    int i = k / n;
    int j = k % n;
        
    if ( (i == 0) && (j >= 1) && (j <= bbn) )
    {
      values[cur] = 1;
      column_index[cur] = i * n + j;
      cur++;
      row_index[k + 1] = row_index[k] + 1; 
    }  
    
    if ( (j == bn) && (i >= 1) && (i < bn) )     
    {
      values[cur] = 1;
      column_index[cur] = i * n + j;
      cur++;
      row_index[k + 1] = row_index[k] + 1; 
    }
    
    if ( (j == 0) && (i == 0) )     
    {
      values[cur] = 1;
      column_index[cur] = i * n + j;
      cur++;
      row_index[k + 1] = row_index[k] + 1; 
    }

    if ( (j == bn) && (i == 0) )     
    {
      values[cur] = 1;
      column_index[cur] = i * n + j;
      cur++;
      row_index[k + 1] = row_index[k] + 1; 
    }

    if ( (j == bn) && (i == bn) )     
    {
      values[cur] = 1;
      column_index[cur] = i * n + j;
      cur++;
      row_index[k + 1] = row_index[k] + 1; 
    }


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
  
  /*
  //Printing solution
  cout << "\nSolution:\n\n";
  for (int h = 0; h < n; h++)
  {
    uint32_t hght = h * n;
    cout << "[";
    for (int w = 0; w < n - 1; w++)
      cout << u[hght + w] << ", ";
    
    cout << u[hght + n - 1] << " ";
    cout<< "] \n";
  }
  */


  //Checking Ch-norm
  double Ch = 2;
  for (int i = 1; i < bn; i++)
  {
    uint32_t hght = i * n;
    for (int j = 1; j < bn; j++)
    {
      double val = cos(M_PI * i * h) * cos(M_PI * j * h);
      val = abs(val - u[hght + j]);
      if (val < Ch)
        Ch = val; 
    } 
  }
  cout << "\n\nCh norm = " << Ch << "\n";

  //Checking L2h-norm
  double Lh = 0;
  for (int i = 1; i < bn; i++)
  {
    uint32_t hght = i * n;
    for (int j = 1; j < bn; j++)
    {
      double val = cos(M_PI * i * h) * cos(M_PI * j * h);
      val = u[hght + j] - val;
      val = val * val;
      val = val * sh;

      if ((i > 1) && (j > 1) && (i < bbn) && (j < bbn))
        Lh += val;

      if ((i == 1) && (j == 1))
        Lh += 1/4 * val;
      if ((i == bbn) && (j == bbn))
        Lh += 1/4 * val;
      if ((i == bbn) && (j == 1))
        Lh += 1/4 * val;
      if ((i == 1) && (j == bbn))
        Lh += 1/4 * val;

      if ((i > 1) && (i < bbn) && (j == 1))
        Lh += 1/2 * val;
      if ((i > 1) && (i < bbn) && (j == bbn))
        Lh += 1/2 * val;    
      if ((j > 1) && (j < bbn) && (i == 1))
        Lh += 1/2 * val;
      if ((j > 1) && (j < bbn) && (i == bbn))
        Lh += 1/2 * val;           

    } 
  
  }
  Lh = sqrt(Lh);

  cout << "Lh norm = " << Lh << "\n";


}
