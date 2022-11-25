#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
//#include <mpi.h>

#define lx 1.0
#define ly 1.0
#define lz 1.0
#define Tau 1.0
#define N 127
#define K 220 

double A[N + 1][N + 1][N + 1];
double B[N + 1][N + 1][N + 1];
double C[N + 1][N + 1][N + 1];

double fi (double x, double y, double z){
  return sin(M_PI / lx * ((lx / N) * x)) * sin(2 * M_PI / ly * ((ly / N) * y)) * sin(3 * M_PI / lz * ((lz / N) * z));
}

double u (double x, double y, double z, double t){
  return fi(x, y ,z) * cos(sqrt(1/(lx * lx) + 4/(ly * ly) + 9/(lz * lz)) * ((Tau / K) * t));
} 

void u_analytic_init (double (*a)[N + 1][N + 1], int n, double t){
  for (int i = 0; i <= n; i++){
    for (int j = 0; j <= n; j++){
      for (int k = 0; k <= n; k++){
        //printf("%8.4lf, %8.4lf, %8.4lf, %8.4lf ||| %8.4lf , %8.4lf, %8.4lf, %8.4lf\n",(lx / n) * i,(ly / n) * j,(lz / n) * k, t, sin(M_PI / lx * (lx / n) * i), sin(2 * M_PI / ly * (ly / n) * j), sin(3 * M_PI / lz * (lz / n) * k), cos(sqrt(1/(lx * lx) + 4/(ly * ly) + 9/(lz * lz)) * t));
        a[i][j][k] = u(i, j, k, t);
      }
    }
  }
}

void u_0_init (double (*a)[N + 1][N + 1], int n){
  for (int i = 0; i <= n; i++){
    for (int j = 0; j <= n; j++){
      for (int k = 0; k <= n; k++){
        a[i][j][k] = fi(i, j, k);
      }
    }
  }
}

void u_1_init (double a[N + 1][N + 1][N + 1], double pa[N + 1][N + 1][N + 1], int n){
  for (int i = 0; i <= n; i++){
    for (int j = 0; j <= n; j++){
      for (int k = 0; k <= n; k++){
        if (i == 0 || i == n){
          a[i][j][k] = 0;
        } else if (k == 0 || k == n){
          a[i][j][k] = 0;
        } else {
          a[i][j][k] = pa[i][j][k] + ((Tau / K) * (Tau / K) / 2) * (
            ((fi(i-1, j, k) - 2 * (fi(i, j, k)) + fi(i+1,j,k)) / ((lx / n) * (lx / n))) +
            ((fi(i, j-1, k) - 2 * (fi(i, j, k)) + fi(i,j+1,k)) / ((ly / n) * (ly / n))) +
            ((fi(i, j, k-1) - 2 * (fi(i, j, k)) + fi(i,j,k+1)) / ((lz / n) * (lz / n))) );
        }
      }
    }
  }
  return;
}

void u_step (double (*a)[N + 1][N + 1], double (*pa)[N + 1][N + 1], double (*ppa)[N + 1][N + 1], int n){

  for (int i = 0; i <= n; i++){
    for (int j = 0; j <= n; j++){
      a[i][j][0] = 0;
      a[i][j][n] = 0;
    }
  }

  for (int j = 0; j <= n; j++){
    for (int k = 0; k <= n; k++){
      a[0][j][k] = 0;
      a[n][j][k] = 0;
    }
  }

  for (int i = 1; i < n; i++){
    for (int k = 1; k < n; k++){ 
      for (int j = 1; j < n; j++){
        a[i][j][k] = 2 * pa[i][j][k] - ppa[i][j][k] + ((Tau / K) * (Tau / K)) * (
          ((pa[i-1][j][k] - 2 * pa[i][j][k] + pa[i+1][j][k]) / ((lx / n) * (lx / n))) +
          ((pa[i][j-1][k] - 2 * pa[i][j][k] + pa[i][j+1][k]) / ((ly / n) * (ly / n))) +
          ((pa[i][j][k-1] - 2 * pa[i][j][k] + pa[i][j][k+1]) / ((lz / n) * (lz / n)))
        );
      }
      a[i][n][k] = 2 * pa[i][n][k] - ppa[i][n][k] + ((Tau / K) * (Tau / K)) * (
          ((pa[i-1][n][k] - 2 * pa[i][n][k] + pa[i+1][n][k]) / ((lx / n) * (lx / n))) +
          ((pa[i][n-1][k] - 2 * pa[i][n][k] + pa[i][1][k]) / ((ly / n) * (ly / n))) +
          ((pa[i][n][k-1] - 2 * pa[i][n][k] + pa[i][n][k+1]) / ((lz / n) * (lz / n)))
        );
      a[i][0][k] = a[i][n][k];
    }
  }

  return;
}

double u_err(double (*a)[N + 1][N + 1], int n, double t){
  double max = 0;
  double tmp;
   for (int i = 0; i <= n; i++){
    for (int j = 0; j <= n; j++){
      for (int k = 0; k <= n; k++){
        tmp = fabs(a[i][j][k] - u(i,j,k,t));
        if (tmp > max)
          max = tmp;
      }
    }
  }   
  return max;
}



void print_grid (double (*a)[N + 1][N + 1], int n){
  printf("-----------------------------------------------------\n");
  for (int i = 0; i <= n; i++){
    printf("x = %d \n", i);
    for (int j = 0; j <= n; j++){
      for (int k = 0; k <= n; k++){
        printf("%8.4lf, ",a[i][j][k]);
      }
      printf("\n");
    }
  }
  printf("-----------------------------------------------------\n");
}


int main(){
  double err = 0;

  // INIT A, t = 0
  u_0_init(A, N);
  err = u_err(A, N, 0);
  //print_grid(A, N);
  printf("ERROR[%d]::%lf;\n",0,err);
  
  // INIT B, t = 1 
  u_1_init(B, A, N);
  err = u_err(B, N, 1);
  //print_grid(B, N);
  printf("ERROR[%d]::%lf;\n",1,err);

  for(int s = 2; s <= K ; s++){
    if (s % 3 == 2){
      u_step(C, B, A, N);
      err = u_err(C, N, s);
      printf("ERROR[%d]::%lf;\n",s,err);
    } else if (s % 3 == 1){
      u_step(B, A, C, N);
      err = u_err(B, N, s);
      printf("ERROR[%d]::%lf;\n",s,err);
    } else if (s % 3 == 0){
      u_step(A, C, B, N);
      err = u_err(A, N, s);
      printf("ERROR[%d]::%lf;\n",s,err);
    }
  }

  return 0;
}
