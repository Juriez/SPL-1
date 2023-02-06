
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ROW 3
#define COL 3

int main (void) {
  int i, j, k;
  double A[ROW][COL] ;//= { {1,3,2}, {1,2,3}, {2,-1,4}};
  double b[ROW] ;//= {17,16,13};
  printf("Enter the matrix :\n");
  for(i =0;i<ROW;i++)
  {
      for(j=0;j<COL;j++)
      {
          scanf("%lf",&A[i][j]);
      }
  }
  printf("Enter the column matrix :\n");
  for(i=0;i<ROW;i++)
  {
      scanf("%lf",&b[i]);
  }
  double Ab[ROW][COL+1];
  double M, L;

  for (i=0; i< ROW; i++) {
    for (j=0; j< COL; j++) {
      Ab[i][j] = A[i][j];
    }
    Ab[i][j] = b[i];
  }

  printf("[ [A] [b] ] :\n");
  for (i=0; i< ROW; i++) {
    for (j=0; j< COL+1; j++) {
      printf("%9.3lf", Ab[i][j] );
    }
    printf("\n");
  }


  for (i=0; i<ROW; i++) {
    L = Ab[i][i];
    for (j=i; j<COL+1; j++) {
        Ab[i][j] = Ab[i][j]/L;
    }

    for (k=0; k<ROW; k++) {
      if (k != i) {
        M = -Ab[k][i];
        for (j=i; j<COL+1; j++) {
          Ab[k][j] = Ab[k][j] + M*Ab[i][j];
        }
      }
    }
  }


  printf("[ [A'] [b'] ] :\n");
  for (i=0; i< ROW; i++) {
    for (j=0; j< COL+1; j++) {
      printf("%9.3lf", Ab[i][j] );
    }
    printf("\n");
  }
  printf("So,tha value of three variables X , Y , Z :");

for (i=0; i< ROW; i++) {
    printf("%9.3lf ", Ab[i][COL] );
  }
  return EXIT_SUCCESS;
}

