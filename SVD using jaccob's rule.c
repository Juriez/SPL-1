
#include <stdio.h>
#include <math.h>
double A[100][100],A_t[100][100],S[100][100],SI[100][100],eig_val[100],V[100][100],Vt[100][100],U[100][100],VtSi[100][100],Ut[100][100],SV[100][100],SVD[100][100],SVDt[100][100];

int l=0; int k=0;
int m,n;
void main() {
    int  i, j, p, q, flag;
    double d[100][100], s[100][100], s1[100][100], mul[100][100], s1t[100][100];
    double temp[100][100], theta, zero=1e-5, max, pi=3.141592654;

  printf("Enter the rows & columns of the matrix ");
  scanf("%d %d",&n,&m);

  printf("Enter the elements row wise \n");
  for(i=0; i<n; i++) {
    for(j=0; j<m; j++)  scanf ("%lf", &A[i][j]);
  }

  printf("The given matrix is\n");
  for(i=0; i<n; i++) {
    for(j=0; j<m; j++)  printf ("%8.4f ", A[i][j]);
    printf ("\n");
  }
  printf ("\n");
  printf("\nTranspose of the matrix :\n");
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
        A_t[j][i]=A[i][j];
    }
  }

  for(i=0; i<m; i++) {
    for(j=0; j<n; j++)  printf ("%lf ", A_t[i][j]);
    printf ("\n");

  }printf ("\n");
  multiplication(A_t,A,mul,m,n);
  printf("\nMultiplication :\n");
  for( i=0;i<m;i++)
         {
          for( j=0;j<m;j++)
           {
            printf("%lf\t",mul[i][j]);
           }
           printf("\n");
         }

  for(i=0; i<m; i++) {
    for(j=0; j<m; j++) {
      d[i][j]= mul[i][j];
      if(i==j)
        s[i][j]= 1;
      else
        s[i][j]= 0;
    }
  }

  do {
    flag= 0;
    p=0; q=1;
    max= fabs(d[p][q]);

    for(i=0; i<m; i++) {
      for(j=0; j<m; j++) {
        if(i!=j) {
          if (max < fabs(d[i][j])) {
            max= fabs(d[i][j]);
            p= i;
            q= j;
          }
        }
      }
    }

    if(d[p][p]==d[q][q]) {
      if (d[p][q] > 0)
        theta= pi/4;
      else
        theta= -pi/4;
    }
    else {
      theta=0.5*atan(2*d[p][q]/(d[p][p]-d[q][q]));
    }

    for(i=0; i<m; i++) {
      for(j=0; j<m; j++) {
        if(i==j) {
         s1[i][j]= 1;
         s1t[i][j]= 1;
        }
        else {
          s1[i][j]= 0;
          s1t[i][j]= 0;
        }
      }
    }

    s1[p][p]= cos(theta);
    s1t[p][p]= s1[p][p];

    s1[q][q]= cos(theta);
    s1t[q][q]= s1[q][q];

    s1[p][q]= -sin(theta);
    s1[q][p]= sin(theta);

    s1t[p][q]= s1[q][p];
    s1t[q][p]= s1[p][q];

    for(i=0; i<m; i++) {
      for(j=0; j<m; j++) {
        temp[i][j]= 0;
        for(p=0; p<n; p++)  temp[i][j]+= s1t[i][p]*d[p][j];
      }
    }

    for(i=0; i<m; i++) {
      for(j=0; j<m; j++) {
        d[i][j]= 0;
        for(p=0; p<n; p++)  d[i][j]+= temp[i][p]*s1[p][j];
      }
    }

    for(i=0; i<m; i++) {
      for(j=0; j<m; j++) {
        temp[i][j]= 0;
        for(p=0; p<n; p++)  temp[i][j]+= s[i][p]*s1[p][j];
      }
    }

    for(i=0; i<m; i++) {
      for(j=0; j<m; j++)  s[i][j]= temp[i][j];
    }

    for(i=0; i<m; i++) {
      for(j=0; j<m; j++) {
        if(i!=j)
          if(fabs(d[i][j]) > zero) flag= 1;
      }
    }
  } while(flag==1);

  printf("The eigenvalues are \n");
  for(i=0; i<m; i++){
    eig_val[i]=d[i][i];
    printf("%8.4lf ",d[i][i]);

  }
    for(i=0;i<m;i++){
        for(j=i+1;j<m;j++){
            if(eig_val[i]<eig_val[j]){
                double temp =0;
                temp = eig_val[i];
                eig_val[i]=eig_val[j];
                eig_val[j]=temp;
            }
        }
    }

    printf("\n Sorted eigen values are :\n");
    for(i=0; i<m; i++){
    printf("%8.4lf ",eig_val[i]);}
    for(i=0;i<m;i++){
        S[i][i] = sqrt(eig_val[i]);
        SI[i][i] = 1/S[i][i];
    }


       while(k < m) {
            double w2[100][100];
    for(int i=0; i<m; i++)
    {
        for(int j=0; j<m; j++)
        {
            w2[i][j]=mul[i][j];
        }
    }

    for(int j=0; j<m; j++)
    {
        w2[j][j] = mul[j][j]-eig_val[k];
    }
    gaussianElimination(w2);
    k++;
    }
    for(i=0;i<m;i++){
        for(j=0;j<m;j++){
            Vt[j][i]=V[i][j];
        }
    }


  multiplication(Vt,SI,VtSi,m,m);
  multiplication(A,VtSi,U,m,m);
  for(i=0;i<m;i++){
        for(j=0;j<m;j++){
            Ut[j][i]=U[i][j];
        }
    }
  multiplication(S,V,SV,m,m);
  multiplication(U,SV,SVDt,m,m);
  for(i=0;i<m;i++){
        for(j=0;j<m;j++){
            SVD[j][i]=SVDt[i][j];
        }
    }


  printf("\nThe corresponding U is: \n");

  for(j=0; j<m; j++) {

    for(i=0; i<m; i++){

      printf("% 8.4lf ",U[i][j]);
      }
      printf("\n");

  }
  printf("\nSingular matrix:\n");

    for(i=0;i<m;i++){
        for(j=0;j<m;j++){
            printf("%8.4lf ",S[i][j]);
        }printf("\n");
    }
  printf("\nThe corresponding v is: \n");

  for(j=0; j<m; j++) {

    for(i=0; i<m; i++){

      printf("% 8.4lf ",V[i][j]);
      }
      printf("\n");

  }
  printf("\nThe Complete SVD is: \n");

  for(j=0; j<m; j++) {

    for(i=0; i<m; i++){

      printf("% 8.4lf ",SVD[i][j]);
      }
      printf("\n");

  }

}
void multiplication(double A[100][100],double B[100][100],double C[100][100],int r,int c){
    for(int i=0;i<r;i++)
        {
        for(int j=0;j<c;j++)
        {
        C[i][j]=0;
        for(int k=0;k<c;k++)
        {
         C[i][j]+=A[i][k]*B[k][j];
        }
        }
        }

}
int forwardElim(double w3[100][100]);

// function to calculate the values of the unknowns
void backSub(double w3[100][100]);

// function to get matrix content
void gaussianElimination(double w3[100][100])
{
    /* reduction into r.e.f. */
    int singular_flag = forwardElim(w3);

    /* if matrix is singular */
    if (singular_flag != -1)
    {
        printf("Singular Matrix.\n");

        /* if the RHS of equation corresponding to
           zero row  is 0, * system has infinitely
           many solutions, else inconsistent*/
        if (w3[singular_flag][100])
            printf("Inconsistent System.");
        else
            printf("May have infinitely many "
                   "solutions.");

        return;
    }


    backSub(w3);
}


void swap_row(double w3[100][100], int i, int j)
{


    for (int k=0; k<=m; k++)
    {
        double temp = w3[i][k];
        w3[i][k] = w3[j][k];
        w3[j][k] = temp;
    }
}


// function to reduce matrix to r.e.f.
int forwardElim(double w3[100][100])
{
    for (int k=0; k<m; k++)
    {
        // Initialize maximum value and index for pivot
        int i_max = k;
        int v_max = w3[i_max][k];

        /* find greater amplitude for pivot if any */
        for (int i = k+1; i < m; i++)
            if (abs(w3[i][k]) > v_max)
                v_max = w3[i][k], i_max = i;

        /* if a principal diagonal element  is zero,
         * it denotes that matrix is singular, and
          */
        if (!w3[k][i_max])
            return k; // Matrix is singular

        /* Swap the greatest value row with current row */
        if (i_max != k)
            swap_row(w3, k, i_max);


        for (int i=k+1; i<m; i++)
        {
            /* factor f to set current row kth element to 0,
             * and subsequently remaining kth column to 0 */
            double f = w3[i][k]/w3[k][k];

            /* subtract fth multiple of corresponding kth
               row element*/
            for (int j=k+1; j<=m; j++)
                w3[i][j] -= w3[k][j]*f;

            /* filling lower triangular matrix with zeros*/
            w3[i][k] = 0;
        }


    }

    return -1;
}


void backSub(double w3[100][100])
{
    double x[100];

    for (int i = m-2; i >= 0; i--)
    {
        x[m-1]=1;

        for (int j=i+1; j<m; j++)
        {

            x[i] =x[i]- w3[i][j]*x[j];
        }


        x[i] = x[i]/w3[i][i];
    }

    //cout<<"Eigenvector for Corresponding eigenvalue :"<<endl;
    printf("\n Eigenvector for Corresponding eigen vector :\n");
    for (int i=0; i<m; i++)
        printf("%lf\n", x[i]);
    //cout<<endl;printf("\n");
    double sum=0;

    //int sum=0;
    for(int j=0; j<m; j++)
    {
        sum = sum+ pow(x[j],2);
    }
    sum = sqrt(sum);

    for(int i=0; i<m; i++)
    {
        x[i]=x[i]/sum;
        V[i][l] = x[i];
    }
    l++;

    printf("\nSolution for the system:\n");
    for (int i=0; i<m; i++)
        printf("%lf\n", x[i]);
}


