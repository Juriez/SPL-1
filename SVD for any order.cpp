#include<bits/stdc++.h>
using namespace std;
#define Dimen 150
double A[Dimen + 1][Dimen + 1],WM[Dimen + 1][Dimen + 1], M[Dimen+1][Dimen+1]={0}, A_t[Dimen+1][Dimen+1],w[Dimen+1][Dimen+1],U[Dimen+1][Dimen+1],Ut[Dimen+1][Dimen+1],SVDt[Dimen+1][Dimen+1], SiUt[Dimen+1][Dimen+1], PI[Dimen+1][Dimen+1],Vt[Dimen+1][Dimen+1],V[Dimen+1][Dimen+1],S[Dimen+1][Dimen+1],SI[Dimen+1][Dimen+1],SV[Dimen+1][Dimen+1],SVD[Dimen+1][Dimen+1], VtSi[Dimen+1][Dimen+1];
double ans = 0, h;
int k = 0;
double coefficients[Dimen+1];
void multiplication(double p[Dimen+1][Dimen+1],double q[Dimen+1][Dimen+1],double r[Dimen+1][Dimen+1]);
void matrixProjection(double P[Dimen+1][Dimen+1],double Q[Dimen+1][Dimen+1],double R[Dimen+1][Dimen+1],int c);
#define e 2.72828

#define sigma 4.0

#define mu 10.0

#define initial (0.398922804/sigma)

#define N 100

#define SIZE 100000

#define phi 1e-12

double Root;
double Root1;
double Root2;
double eig_val[Dimen];
int n, l = 0;

double *a,*b,*c, r, s, old_r, old_s,dr,ds,root_p, root_q;

bool last = false, f = false;
//Finding Eigen value using Bairstow Algorithm.

void new_arr()

{

    a = new double[N];

    b = new double[N];

    c = new double[N];

}



void del_arr()

{

    delete a;

    delete b;

    delete c;



}

double absolute(double x)

{
    //return abs value

    if(x<0)

    {

        x *= -1;

    }

    return x;

}


double remove_eror(double val)

{

    //floating point eror

    int integer = val;

    if(absolute(integer - val) <= phi)

    {

        val = (double)integer;

    }



    return val;

}


void printRoot(double x,double p, double q)

{



    if(!last)

    {

        //for r and s , we need to change sign of r and s

        p = (-1)*remove_eror(p);

        q =(-1)*remove_eror(q);

    }

    else

    {

        // we find solution from rest of equation

        p = remove_eror(p);

        q = remove_eror(q);

    }



    double determine = (p*p) - (4*x*q);





    if(determine<0)

    {

        determine *= -1;



        determine = sqrt(determine);



        if(p == 0)

        {

            //pure imaginary number

            if((determine/(2*x)) == 1)

            {

                //coefficient 1

                cout<<"\tRoot: "<<"i"<<endl;

                cout<<"\tRoot: "<<"-i"<<endl;

            }

            else

            {



                cout<<"\tRoot: "<<(determine/(2*x))<<"i"<<endl;

                cout<<"\tRoot: "<<(determine/(2*x))<<"-i"<<endl;

            }

        }

        else

        {

            // not a pure imaginary number

            if((determine/(2*x)) == 1)

            {

                //coefficient 1 , not necessary to print

                cout<<"\tRoot: "<<((-p)/(2*x))<<" + "<<"i"<<endl;

                cout<<"\tRoot: "<<((-p)/(2*x))<<" - "<<"i"<<endl;

            }

            else

            {

                //there are coefficient

                cout<<"\tRoot: "<<((-p)/(2*x))<<" + "<<(determine/(2*x))<<"i"<<endl;//

                cout<<"\tRoot: "<<((-p)/(2*x))<<" - "<<(determine/(2*x))<<"i"<<endl;

            }

        }

    }

    else

    {

        //cout<<determine<<endl;

        determine = sqrt(determine);

        //cout<<determine<<endl;

        double first = remove_eror(((-p) - determine)/(2*x));

        double second =  remove_eror(((-p) + determine)/(2*x));



        cout<<"\tRoot: "<<first<<endl;
        Root1=first;

        cout<<"\tRoot: "<<second<<endl;
        Root2=second;

        eig_val[k] = Root1;
        eig_val[k + 1] = Root2;
        k += 2;
        f = true;
    }

}

void print_rootOne(double x, double y)

{

    // If existance equation has only one solution

    //x *= -1;

    //y *= -1;
    //double roots[Dimen];
    int i=0;



    double root = -(y/x);

    //if(!f){
    Root = root;
    eig_val[k] = Root;
    k++;
   // f = false;
    //}
    cout<<"\tRoot: "<<root<<endl;
    //roots[i]=root;
    //i++;

}




void cal_r_s()

{

    //for iteration we need to find r and s



    dr = (b[0]*c[3] - b[1]*c[2]) / (c[2]*c[2] - c[1]*c[3]);

    ds = (b[1]*c[1] - b[0]*c[2]) / (c[2]*c[2] - c[1]*c[3]);



    old_r = r ;

    old_s = s;

    r += dr;

    s += ds;

}

void cal_col(double p[], double q[])

{

    //iterate column to find solution

    q[n] = p[n];

    q[n-1] = p[n-1] + q[n]*r;



    for(int i=n-2; i>=0; i--)

    {

        q[i] = p[i] + (q[i+1]*r) + (q[i+2]*s);

    }

}

void reduce_eqn()

{

    //After iteration, found two solution and the equations's power reduce by two.

    //Replace a[] by b[]

    for(int i=0; i<n-1; i++)

    {

        a[i] = b[i+2];

    }



    n -= 2;

}





void RootFinding()

{

    //iteration until equation goes to power 1 or 2

    double ratio_s, ratio_r;


    if(n == 0)

    {

        cout<<"No such variable.\n Wrong input\n\n";

        exit(0);

    }

    else if(n == 1)

    {

        print_rootOne(a[n], a[n-1]);

    }

    else if(n == 2)

    {

        last = true;

        printRoot(a[n], a[n-1], a[n-2]);

    }

    else

    {

        while(1)

        {

            cal_col(a,b);

            cal_col(b,c);

            cal_r_s();

            ratio_s = ds/old_s;

            ratio_r = dr/old_r;


            if(((absolute(b[0]) <= phi) && (absolute(b[1]) <= phi)) || ((absolute(ratio_r) <= phi) || (absolute(ratio_s) <= phi)))

            {

                printRoot(1,r,s);

                if(n == 4)

                {

                    last = true;

                    printRoot(b[n],b[n-1],b[n-2]);

                    break;

                }
                if(n == 3)

                {

                    print_rootOne(b[n], b[n-1]);

                    break;

                }



                reduce_eqn();

            }

        }

    }

}

void Get_Coefficient(double pass[Dimen+1], int total)

{

    new_arr();

    for(int i= total; i>=0; i--)

    {

        a[i] = pass[i];

    }

    n = total;





    r = a[n-1]/a[n];

    s = a[n-2]/a[n];



    if(r == 0)

    {

        r = 0.1;

    }



    if(s ==0)

    {

        s = 0.1;

    }



    RootFinding();



    del_arr();         //deallocate memory

    return;

}


// function to reduce matrix to r.e.f.  Returns a value to
// indicate whether matrix is singular or not
int forwardElim(double w3[Dimen+1][Dimen+1]);

// function to calculate the values of the unknowns
void backSub(double w3[Dimen+1][Dimen+1]);

// function to get matrix content
void gaussianElimination(double w3[Dimen+1][Dimen+1])
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
        if (w3[singular_flag][Dimen])
            printf("Inconsistent System.");
        else
            printf("May have infinitely many "
                   "solutions.");

        return;
    }


    backSub(w3);
}


void swap_row(double w3[Dimen+1][Dimen+1], int i, int j)
{


    for (int p=0; p<=Dimen; p++)
    {
        double temp = w3[i][p];
        w3[i][p] = w3[j][p];
        w3[j][p] = temp;
    }
}


// function to reduce matrix to r.e.f.
int forwardElim(double w3[Dimen+1][Dimen+1])
{
    for (int p=0; p<Dimen; p++)
    {
        // Initialize maximum value and index for pivot
        int i_max = p;
        int v_max = w3[i_max][p];

        /* find greater amplitude for pivot if any */
        for (int i = p+1; i < Dimen; i++)
            if (abs(w3[i][p]) > v_max)
                v_max = w3[i][p], i_max = i;

        /* if a principal diagonal element  is zero,
         * it denotes that matrix is singular, and
          */
        if (!w3[p][i_max])
            return p; // Matrix is singular

        /* Swap the greatest value row with current row */
        if (i_max != p)
            swap_row(w3, p, i_max);


        for (int i=p+1; i<Dimen; i++)
        {
            /* factor f to set current row kth element to 0,
             * and subsequently remaining kth column to 0 */
            double f = w3[i][p]/w3[p][p];

            /* subtract fth multiple of corresponding kth
               row element*/
            for (int j=p+1; j<=Dimen; j++)
                w3[i][j] -= w3[p][j]*f;

            /* filling lower triangular matrix with zeros*/
            w3[i][p] = 0;
        }


    }

    return -1;
}


void backSub(double w3[Dimen+1][Dimen+1])
{
    double x[Dimen];

    for (int i = Dimen-2; i >= 0; i--)
    {
        x[Dimen-1]=1;

        for (int j=i+1; j<Dimen; j++)
        {

            x[i] =x[i]- w3[i][j]*x[j];
        }


        x[i] = x[i]/w3[i][i];
    }

    cout<<"Eigenvector for Corresponding eigenvalue :"<<endl;
    for (int i=0; i<Dimen; i++)
        printf("%lf\n", x[i]);
    cout<<endl;
    double sum=0;

    //int sum=0;
    for(int j=0; j<Dimen; j++)
    {
        sum = sum+ pow(x[j],2);
    }
    sum = sqrt(sum);

    for(int i=0; i<Dimen; i++)
    {
        x[i]=x[i]/sum;

        Vt[i][l] = x[i];
        //V[l][i] = Vt[i][l];


    }
    l++;

    printf("\nSolution for the system:\n");
    for (int i=0; i<Dimen; i++)
        printf("%lf\n", x[i]);
}


//Finding coefficients using Faddeev-LeVerrier algorithm.

void coefficient_calcul()
{
     for(int j=1;j<=Dimen;j++)
    {
        for(int i=0;i<Dimen;i++){
        for(int p=0;p<Dimen;p++){
            WM[i][p]=0;
        }
      }
        double trace =0;
        multiplication(w,M,WM);
        for(int i=0;i<Dimen;i++){
            trace += WM[i][i];
        }
        trace = ((-1)*trace)/j;
        coefficients[j]=trace;
        //trace = 0;
        for(int i=0;i<Dimen;i++){
            //
            for(int p=0;p<Dimen;p++){
                if(i == p)
                    {M[i][p] = WM[i][p] + coefficients[j];}
                else
                    {M[i][p] = WM[i][p];}
            }
        }



    }
   /* for(int i=0;i<Dimen;i++){
        for(int j=0;j<Dimen;j++){
            cout<<WM[i][j]<<" ";
        }cout<<endl;
    }*/
}

double dtermant(double A[Dimen+1][Dimen+1],double k)
{
    double s = 1,  b[Dimen+1][Dimen+1];
    double det=0;
    int i, j, m, n, c;
    if (k == 1)
    {
        return (A[0][0]);
    }
    else
    {
        det = 0;
        for (c = 0; c < k; c++)
        {
            m = 0;
            n = 0;
            for (i = 0; i < k; i++)
            {
                for (j = 0 ; j < k; j++)
                {
                    b[i][j] = 0;
                    if (i != 0 && j != c)
                    {
                        b[m][n] = A[i][j];
                        if (n < (k - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det = det + s * (A[0][c] * dtermant(b, k - 1));
            s = -1 * s;
        }

    }
    //cout<<endl<<"Determinant of the matrix is : "<<det<<endl;
    return det;

}

void multiplication(double p[Dimen+1][Dimen+1],double q[Dimen+1][Dimen+1],double r[Dimen+1][Dimen+1])
{
    double sum =0;
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            for(int k=0; k<Dimen; k++)
            {
                r[i][j]=p[i][k]*q[k][j];
                sum=sum+r[i][j];

            }
            r[i][j]=sum;
            sum=0;
        }
    }

}

void transpose(double P[Dimen+1][Dimen+1],double Q[Dimen+1][Dimen+1])

{
    //double trans_val;
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            Q[i][j] = P[j][i];

        }
        //cout<<endl;
    }

    /*cout<<endl<<"So,the transpose is :"<<endl;

    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<Q[i][j]<<" ";
        }
        cout<<endl;
    }*/
}

void SVDforHigherOrder()
{
      int  i, j, p, q, flag;
    double d[Dimen+1][Dimen+1], s[Dimen+1][Dimen+1], s1[Dimen+1][Dimen+1], s1t[Dimen+1][Dimen+1];
    //double W[Dimen+1][Dimen+1];
    double temp[Dimen+1][Dimen+1], theta, zero=1e-5, max, pi=3.141592654;

  for(i=0; i<Dimen; i++) {
    for(j=0; j<Dimen; j++) {
      d[i][j]= w[i][j];
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

    for(i=0; i<Dimen; i++) {
      for(j=0; j<Dimen; j++) {
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

    for(i=0; i<Dimen; i++) {
      for(j=0; j<Dimen; j++) {
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

    for(i=0; i<Dimen; i++) {
      for(j=0; j<Dimen; j++) {
        temp[i][j]= 0;
        for(p=0; p<Dimen; p++)  temp[i][j]+= s1t[i][p]*d[p][j];
      }
    }

    for(i=0; i<Dimen; i++) {
      for(j=0; j<Dimen; j++) {
        d[i][j]= 0;
        for(p=0; p<Dimen; p++)  d[i][j]+= temp[i][p]*s1[p][j];
      }
    }

    for(i=0; i<Dimen; i++) {
      for(j=0; j<Dimen; j++) {
        temp[i][j]= 0;
        for(p=0; p<Dimen; p++)  temp[i][j]+= s[i][p]*s1[p][j];
      }
    }

    for(i=0; i<Dimen; i++) {
      for(j=0; j<Dimen; j++)  s[i][j]= temp[i][j];
    }

    for(i=0; i<Dimen; i++) {
      for(j=0; j<Dimen; j++) {
        if(i!=j)
          if(fabs(d[i][j]) > zero) flag= 1;
      }
    }
  } while(flag==1);

  printf("\nThe eigenvalues are \n");
  for(i=0; i<Dimen; i++){
    eig_val[i]=d[i][i];
    printf("%8.4lf ",d[i][i]);

  }
    for(i=0;i<Dimen;i++){
        for(j=i+1;j<Dimen;j++){
            if(eig_val[i]<eig_val[j]){
                double temp =0;
                temp = eig_val[i];
                eig_val[i]=eig_val[j];
                eig_val[j]=temp;
            }
        }
    }

    printf("\n Sorted eigen values are :\n");
    for(i=0; i<Dimen; i++){
    printf("%8.4lf ",eig_val[i]);}
    for(i=0;i<Dimen;i++){
        S[i][i] = sqrt(eig_val[i]);
        SI[i][i] = 1/S[i][i];
    }
   // double W[Dimen+1][Dimen+1];
    //multiplication(A_t,A,W);
     double w2[Dimen][Dimen+1];

    int k = 0;

    while(k < Dimen) {
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            w2[i][j]=w[i][j];
        }
    }

    for(int j=0; j<Dimen; j++)
    {
        w2[j][j] = w[j][j]-eig_val[k];
    }
    /*for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<w2[i][j]<<" ";
        }cout<<endl;
    }*/
    cout<<endl<<endl;
    gaussianElimination(w2);
    k++;
    }

    /*for(int i=0;i<Dimen;i++){
        for(int j=0;j<Dimen;j++){
            V[i][j]=Vt[j][i];
        }
    }*/
    transpose(Vt,V);
    cout<<"Vt :"<<endl;

    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<Vt[i][j]<<" ";
        }cout<<endl;
    }


  multiplication(Vt,SI,VtSi);
  multiplication(A,VtSi,U);
  /*for(i=0;i<Dimen;i++){
        for(j=0;j<Dimen;j++){
            Ut[j][i]=U[i][j];
        }
    }*/
    transpose(U,Ut);
     /*cout<<"Ut :"<<endl;
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<Ut[i][j]<<" ";
        }cout<<endl;
    }*/
  multiplication(S,V,SV);
  multiplication(U,SV,SVD);
  /*for(i=0;i<Dimen;i++){
        for(j=0;j<Dimen;j++){
            SVD[j][i]=SVDt[i][j];
        }
    }*/
    transpose(SVD,SVDt);


  printf("\nThe corresponding U is: \n");

  for(j=0; j<Dimen; j++) {

    for(i=0; i<Dimen; i++){

      printf("% 8.4lf ",U[j][i]);
      }
      printf("\n");

  }
  printf("\nSingular matrix:\n");

    for(i=0;i<Dimen;i++){
        for(j=0;j<Dimen;j++){
            printf("%8.4lf ",S[i][j]);
        }printf("\n");
    }
  printf("\nThe corresponding v is: \n");

  for(j=0; j<Dimen; j++) {

    for(i=0; i<Dimen; i++){

      printf("% 8.4lf ",V[j][i]);
      }
      printf("\n");

  }
  printf("\nThe Complete SVD is: \n");

  for(j=0; j<Dimen; j++) {

    for(i=0; i<Dimen; i++){

      printf("% 8.4lf ",SVD[j][i]);
      }
      printf("\n");

  }
   cout<<endl<<"Error between main matrix & SVD is : ";
    double sum = 0;
    for(int i=0;i<Dimen;i++)
    {
        for(int j=0;j<Dimen;j++)
        {
            sum = sum+ pow((A[i][j] - SVD[i][j]),2);
        }
    }
    sum = sqrt(sum);
    cout<< sum<<endl;

    // calculating Pseudo inverse
    multiplication(SI,Ut,SiUt);
    multiplication(Vt,SiUt,PI);

    cout<<endl<<"Pseudo Inverse of A is :"<<endl;
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<PI[i][j]<<" ";
        }
        cout<<endl;
    }
}
void LeverageScoreSampling(){
        double D[Dimen+1][Dimen+1]= {0};
        double lvs[Dimen],total_sum = 0, rowwise_sum = 0;
        int indx_of_high_lvs[Dimen+1];
        double arr[Dimen+1];
        int d = Dimen-1, count = 0;
        for(int i=0;i<Dimen;i++){
            D[i][i] = eig_val[i]/d;
        }


        cout<<"D :"<<endl<<endl;
        for(int i=0;i<Dimen;i++){
            for(int j=0;j<Dimen;j++){
                cout<<D[i][j]<<" ";
            }cout<<endl;
        }


        for(int i=0;i<Dimen;i++){
                double sum = 0;
            for(int j=0;j<Dimen;j++){
                sum = sum + pow(U[i][j],2);
            }
             lvs[i] = sum;

             lvs[i] = lvs[i] / D[i][i];
             arr[i] = lvs[i];
             total_sum =total_sum+lvs[i];
             sum = 0;
        }

        /*cout<<"Leverage Scores :"<<endl;
        for(int i=0 ;i<Dimen;i++){
            cout<<arr[i] << "  ";
        }cout<<endl;*/

        for(int i=0;i<Dimen;i++){
            for(int j=i+1;j<Dimen;j++){
                double temp = 0;
                if(lvs[i]<lvs[j]){
                    temp = lvs[j];
                    lvs[j]=lvs[i];
                    lvs[i]=temp;
                }
            }
        }


        cout<<endl<<"Sorted Leverage Scores :"<<endl;
        for(int i=0 ;i<Dimen;i++){
            cout<<lvs[i] << "  ";
        }cout<<endl;

        cout<<endl<<"Total_sum : "<<total_sum<<endl;
        for(int i=0;i<Dimen;i++){
            rowwise_sum = rowwise_sum + lvs[i];
            if(rowwise_sum/total_sum >= 0.95){
                count++;
                break;
            }
            else
                count++;
        }
        cout<<endl<<"Count : "<<count<<endl;
        for(int i=0;i<count;i++){
            for(int j=0;j<Dimen;j++){
                if(lvs[i] == arr[j])
                    indx_of_high_lvs[i] = j;
            }
        }
        cout<<endl<<"Index of High Leverage scores : "<<endl;
        for(int i=0;i<count;i++){
            cout<<indx_of_high_lvs[i]<<"  ";
        }
        cout<<endl;

        double A_k[Dimen+1][Dimen+1],U_k[Dimen+1][Dimen+1],S_k[Dimen+1][Dimen+1] = {0},Vt_k[Dimen+1][Dimen+1],V_k[Dimen+1][Dimen+1];
        double SkVkt[Dimen+1][Dimen+1],sampl[Dimen+1][Dimen+1];
        int ind = count-1, ind2 = count-1,ind3 = 0;
        for(int i=0;i<count;i++){
            for(int j=0;j<Dimen;j++){
                U_k[j][i] = U[indx_of_high_lvs[ind]][j];
                V_k[j][i] = V[indx_of_high_lvs[ind]][j];
                //S_k[i][j] = S[indx_of_high_lvs[ind]][j];
                sampl[i][j] = A[indx_of_high_lvs[ind3]][j];

            }
            ind--;
            ind3++;
        }
        transpose(V_k,Vt_k);
        for(int i=0;i<count;i++){
            S_k[i][i] = S[indx_of_high_lvs[ind2]][indx_of_high_lvs[ind2]];
            ind2--;
        }

        cout<<endl<<"U_k :"<<endl;
        for(int i=0;i<Dimen;i++){
            for(int j=0;j<count;j++){
                cout<<U_k[i][j]<<" ";
            }cout<<endl;
        }
        cout<<endl<<"S_k :"<<endl;
        for(int i=0;i<count;i++){
            for(int j=0;j<count;j++){
                cout<<S_k[i][j]<<" ";
            }cout<<endl;
        }
        cout<<endl<<"Vt_k :"<<endl;
        for(int i=0;i<count;i++){
            for(int j=0;j<Dimen;j++){
                cout<<Vt_k[i][j]<<" ";
            }cout<<endl;
        }


        int g=0;
        for(int i=0;i<count;i++){
            for(int j=0;j<Dimen;j++){
                SkVkt[i][j] = S_k[i][i]*Vt_k[i][j];
            }
        }

        /*cout<<"SKVKT :"<<endl;
        for(int i=0;i<count;i++){
            for(int j=0;j<Dimen;j++){
                cout<<SkVkt[i][j]<<" ";
            }cout<<endl;
        }*/
        matrixProjection(U_k,SkVkt,A_k,count);
        cout<<endl<<"Reconstracted A_K :"<<endl;
        for(int i=0;i<Dimen;i++){
            for(int j=0;j<Dimen;j++){
                cout<<A_k[i][j]<<" ";
            }cout<<endl;
        }
        //int ind3 = count-1;
        cout<<endl<<"Sample data Matrix :"<<endl;
        for(int i=0;i<count;i++){
            for(int j=0;j<Dimen;j++){
                cout<<sampl[i][j]<<" ";
            }
            cout<<endl;
        }





}
void matrixProjection(double P[Dimen+1][Dimen+1],double Q[Dimen+1][Dimen+1],double R[Dimen+1][Dimen+1],int c){
    for(int i=0;i<Dimen;i++)
        {
        for(int j=0;j<Dimen;j++)
        {
        R[i][j]=0;
        for(int g=0;g<c;g++)
        {
         R[i][j]+=P[i][g]*Q[g][j];
        }
        }
        }

}


int main()
{
    int maxn=200,minn=0;
    int ran,number;
    for(int i=0;i<1000;i++){
        for(int j=0;j<1000;j++){
    ran = rand();
    number=ran % (maxn-minn)+minn ;
    A[i][j]=number;
            }
    }
    FILE *fp;
    fp = fopen("spl_task2.txt","w");
    for(int i=0;i<Dimen;i++){
    for(int j=0;j<Dimen;j++){
        fprintf(fp,"%lf ",A[i][j]);
    }fprintf(fp,"\n");
    }
    /*FILE *file;
    file=fopen("spl_task2.txt","r");
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            fscanf(file,"%lf",&A[i][j]);
        }
    }*/

    //cout<<endl<<"kdbj"<<coun<<endl;
    cout<<endl<<"Given matrix :"<<endl;
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }

    transpose(A,A_t);
    cout<<endl<<"So,the transpose is :"<<endl;

    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<A_t[i][j]<<" ";
        }
        cout<<endl;
    }

    multiplication(A_t,A,w);
    cout<<endl<<"So, multiplication between A_t & A is :"<<endl;

    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<w[i][j]<<" ";
        }
        cout<<endl;
    }

    //double det = dtermant(A,Dimen);
    //cout<<endl<<"Determinant of the matrix is : "<<det<<endl;

    for(int i=0;i<Dimen;i++)
    {
        M[i][i]=1;
    }
    if(Dimen>10){
        SVDforHigherOrder();
    }
    else
    {
    coefficient_calcul();
    cout<<"Coefficients of the characteristic equation are :"<<endl;
    coefficients[0]=1;
    for(int i=0;i<=Dimen;i++){
        cout<<coefficients[i]<<"  ";
    }
    double co[Dimen+1];
    for(int i=0; i<=Dimen; i++)
    {

        co[Dimen-i]=coefficients[i];


    }


    cout<<endl<<"Solving polynomial equation :"<<endl;


    Get_Coefficient(co, Dimen);
    cout<<endl<<"Eigen values are :"<<endl;

    for(int i=0; i<Dimen; i++)
    {
        if(eig_val[i]) cout<<eig_val[i] << ' ';
    }

    double temp =0 ;

    for(int i=0;i<Dimen;i++)
    {
        for(int j=i+1;j<Dimen;j++)
        {
            if(eig_val[i]<eig_val[j])
            {
                temp = eig_val[i];
                eig_val[i]=eig_val[j];
                eig_val[j]=temp;
            }
        }
    }

    cout<<endl<<"Sorted Eigen values are :"<<endl;

    for(int i=0; i<Dimen; i++)
    {
        if(eig_val[i]) cout<<eig_val[i] << ' ';
    }
    cout<<endl;

    for(int i=0;i<Dimen;i++)
    {
        S[i][i]=sqrt(eig_val[i]);
        SI[i][i]= 1/S[i][i];
    }

    /*cout<<endl<<"Singular matrix Inverse:"<<endl;
    for(int i=0;i<Dimen;i++)
    {
        for(int j=0;j<Dimen;j++)
        {
            cout<<SI[i][j]<<" ";
        }cout<<endl;
    }*/

    double w2[Dimen][Dimen+1];

    int k = 0;

    while(k < Dimen) {
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            w2[i][j]=w[i][j];
        }
    }

    for(int j=0; j<Dimen; j++)
    {
        w2[j][j] = w[j][j]-eig_val[k];
    }
    gaussianElimination(w2);
    k++;
    }
    cout<<"VT :"<<endl;
    for(int i=0;i<Dimen;i++)
    {
        for(int j=0;j<Dimen;j++)
        {
            cout<<Vt[i][j]<<" ";
        }cout<<endl;
    }

    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            V[i][j] = Vt[j][i];

        }

    }

    multiplication(Vt,SI,VtSi);  // start to calculate U ;
    multiplication(A,VtSi,U);  // end to calculate U ;

    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            Ut[i][j] = U[j][i];

        }

    }

    cout<<endl<<"So, the Left Singular Vector U is :"<<endl;

    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<U[i][j]<<" ";
        }
        cout<<endl;
    }

    cout<<endl<<"Singular matrix S:"<<endl;

    for(int i=0;i<Dimen;i++)
    {
        for(int j=0;j<Dimen;j++)
        {
            cout<<S[i][j]<<" ";
        }cout<<endl;
    }

    cout<<endl;
    cout<<"Right Singular Vector V is :"<<endl;

    for(int i=0;i<Dimen;i++)
    {
        for(int j=0;j<Dimen;j++)
        {
            cout<<V[i][j]<<" ";
        }cout<<endl;
    }
    cout<<endl;

    multiplication(S,V,SV); // calculating svd
    multiplication(U,SV,SVD);

    cout<<"So, the complete SVD is :"<<endl;
    for(int i=0;i<Dimen;i++)
    {
        for(int j=0;j<Dimen;j++)
        {
            cout<<SVD[i][j]<<" ";
        }cout<<endl;
    }


    cout<<endl<<"Error between main matrix & SVD is : ";
    double sum = 0;
    for(int i=0;i<Dimen;i++)
    {
        for(int j=0;j<Dimen;j++)
        {
            sum = sum+ pow((A[i][j] - SVD[i][j]),2);
        }
    }
    sum = sqrt(sum);
    cout<< sum<<endl;

    // calculating Pseudo inverse
    multiplication(SI,Ut,SiUt);
    multiplication(Vt,SiUt,PI);

    cout<<endl<<"Pseudo Inverse of A is :"<<endl;
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<PI[i][j]<<" ";
        }
        cout<<endl;
    }
    }
    LeverageScoreSampling();


}
/*
-9 7 10 -3 2
7 6 5 4 1
10 5 7 2 0
-3 4 2 1 3
2 1 0 3 7
*/
/*
1 5 -7 3 0
5 6 2 4 -3
-7 2 0 1 9
3 4 1 5 4
0 -3 9 4 2

*/
/*
1 3 -6
0 2 2
0 -1 5
*/
/*
3 1 5
3 3 1
4 6 4
*/
/*
 8 11 2 8
 0 -7 2 -1
-3 -7 2 1
 1 1 2 4
 */
 /*
 2 5 -3 3 1
5 9 4 -1 2
-3 4 7 1 6
3 -1 1 -3 2
1 2 6 2 -5
*/
/*
2 5 -3 3 1 9 8
5 9 4 -1 2 10 2
-3 4 7 1 6 0 -9
3 -1 1 -3 2 11 6
1 2 6 2 -5 9 -7
9 10 0 11 9 -7 4
8 2 -9 6 -7 4 7
*/
/*
3 -10 0 11 9 5 7 -3 0 1
8 -23 0 11 4 9 6 -7 0 4
2 -11 0 3 -6 8 9 8 3 7
0 -4 11 5 2 9 10 8 3 2
1 5 -2 7 9 8 -2 0 8 -11
-21 0 -9 4 3 8 9 0 -2 -3
-3 2 -6 29 8 -4 0 -9 3 -2
-3 0 1 10 9 -4 -6 -7 4 6
10 3 6 8 9 -4 3 0 -3 1
-9 0 1 2 11 9 -8 6 5 -3
*/
/*
3 -10 0 11 9 5 7 -3 0 1 -9 11
8 -23 0 11 4 9 6 -7 0 4 0 3
2 -11 0 3 -6 8 9 8 3 7 5 -6
0 -4 11 5 2 9 10 8 3 2 3 -1
1 5 -2 7 9 8 -2 0 8 -11 2 0
-21 0 -9 4 3 8 9 0 -2 -3 -8 7
-3 2 -6 29 8 -4 0 -9 3 -2 10 11
-3 0 1 10 9 -4 -6 -7 4 6 -1 3
10 3 6 8 9 -4 3 0 -3 1 -11 3
-9 0 1 2 11 9 -8 6 5 -3 -5 8
10 2 -8 0 -9 3 8 11 -3 7 -21 2
-15 0 1 9 -1 8 9 21 -4 2 7 -3
*/
/*
3 -10 0 11 9 5 7 -3 0 1 -9 11 -21 19 10
8 -23 0 11 4 9 6 -7 0 4 0 3 0 -22 -11
2 -11 0 3 -6 8 9 8 3 7 5 -6 18 -4 4
0 -4 11 5 2 9 10 8 3 2 3 -1 -2 8 -5
1 5 -2 7 9 8 -2 0 8 -11 2 0 -26 -1 0
-21 0 -9 4 3 8 9 0 -2 -3 -8 7 9 8 -2
-3 2 -6 29 8 -4 0 -9 3 -2 10 11 2 1 11
-3 0 1 10 9 -4 -6 -7 4 6 -1 3 0 -5 7
10 3 6 8 9 -4 3 0 -3 1 -11 3 11 3 -9
-9 0 1 2 11 9 -8 6 5 -3 -5 8 7 0 2
10 2 -8 0 -9 3 8 11 -3 7 -21 2 -5 3 1
-15 0 1 9 -1 8 9 21 -4 2 7 -3 -10 9 -5
-11 3 20 21 0 9 -4 7 6 2 -1 -21 0 16 10
-15 2 8 5 -2 8 6 11 -15 3 16 1 0 -7 20
10 -21 0 -1 2 9 8 6 11 -2 0 -1 13 14 -15
*/


