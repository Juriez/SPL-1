#include<bits/stdc++.h>
using namespace std;
#define Dimen 4
double A[Dimen + 1][Dimen + 1],AM[Dimen + 1][Dimen + 1], M[Dimen+1][Dimen+1], A_t[Dimen+1][Dimen+1], Inverse[Dimen+1][Dimen+1],w[Dimen+1][Dimen+1];
int ans = 0, h;

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
int n, k = 0;

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

                //cout<<"\tRoot: "<<"i"<<endl;

                //cout<<"\tRoot: "<<"-i"<<endl;

            }

            else

            {



                //cout<<"\tRoot: "<<(determine/(2*x))<<"i"<<endl;

                //cout<<"\tRoot: "<<(determine/(2*x))<<"-i"<<endl;

            }

        }

        else

        {

            // not a pure imaginary number

            if((determine/(2*x)) == 1)

            {

                //coefficient 1 , not necessary to print

                //cout<<"\tRoot: "<<((-p)/(2*x))<<" + "<<"i"<<endl;

                //cout<<"\tRoot: "<<((-p)/(2*x))<<" - "<<"i"<<endl;

            }

            else

            {

                //there are coefficient

                //cout<<"\tRoot: "<<((-p)/(2*x))<<" + "<<(determine/(2*x))<<"i"<<endl;//

                //cout<<"\tRoot: "<<((-p)/(2*x))<<" - "<<(determine/(2*x))<<"i"<<endl;

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



        //cout<<"\tRoot: "<<first<<endl;
        Root1=first;

        //cout<<"\tRoot: "<<second<<endl;
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

    //if(!f)
    Root = root;
    eig_val[k] = Root;
    k++;
    //f = false;
    //}
    //cout<<"\tRoot: "<<root<<endl;
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

//Finding coefficients using Faddeev-LeVerrier algorithm.

double coefficient_calcul(int k) {
    for(int i=0;i<Dimen;i++) {
        for(int j=0;j<Dimen;j++) {
            if(k - 1 == 0) h = 0;
            else  h = ans/(k - 1);
            if(i == j)M[i][j] = AM[i][j] - h;
            else M[i][j] = AM[i][j];
            //cout << M[i][j] << ' ';
            Inverse[i][j]=M[i][j];
        }
        //cout << '\n';
    }

    for(int i=0;i<Dimen;i++) {
       for(int j=0;j<Dimen;j++) {
           int sum = 0;
           for(int l=0;l<Dimen;l++){
              sum += A[i][l]*M[l][j];
           }
           AM[i][j]=sum;
           //cout << sum << ' ';
       }
       //cout << '\n';
    }
    ans = 0;

    for(int i=0, j = 0;i<Dimen;i++,j++){
        ans += AM[i][j];
    }
    return (-1)*ans/k;
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
        for (i = 0;i < k; i++)
          {
            for (j = 0 ;j < k; j++)
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

void multiplication(double A[Dimen+1][Dimen+1],double A_t[Dimen+1][Dimen+1])
{
    double sum =0;
    for(int i=0;i<Dimen;i++)
   {
       for(int j=0;j<Dimen;j++)
       {
           for(int k=0;k<Dimen;k++)
           {
               w[i][j]=A[i][k]*A_t[k][j];
               sum=sum+w[i][j];

           }
           w[i][j]=sum;
           sum=0;
       }
   }
   cout<<endl<<"So, multiplication of to matrix is :"<<endl;
   for(int i=0;i<Dimen;i++)
   {
       for(int j=0;j<Dimen;j++)
       {
           cout<<w[i][j]<<" ";
       }cout<<endl;
   }
}

void transpose(double A[Dimen+1][Dimen+1])

{
    //double trans_val;
    for(int i=0;i<Dimen;i++)
    {
        for(int j=0;j<Dimen;j++)
        {
           A_t[i][j] = A[j][i];

        }
        //cout<<endl;
    }

    cout<<endl<<"So,the transpose is :"<<endl;

    for(int i=0;i<Dimen;i++)
    {
        for(int j=0;j<Dimen;j++)
        {
            cout<<A_t[i][j]<<" ";
        }
        cout<<endl;
    }
}

int main() {
    FILE *file;
    file=fopen("spl_task2.txt","r");
    for(int i=0;i<Dimen;i++){
        for(int j=0;j<Dimen;j++) {
            fscanf(file,"%lf",&A[i][j]);
        }
    }
    cout<<endl<<"Given matrix :"<<endl;
    for(int i=0;i<Dimen;i++)
    {
        for(int j=0;j<Dimen;j++)
        {
            cout<<A[i][j]<<" ";
        }cout<<endl;
    }

    transpose(A);

    multiplication(A,A_t);

    double det = dtermant(A,Dimen);
    cout<<endl<<"Determinant of the matrix is : "<<det<<endl;

    for(int i = 0,j = 0; i < Dimen; i++,j++) AM[i][j] = 1;
    double a[Dimen + 1],coefficients[Dimen+1];
    for(int k = 1; k <= Dimen;k++)a[k] = coefficient_calcul(k);
    a[0] = 1;
    cout<<endl<<"Coefficients for the polynomial function are..."<<endl;
    for(int k = 0; k <= Dimen;k++) cout << a[k] << " ";
    //cout<<endl<<"Determinant of the matrix : "<< (-1)*a[Dimen]<<endl;

    cout<<endl<<"Solving polynomial equation :"<<endl;


    for(int i=0;i<=Dimen;i++){

        coefficients[Dimen-i]=a[i];
        //cin>>a[n-i];

    }
    //cout<<coefficients[0];
    for(int i=0;i<Dimen;i++)
    {
        for(int j=0;j<Dimen;j++)
        {
            Inverse[i][j]=(-1 )* Inverse[i][j]/coefficients[0];
        }
        //cout<<endl;
    }
    cout<<endl<<"So the Inverse of the given matrix is :"<<endl;
    for(int i=0;i<Dimen;i++)
    {
        for(int j=0;j<Dimen;j++)
        {
           cout<<Inverse[i][j]<<" ";
        }
        cout<<endl;
    }


    Get_Coefficient(coefficients,Dimen);
    cout<<endl<<"Eigen values are :"<<endl;

    for(int i=0;i<Dimen;i++)
    {
        if(eig_val[i]) cout<<eig_val[i] << ' ';
    }

}
/*
90 23 67 87 89
78 23 10 45 67
19 56 87 98 20
78 45 23 91 96
69 73 86 0 92
*/
/*
-3 2 0 0
-3 4 0 0
0 0 -5 -4
0 0 -2 2
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
