#include<bits/stdc++.h>
using namespace std;
#define Dimen 5
double A[Dimen + 1][Dimen + 1],WM[Dimen + 1][Dimen + 1], M[Dimen+1][Dimen+1], A_t[Dimen+1][Dimen+1], Inverse[Dimen+1][Dimen+1],w[Dimen+1][Dimen+1],U[Dimen+1][Dimen+1],Ut[Dimen+1][Dimen+1], SiUt[Dimen+1][Dimen+1], PI[Dimen+1][Dimen+1],Vt[Dimen+1][Dimen+1],V[Dimen+1][Dimen+1],S[Dimen+1][Dimen+1],SI[Dimen+1][Dimen+1],SV[Dimen+1][Dimen+1],SVD[Dimen+1][Dimen+1], VtSi[Dimen+1][Dimen+1];
int ans = 0, h, k = 0;

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


// function to reduce matrix to r.e.f.  Returns a value to
// indicate whether matrix is singular or not
int forwardElim(double w3[Dimen][Dimen+1]);

// function to calculate the values of the unknowns
void backSub(double w3[Dimen][Dimen+1]);

// function to get matrix content
void gaussianElimination(double w3[Dimen][Dimen+1])
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


void swap_row(double w3[Dimen][Dimen+1], int i, int j)
{


    for (int k=0; k<=Dimen; k++)
    {
        double temp = w3[i][k];
        w3[i][k] = w3[j][k];
        w3[j][k] = temp;
    }
}


// function to reduce matrix to r.e.f.
int forwardElim(double w3[Dimen][Dimen+1])
{
    for (int k=0; k<Dimen; k++)
    {
        // Initialize maximum value and index for pivot
        int i_max = k;
        int v_max = w3[i_max][k];

        /* find greater amplitude for pivot if any */
        for (int i = k+1; i < Dimen; i++)
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


        for (int i=k+1; i<Dimen; i++)
        {
            /* factor f to set current row kth element to 0,
             * and subsequently remaining kth column to 0 */
            double f = w3[i][k]/w3[k][k];

            /* subtract fth multiple of corresponding kth
               row element*/
            for (int j=k+1; j<=Dimen; j++)
                w3[i][j] -= w3[k][j]*f;

            /* filling lower triangular matrix with zeros*/
            w3[i][k] = 0;
        }


    }

    return -1;
}


void backSub(double w3[Dimen][Dimen+1])
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
    }
    l++;

    printf("\nSolution for the system:\n");
    for (int i=0; i<Dimen; i++)
        printf("%lf\n", x[i]);
}


//Finding coefficients using Faddeev-LeVerrier algorithm.

double coefficient_calcul(int k)
{
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            if(k - 1 == 0) h = 0;
            else  h = ans/(k - 1);
            if(i == j)M[i][j] = WM[i][j] - h;
            else M[i][j] = WM[i][j];
            //cout << M[i][j] << ' ';
            Inverse[i][j]=M[i][j];
        }
        //cout << '\n';
    }

    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            int sum = 0;
            for(int l=0; l<Dimen; l++)
            {
                sum += w[i][l]*M[l][j];
            }
            WM[i][j]=sum;
            //cout << sum << ' ';
        }
        //cout << '\n';
    }
    ans = 0;

    for(int i=0, j = 0; i<Dimen; i++,j++)
    {
        ans += WM[i][j];
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

void multiplication(double a[Dimen+1][Dimen+1],double b[Dimen+1][Dimen+1],double c[Dimen+1][Dimen+1])
{
    double sum =0;
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            for(int k=0; k<Dimen; k++)
            {
                c[i][j]=a[i][k]*b[k][j];
                sum=sum+c[i][j];

            }
            c[i][j]=sum;
            sum=0;
        }
    }

}

void transpose(double A[Dimen+1][Dimen+1])

{
    //double trans_val;
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            A_t[i][j] = A[j][i];

        }
        //cout<<endl;
    }

    cout<<endl<<"So,the transpose is :"<<endl;

    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<A_t[i][j]<<" ";
        }
        cout<<endl;
    }
}

int main()
{
    FILE *file;
    file=fopen("spl_task2.txt","r");
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            fscanf(file,"%lf",&A[i][j]);
        }
    }
    cout<<endl<<"Given matrix :"<<endl;
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }

    transpose(A);

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

    double det = dtermant(A,Dimen);
    cout<<endl<<"Determinant of the matrix is : "<<det<<endl;

    for(int i = 0,j = 0; i < Dimen; i++,j++) WM[i][j] = 1;
    double a[Dimen + 1],coefficients[Dimen+1];
    for(int k = 1; k <= Dimen; k++)a[k] = coefficient_calcul(k);
    a[0] = 1;
    cout<<endl<<"Coefficients for the polynomial function are..."<<endl;
    for(int k = 0; k <= Dimen; k++) cout << a[k] << " ";

    cout<<endl<<"Solving polynomial equation :"<<endl;


    for(int i=0; i<=Dimen; i++)
    {

        coefficients[Dimen-i]=a[i];
        //cin>>a[n-i];

    }
    //cout<<coefficients[0];
    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            Inverse[i][j]=(-1 )* Inverse[i][j]/coefficients[0];
        }
        //cout<<endl;
    }
    cout<<endl<<"So the Inverse of the given matrix is :"<<endl;

    for(int i=0; i<Dimen; i++)
    {
        for(int j=0; j<Dimen; j++)
        {
            cout<<Inverse[i][j]<<" ";
        }
        cout<<endl;
    }


    Get_Coefficient(coefficients,Dimen);
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


