
#include<bits/stdc++.h>
using namespace std;
double  left_singular_vector(int a[2][2],int b[2][2],double Du[2][2],double u[2][2])
{

    double mat2[2][2],sum=0;
    for(int i=0;i<2;i++)
   {
       for(int j=0;j<2;j++)
       {
           for(int k=0;k<2;k++)
           {
               mat2[i][j]=a[i][k]*b[k][j];
               sum=sum+mat2[i][j];

           }
           mat2[i][j]=sum;
           sum=0;
       }

   }

      int ev1,ev2;
    if(mat2[0][0]+mat2[1][1]>0)
    {
      ev1=((mat2[0][0]+mat2[1][1])+sqrt((mat2[0][0]+mat2[1][1])*(mat2[0][0]+mat2[1][1]) - 4*(mat2[0][0]*mat2[1][1] - mat2[0][1]*mat2[1][0])))/2;
      //cout<<"Ev1 ="<<ev1<<endl;
      ev2=((mat2[0][0]+mat2[1][1])-sqrt((mat2[0][0]+mat2[1][1])*(mat2[0][0]+mat2[1][1]) - 4*(mat2[0][0]*mat2[1][1] - mat2[0][1]*mat2[1][0])))/2;
      //cout<<"Ev2 ="<<ev2<<endl;
    }
    else
    {
       ev1=(-(mat2[0][0]+mat2[1][1])+sqrt((mat2[0][0]+mat2[1][1])*(mat2[0][0]+mat2[1][1]) - 4*(mat2[0][0]*mat2[1][1] - mat2[0][1]*mat2[1][0])))/2;
       //cout<<"Ev1 ="<<ev1<<endl;
       ev2=(-(mat2[0][0]+mat2[1][1])-sqrt((mat2[0][0]+mat2[1][1])*(mat2[0][0]+mat2[1][1]) - 4*(mat2[0][0]*mat2[1][1] - mat2[0][1]*mat2[1][0])))/2;
       //cout<<"Ev2 ="<<ev2<<endl;
    }
     //cout<<"Vector :"<<endl;
     double x1=1;
     double y1=-1*(mat2[0][1])/(mat2[0][0]-ev1);
     double x2=1;
     double y2=-1*(mat2[0][1])/(mat2[0][0]-ev2);

     //singular value calculate;
      //cout<<y2<<" "<<x2<<endl;
      double eigen_vector_mat[2][2]={{y1,x1},{y2,x2}};
      double eigen_vector_mat1[1][2];//={y1,x1};
      eigen_vector_mat1[0][0]=y1;
      eigen_vector_mat1[0][1]=x1;
      double eigen_vector_mat2[1][2];//={y2,x2};
      eigen_vector_mat2[0][0]=y2;
      eigen_vector_mat2[0][1]=x2;

     double s1=sqrt(fabs(ev1));
     double s2=sqrt(fabs(ev2));
     //double Du[2][2]={{s1,0},{0,s2}};
     Du[0][0]=s1;
     Du[0][1]=0;
     Du[1][0]=0;
     Du[1][1]=s2;

    double eigen_vector_mat1_trans[2][1];
    eigen_vector_mat1_trans[0][0]=eigen_vector_mat1[0][0];
    eigen_vector_mat1_trans[1][0]=eigen_vector_mat1[0][1];
    //cout<<eigen_vector_mat1_trans[0][0]<<" "<<eigen_vector_mat1_trans[1][0];
    double denominator1=sqrt(eigen_vector_mat1[0][0]*eigen_vector_mat1_trans[0][0] + eigen_vector_mat1_trans[1][0]*eigen_vector_mat1[0][1]);
    //cout<<denominator1;
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<1;j++)
        {
            eigen_vector_mat1_trans[i][j]=eigen_vector_mat1_trans[i][j]/denominator1;
        }

    }

    double eigen_vector_mat2_trans[2][1];
    eigen_vector_mat2_trans[0][0]=eigen_vector_mat2[0][0];
    eigen_vector_mat2_trans[1][0]=eigen_vector_mat2[0][1];

    double denominator2=sqrt(eigen_vector_mat2[0][0]*eigen_vector_mat2_trans[0][0] +eigen_vector_mat2_trans[1][0]*eigen_vector_mat2[0][1]);
    //cout<<denominator1;
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<1;j++)
        {
            eigen_vector_mat2_trans[i][j]=eigen_vector_mat2_trans[i][j]/denominator2;
        }

    }
    //cout<<"Orthonormal Matrix U:"<<endl;

    u[0][0]=eigen_vector_mat1_trans[0][0];
    u[0][1]=eigen_vector_mat2_trans[0][0];
    u[1][0]=eigen_vector_mat1_trans[1][0];
    u[1][1]=eigen_vector_mat2_trans[1][0];


}

double  right_singular_vector(int a[2][2],int b[2][2],double Dv[2][2],double v[2][2])
{

    double mat2[2][2],sum=0;
    for(int i=0;i<2;i++)
   {
       for(int j=0;j<2;j++)
       {
           for(int k=0;k<2;k++)
           {
               mat2[i][j]=a[i][k]*b[k][j];
               sum=sum+mat2[i][j];

           }
           mat2[i][j]=sum;
           sum=0;
       }

   }

      int ev1,ev2;
    if(mat2[0][0]+mat2[1][1]>0)
    {
      ev1=((mat2[0][0]+mat2[1][1])+sqrt((mat2[0][0]+mat2[1][1])*(mat2[0][0]+mat2[1][1]) - 4*(mat2[0][0]*mat2[1][1] - mat2[0][1]*mat2[1][0])))/2;
      //cout<<"Ev1 ="<<ev1<<endl;
      ev2=((mat2[0][0]+mat2[1][1])-sqrt((mat2[0][0]+mat2[1][1])*(mat2[0][0]+mat2[1][1]) - 4*(mat2[0][0]*mat2[1][1] - mat2[0][1]*mat2[1][0])))/2;
      //cout<<"Ev2 ="<<ev2<<endl;
    }
    else
    {
       ev1=(-(mat2[0][0]+mat2[1][1])+sqrt((mat2[0][0]+mat2[1][1])*(mat2[0][0]+mat2[1][1]) - 4*(mat2[0][0]*mat2[1][1] - mat2[0][1]*mat2[1][0])))/2;
       //cout<<"Ev1 ="<<ev1<<endl;
       ev2=(-(mat2[0][0]+mat2[1][1])-sqrt((mat2[0][0]+mat2[1][1])*(mat2[0][0]+mat2[1][1]) - 4*(mat2[0][0]*mat2[1][1] - mat2[0][1]*mat2[1][0])))/2;
       //cout<<"Ev2 ="<<ev2<<endl;
    }
     //cout<<"Vector :"<<endl;
     double x1=1;
     double y1=-1*(mat2[0][1])/(mat2[0][0]-ev1);
     double x2=1;
     double y2=-1*(mat2[0][1])/(mat2[0][0]-ev2);

     //singular value calculate;
      //cout<<y1<<" "<<x1<<endl;
      double eigen_vector_mat[2][2]={{y1,x1},{y2,x2}};
      double eigen_vector_mat1[1][2];//={y1,x1};
      eigen_vector_mat1[0][0]=y1;
      eigen_vector_mat1[0][1]=x1;
      double eigen_vector_mat2[1][2];//={y2,x2};
      eigen_vector_mat2[0][0]=y2;
      eigen_vector_mat2[0][1]=x2;

     double s1=sqrt(fabs(ev1));
     double s2=sqrt(fabs(ev2));
     //double Dv[2][2]={{s1,0},{0,s2}};
     Dv[0][0]=s1;
     Dv[0][1]=0;
     Dv[1][0]=0;
     Dv[1][1]=s2;

    double eigen_vector_mat1_trans[2][1];
    eigen_vector_mat1_trans[0][0]=eigen_vector_mat1[0][0];
    eigen_vector_mat1_trans[1][0]=eigen_vector_mat1[0][1];
    //cout<<eigen_vector_mat1_trans[0][0]<<" "<<eigen_vector_mat1_trans[1][0];
    double denominator1=sqrt(eigen_vector_mat1[0][0]*eigen_vector_mat1_trans[0][0] +eigen_vector_mat1_trans[1][0]*eigen_vector_mat1[0][1]);
    //cout<<denominator1;
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<1;j++)
        {
            eigen_vector_mat1_trans[i][j]=eigen_vector_mat1_trans[i][j]/denominator1;
        }

    }

    double eigen_vector_mat2_trans[2][1];
    eigen_vector_mat2_trans[0][0]=eigen_vector_mat2[0][0];
    eigen_vector_mat2_trans[1][0]=eigen_vector_mat2[0][1];
    double denominator2=sqrt(eigen_vector_mat2[0][0]*eigen_vector_mat2_trans[0][0] +eigen_vector_mat2_trans[1][0]*eigen_vector_mat2[0][1]);
    //cout<<denominator1;
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<1;j++)
        {
            eigen_vector_mat2_trans[i][j]=eigen_vector_mat2_trans[i][j]/denominator2;
        }

    }

    v[0][0]=eigen_vector_mat1_trans[0][0];
    v[0][1]=eigen_vector_mat2_trans[0][0];
    v[1][0]=eigen_vector_mat1_trans[1][0];
    v[1][1]=eigen_vector_mat2_trans[1][0];


}
int main()
{
    int a[2][2],b[2][2],multi[2][2];
    double Du[2][2],Dv[2][2],u[2][2],v[2][2],vt[2][2],D[2][2],Ut[2][2],Di[2][2],V_DI[2][2],Inverse[2][2];
    cout<<"Enter the matrix:"<<endl;
    FILE *file;
    file=fopen("spl_task1.txt","r");
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            fscanf(file,"%d",&a[i][j]);
        }
    }
    /*for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            cin>>a[i][j];
        }

    }*/
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            b[i][j]=a[j][i];
        }

    }

    //cout<<"After multiply: KU"<<endl;
    left_singular_vector(a,b,Du,u);
    //cout<<"After 2nd multiply"<<endl;
    right_singular_vector(b,a,Dv,v);

    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            vt[i][j]=v[j][i];
        }

    }


    D[0][0]=Du[0][0]*Dv[0][0] + Du[0][1]*Dv[1][0];
    D[0][1]=Du[0][0]*Dv[0][1] + Du[0][1]*Dv[1][1];
    D[1][0]=Du[1][0]*Dv[0][0] + Du[1][1]*Dv[1][0];
    D[1][1]=Du[1][0]*Dv[0][1] + Du[1][1]*Dv[1][1];

    cout<<"So the SVD of given matrix is the form of UV^tD is:"<<endl<<endl;
    cout<<"U :";
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            cout<<"     "<<u[i][j]<<" ";
        }cout<<"  "<<endl;
    }
    cout<<"V^t :";
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            cout<<"     "<<vt[i][j]<<" ";
        }cout<<"  "<<endl;
    }
    cout<<"D :";
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            cout<<"     "<<D[i][j]<<" ";
        }cout<<"  "<<endl;
    }
    //pseudo inverse;
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            if(D[i][j]!=0)
             Di[i][j]=1/D[j][i];
             else
               Di[i][j]=0;
        }

    }
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            Ut[i][j]=u[j][i];
        }

    }
    V_DI[0][0]=v[0][0]*Di[0][0] + v[0][1]*Di[1][0];
    V_DI[0][1]=v[0][0]*Di[0][1] + v[0][1]*Di[1][1];
    V_DI[1][0]=v[1][0]*Di[0][0] + v[1][1]*Di[1][0];
    V_DI[1][1]=v[1][0]*Di[0][1] + v[1][1]*Di[1][1];

    Inverse[0][0]=V_DI[0][0]*Ut[0][0] + V_DI[0][1]*Ut[1][0];
    Inverse[0][1]=V_DI[0][0]*Ut[0][1] + V_DI[0][1]*Ut[1][1];
    Inverse[1][0]=V_DI[1][0]*Ut[0][0] + V_DI[1][1]*Ut[1][0];
    Inverse[1][1]=V_DI[1][0]*Ut[0][1] + V_DI[1][1]*Ut[1][1];
    cout<<endl<<"Pseudo Inverse for given matrix is the form of V.D^-1.U^t is :"<<endl;
    for(int i=0;i<2;i++)
      {
        for(int j=0;j<2;j++)
        {
            cout<<Inverse[i][j]<<" ";
        }
        cout<<endl;

      }

}
