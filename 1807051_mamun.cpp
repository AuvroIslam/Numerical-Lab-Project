#include<bits/stdc++.h>
using namespace std;


bool LU_Factorization(vector<vector<double>>&coefficient,vector<vector<double>>&L,vector<vector<double>>&U)
{
    int n=coefficient.size();
    L=vector<vector<double>>(n,vector<double>(n,0.0));
    U=vector<vector<double>>(n,vector<double>(n,0.0));
    for(int i=0;i<n;i++)
    {
        for(int k=i;k<n;k++)
        {
            double sum=0;
            for(int j=0;j<i;j++)
            {
                sum+=(L[i][j]*U[j][k]);
            }
            U[i][k]=coefficient[i][k]-sum;
        }
        for(int k=i;k<n;k++)
        {
            if(i==k)
                L[i][i]=1;
            else
            {
                double sum=0.0;
                for(int j=0;j<i;j++)
                    sum+=(L[k][j]*U[j][i]);
                if(U[i][i]==0)
                {
                    cout<<"LU Factorization is not possible"<<endl;
                    return false;
                }
                L[k][i]=(coefficient[k][i]-sum)/U[i][i];
            }
        }
    }
    return true;
}
double fun(double x,double y)
{
    return (2*x+3*y);
    //return y;
}

double Runge_Kutta_Method(double y0, double x0, double x, double n)
{
    double xn,yn,k,k1,k2,k3,k4;
    double h=(x-x0)/n;
    for(int i=0;i<n;i++)
    {
        k1=h*(fun(x0,y0));
        k2=h*(fun((x0+h/2),(y0+k1/2)));
        k3=h*(fun((x0+h/2),(y0+k2/2)));
        k4=h*(fun((x0+h),(y0+k3)));
        k=(k1+2*k2+2*k3+k4)/6;
        yn=y0+k;
        x0+=h;
        y0=yn;
    }
    return yn;
}
