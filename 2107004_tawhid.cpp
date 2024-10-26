
#include<bits/stdc++.h>
using namespace std;
const double max_tolerance=0.00001;
const int max_loop=100;
void swap_rows(vector<vector<double>>&matrix, int row1,int row2,int sz)
{
    for(int i=0;i<=sz;i++)
    {
        swap(matrix[row1][i],matrix[row2][i]);
    }
}
void Jacobi_Iterative_Method(vector<vector<double>> &coefficient,vector<double> &rhs,vector<double>&initial)
{
    int n=coefficient.size();
    vector<double> x(n,0.0);
    for(int i=0;i<max_loop;i++)
    {
        for(int j=0;j<n;j++)
        {
            double approx=0.0;
            for(int k=0;k<n;k++)
            {
                if(j!=k)
                {
                    approx=approx+coefficient[j][k]*initial[k];
                }
            }
            x[j]=(rhs[j]-approx)/coefficient[j][j];
        }
        double max_error=0.0;
        for(int k=0;k<n;k++)
        {
            max_error=max(max_error,fabs(x[k]-initial[k]));

        }
        if(max_error<max_tolerance)
        {
            break;
        }
        initial=x;
    }
    for(int i=0;i<n;i++)
    {
        cout<<"x"<<i+1<<"= "<<initial[i]<<endl;
    }
}

void Gauss_Siedel_Method(vector<vector<double>> &coefficient,vector<double> &rhs,vector<double>&initial)
{
    int n=coefficient.size();
    for(int i=0;i<max_loop;i++)
    {
        double max_error=0.0;
        for(int j=0;j<n;j++)
        {
            double approx=0.0;
            for(int k=0;k<n;k++)
            {
                if(j!=k)
                {
                    approx=approx+coefficient[j][k]*initial[k];
                }
            }
            double new_value=(rhs[j]-approx)/coefficient[j][j];
            max_error=max(max_error,fabs(new_value-initial[j]));
            initial[j]=new_value;

        }
        if(max_error<max_tolerance)
        {
            break;
        }
    }
    for(int i=0;i<n;i++)
    {
        cout<<"x"<<i+1<<"= "<<initial[i]<<endl;
    }
}


vector<vector<double>> Gauss_Elimination_Method(int sz,vector<vector<double>>&A1, vector<double> &rhs)
{
    for(int i=0;i<sz;i++)
    {
        A1[i].push_back(rhs[i]);
    }
    for(int i=0;i<sz;i++)
    {
        if(A1[i][i]==0)
        {
            bool flag=false;
            for(int k=i+1;k<sz;k++)
            {
                if(A1[k][i]!=0)
                {
                    swap_rows(A1,i,k,sz);
                    flag=true;
                    break;
                }
            }
            if(!flag)
            {
                cerr<<"No suitable rows found for swapping"<<endl;
                return A1;
            }
        }
        for(int j=sz-1;j>i;j--)
        {
            double vag_korbo=A1[j][i]/A1[i][i];
            if(vag_korbo!=0.0)
            {
                for(int k=0;k<sz+1;k++)
                {
                    A1[j][k]=A1[j][k]-(vag_korbo*A1[i][k]);
                }
            }
        }

    }
    return A1;
}

vector<vector<double>> Gauss_Jordan_Method(int sz, vector<vector<double>> &A2)
{
    for(int i=sz-1;i>=0;i--)
    {
        if(A2[i][i]==0)
        {
            bool flag=false;
            for(int k=i+1;k<sz;k++)
            {
                if(A2[k][i]!=0)
                {
                    swap_rows(A2,i,k,sz);
                    flag=true;
                    break;
                }
            }
            if(!flag)
            {
                cerr<<"No suitable rows found for swapping"<<endl;
                return A2;
            }
        }
        for(int j=0;j<i;j++)
        {
            double vag_korbo=A2[j][i]/A2[i][i];
            if(vag_korbo!=0.0)
            {
                for(int k=0;k<sz+1;k++)
                {
                    A2[j][k]=A2[j][k]-(vag_korbo*A2[i][k]);
                }
            }
        }
    }
    return A2;
}

vector<vector<double>> row_echelon(int sz, vector<vector<double>> &A3)
{
    for(int i=0,j=0,k=sz;i<sz && j<sz;i++,j++)
    {
        double solution=A3[i][k]/A3[i][j];
        A3[i][k]=solution;
        A3[i][j]=1;
    }
    return A3;
}
vector<vector<double>> Matrix_Inversion(vector<vector<double>>&A)
{
    int n=A.size();
    vector<vector<double>> inv(n,vector<double>(n,0.0));
    for(int i=0;i<n;i++)
        inv[i][i]=1.0;
    for(int i=0;i<n;i++)
    {
        double dig=A[i][i];
        for(int j=0;j<n;j++)
        {
            A[i][j]/=dig;
            inv[i][j]/=dig;
        }
        for(int k=0;k<n;k++)
        {
            if(k!=i)
            {
                double factor=A[k][i];
                for(int j=0;j<n;j++)
                {
                    A[k][j]-=factor*A[i][j];
                    inv[k][j]-=factor*inv[i][j];
                }
            }
        }
    }
    return inv;
}

