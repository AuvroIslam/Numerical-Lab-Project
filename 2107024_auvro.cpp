#include<bits/stdc++.h>
const double max_tolerance=0.00001;
using namespace std;
double eqn(vector<double> &coefficients,int degree, double x)
{
    double result=0.0;
    for(int i=0;i<=degree;i++)
    {
        result+=coefficients[i]*(pow(x,degree-i));
    }
    return result;
}

double Bi_Section_Method(int degree,vector<double>&coefficients,double a,double b)
{
    if(eqn(coefficients,degree,a)*eqn(coefficients,degree,b)>=0)
    {
        cout<<"Invalid interval for Bisection Method"<<endl;
        return -1;
    }
    double c=a;
    while((b-a)>=max_tolerance)
    {
        c=(a+b)/2;
        if(eqn(coefficients,degree,c)==0.0)
        {
            return c;
        }
        else if(eqn(coefficients,degree,c)*eqn(coefficients,degree,a)<0)
            b=c;
        else if(eqn(coefficients,degree,c)*eqn(coefficients,degree,b)<0)
            a=c;
    }
    return c;
}

double Newton_Raphson_Method(int degree,vector<double>& coefficients,vector<double>&deriCoefficient,double x0)
{
    double x=x0;
    while(fabs(eqn(coefficients,degree,x))>max_tolerance)
    {
        double fx= eqn(coefficients,degree,x);
        double dfx= eqn(deriCoefficient,degree-1,x);
        if(fabs(dfx)<max_tolerance)
        {
            cout<<"Derivative near zero, may not converge."<<endl;
            return x;
        }
        x=x-(fx/dfx);
    }
    return x;
}

double Secant_Method(int degree,vector<double>& coefficients,double x1, double x2)
{
    double x3;
    while(true)
    {
        double fx1=eqn(coefficients,degree,x1);
        double fx2=eqn(coefficients,degree,x2);
        if(fabs(fx1-fx2)<max_tolerance)
        {
            cout<<"Secant method failed due to small denominator"<<endl;
            return x2;
        }
        x3=x2-(fx2*(x2-x1)/(fx2-fx1));
        if(eqn(coefficients,degree,x3)==0.0)
        {
            return x3;
        }
        if(fabs(x3-x2)<=max_tolerance)
        {
            return x3;
        }
        x1=x2;
        x2=x3;
    }
    return x3;
}

double False_Position_Method(int degree, vector<double>&coefficients,double a,double b)
{
    if(eqn(coefficients,degree,a)*eqn(coefficients,degree,b)>=0)
    {
        cout<<"Invalid interval for False Position method"<<endl;
        return -1;
    }
    double old_x=a;
    while(true)
    {
        double x=(a*eqn(coefficients,degree,b)-b*eqn(coefficients,degree,a))/(eqn(coefficients,degree,b)-eqn(coefficients,degree,a));
        if(eqn(coefficients,degree,x)==0.0)
        {
            return x;
        }
        else if(eqn(coefficients,degree,x)*eqn(coefficients,degree,a)<0)
            b=x;
        else if(eqn(coefficients,degree,x)*eqn(coefficients,degree,b)<0)
            a=x;
        if(fabs(x-old_x)<=max_tolerance)
        {
            return x;
        }
        old_x=x;
    }
    return old_x;
}


