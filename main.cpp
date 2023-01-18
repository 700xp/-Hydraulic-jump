#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>

#define eps 0.001 
#define M 1000

using namespace std;


double *J_in = new double[M];
double *J_out = new double[M];


double f(double x, double y)
{
    return 1.2*y*y*y*x*x-x;
}

double g(double x, double y)
{
    return y-3*y*y*y*y*y*x*x*x*x; 
}

double **RGK(double a, double b, double h, double X0, double Y0)
{
    int n = (b-a) / h;
    double *T = new double[n+1];
    double *X = new double[n+1];
    double *Y = new double[n+1];

    double **Res = new double*[2];
    for(int i = 0; i < 2; i++)
        Res[i] = new double[n+1];

    double K1,K2,K3,K4;
    double L1,L2,L3,L4;

    X[0] = X0;
    Y[0] = Y0;
    T[0] = a;

        for(int i = 1; i < n+1; i++)
    {
        T[i] = a + h*i;

        K1 = g(X[i-1], Y[i-1]);
        L1 = f(X[i-1], Y[i-1]);

        K2 = g(X[i-1] + L1*h/2, Y[i-1] + K1*h/2);
        L2 = f(X[i-1] + L1*h/2, Y[i-1] + K1*h/2);

        K3 = g(X[i-1] + L2*h/2, Y[i-1] + K2*h/2);
        L3 = f(X[i-1] + L2*h/2, Y[i-1] + K2*h/2);

        K4 = g(X[i-1] + L3*h, Y[i-1] + K3*h);
        L4 = f(X[i-1] + L3*h, Y[i-1] + K3*h);


        Y[i] = Y[i-1] + h*(K1 + 2*K2 + 2*K3 + K4)/6;
        X[i] = X[i-1] + h*(L1 + 2*L2 + 2*L3 + L4)/6;

        Res[0][i] = X[i];
        Res[1][i] = Y[i];
    }

    return Res;
}

double *Interpol_in(int m, int l, double R1, double R2, double *X, double *Y)
{
    double* R = new double[m];
    double* Y_Int = new double[m];

    for(int i = 0; i < m; i++)
    {
        R[i] = R1 + i*(R2-R1)/m;
    }
    Y_Int[0] = 0;

    for(int i = 1; i < m-1; i ++)
    {
        if(R[i] < X[l + 1])
            Y_Int[i] = Y[l] + (Y[l + 1] - Y[l]) * (R[i] - X[l]) / (X[l + 1] - X[l]);
        else
        {
            while (R[i] > X[l + 1])
                l ++;
            Y_Int[i] = Y[l] + (Y[l + 1] - Y[l]) * (R[i] - X[l]) / (X[l + 1] - X[l]);
        }
        if(fabs(R[i] - R2) < eps)
            break;
    }
    return Y_Int;
}

double *Interpol_out(int m, int l, double R1, double R2, double *X, double *Y)
{
    double *R = new double[m];
    double *Y_Int = new double[m];

    Y_Int[0] = 0;

    for(int i = 0; i < m; i++)
        R[i] = R2 - i*(R2-R1)/m;

    for(int i = 1; i < m; i ++)
    {
        if(R[i] > X[l + 1])
            Y_Int[i] = Y[l] + (Y[l + 1] - Y[l]) * (R[i] - X[l]) / (X[l + 1] - X[l]);
        else
        {
            while (R[i] < X[l + 1])
                l ++;
            Y_Int[i] = Y[l] + (Y[l + 1] - Y[l]) * (R[i] - X[l]) / (X[l + 1] - X[l]);
        }
        if(fabs(R[i] - R1) < eps)
            break;
    }
    //reverse array Y_int
    double temp = 0;
    for(int i = 0; i < (m/2); i++)
    {
        temp = Y_Int[i];
        Y_Int[i] = Y_Int[m-i-1];
        Y_Int[m-i-1] = temp;
    }
    return Y_Int;
}


double Jump(int m, double *R, double *Y_in, double *Y_out)
{
    double R_jump = 0;
    double *J_in = new double[m];
    double *J_out = new double[m];

    int k = 0;

    for(int i = 1; i < m-1; i++)
    {
        J_in[i] = Y_in[i] / R[i] + 0.5 / (Y_in[i] * Y_in[i] * R[i] * R[i]);
        J_out[i] = Y_out[i] / R[i] + 0.5 / (Y_out[i] * Y_out[i] * R[i] * R[i]);
    }


    for(int i = 1; i < m-1; i ++)
    {
        if((J_in[i] - J_out[i] < 0) & (k == 0))
            {
                R_jump = R[i]; 
                k++;
            }
    }
    return R_jump;
}

int main(void)
{
    //solve for t in [a,b] with step h
    double a = 0;
    double b = 10;
    double h = 0.0001;
    int n = (b-a)/h;
    int m = 1000; //size of interpol arrai

    double *T = new double[n+1];

    double *Xin = new double[n+1];
    double *Yin = new double[n+1];

    double *Xout = new double[n+1];
    double *Yout = new double[n+1];


    double **RK1 = new double*[n+1];
    for(int i = 0; i < 2; i++)
        RK1[i] = new double[n+1];

    double **RK2 = new double*[n+1];
    for(int i = 0; i < 2; i++)
        RK2[i] = new double[n+1];

    double R1, R2; //R1 - радиус для Rout, R2 - радиус для Rin
    int k1, k2, l_in1, l_out1, l_in2, l_out2;
    k1 = 0;
    k2 = 0;
    l_in1 = 0, l_in2 = 0;  
    l_out1 = 0, l_out2 = 0;

    double* R = new double[m]; //сетка на отрезке [R1,R2]

    //start condition
    T[0] = a;

    double Rf = 10;
    double Us = 5;

    double R_jump = 0;
    ofstream mout;
    mout.open("file.txt");

    for(Us = 5; Us < 10; Us+=0.1)
{
    Xout[0] = Rf;
    Yout[0] = pow(Rf,-1./3)*0.94;

    Xin[0] = 0.3;
    Yin[0] = Us;

    RK1 = RGK(a,b,h,Xout[0],Yout[0]);
    RK2 = RGK(a,b,h,Xin[0],Yin[0]);

    RK1[0][0] = Xout[0];
    RK1[1][0] = Yout[0];

    RK2[0][0] = Xin[0];
    RK2[1][0] = Yin[0];

    for(int i = 1; i < n+1; i++)
    {
        Xout[i] = RK1[0][i];
        Yout[i] = RK1[1][i];
    }

    for(int i = 1; i < n+1; i++)
    {
        Xin[i] = RK2[0][i];
        Yin[i] = RK2[1][i];
    }


    for(int i = 2; i < n+1; i++)
    {
        if(((f(Xout[1],Yout[1])*f(Xout[i],Yout[i])) < 0) & (k1 == 0))
        {
            R1 = Xout[i];
            k1 ++;
            l_out1 = i;
        }
        if((f(Xin[1],Yin[1])*f(Xin[i],Yin[i]) < 0) & (k2 == 0))
        {
            R2 = Xin[i];
            k2++;
            l_in2 = i;
        }
    }
    for(int i = 0; i < l_out1; i++)
    {
        if(fabs(Xout[i] - R2) < 0.001)
            l_out2 = i;
    }

    for (int i = 0; i < l_in2; i++)
    {
        if(fabs(Xin[i] - R1) < 0.005)
            l_in1 = i;
    }

    for(int i = 0; i < m; i ++)
    {
        R[i] = R1 + i * (R2 - R1)/m;
    }

    double* YoutINT = new double[m];
    double* YinINT = new double[m];

    YoutINT[0] = Yout[l_out2];
    YinINT[0] = Yin[l_in1];



    YoutINT = Interpol_out(m,l_out2,R1,R2,Xout,Yout);
    for(int i = 1; i < m-1; i++)

    R_jump = Jump(m,R,YinINT,YoutINT);

    mout << Us <<"\t" << R_jump << endl;

}
return 0;
}
