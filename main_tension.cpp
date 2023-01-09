#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>

#define eps 0.0005
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
        //if(f(R[0],YinINT[0])*f(R[i],YinINT[i]) < 0)
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
    //переворачиваем массив Y_int
    double temp = 0;
    for(int i = 0; i < (m/2); i++)
    {
        temp = Y_Int[i];
        Y_Int[i] = Y_Int[m-i-1];
        Y_Int[m-i-1] = temp;
    }
    return Y_Int;
}


double Jump(int m, double *R, double *Y_in, double *Y_out, double Bo)
{
    double R_jump = 0;
    double *J_in = new double[m];
    double *J_out = new double[m];

    int k = 0;

    for(int i = 1; i < m-1; i++)
    {
        J_in[i] = 1.2 * Y_in[i] / R[i] + 0.5 / (Y_in[i] * Y_in[i] * R[i] * R[i]);
        J_out[i] = 1.2 * Y_out[i] / R[i] + 0.5 / (Y_out[i] * Y_out[i] * R[i] * R[i]);
    }

    for(int i = 1; i < m-1; i++)
    {

        if (J_in[i] - J_out[i] < 1/Bo * (Y_in[i] - Y_out[i])/(R[i] * R[i] * Y_in[i] * Y_out[i]) & (k == 0))
        {
            R_jump = R[i-1] - (J_in[i-1] - J_out[i]) * (R[i-1] - R[i]) / (J_in[i-1] - J_out[i]);
            k++;
        }
    }
    return R_jump;
}

int main(void)
{
    //решение ищЄм на отрезке t \in [a,b] с шагом h
    double a = 0;
    double b = 10;
    double h = 0.0001;
    int n = (b-a)/h;
    int m = 1000; //размер массивов дл€ интерпол€ции

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

    double R1, R2; //R1 - радиус дл€ Rout, R2 - радиус дл€ Rin
    int k1, k2, l_in1, l_out1, l_in2, l_out2;
    k1 = 0;
    k2 = 0;
    l_in1 = 0, l_in2 = 0;  //индекс в массиве Xin[i], которому соответствует R2
    l_out1 = 0, l_out2 = 0; //индекс в массиве Xout[i], которому соответствует R1

    double* R = new double[m]; //сетка на отрезке [R1,R2]

    ofstream nout;
    nout.open("natyaz7.txt"); //файл дл€ вывода графика (Bo, R_jump)

    ofstream mout;
    mout.open("critBo8.txt"); //файл дл€ вывода критического значени€ Bo при Us = n

    //начальные услови€
    T[0] = a;

    double Rf = 10;
    double Us = 8;

    double R_jump = 0;

    double Bo = 10;


    for(Bo = 10; Bo > 0; Bo -= 0.02)
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

//    ofstream fout;
//    fout.open("file_in.txt");
//    ofstream gout;
//    gout.open("file_out.txt");

//    gout << Xout[0] <<"\t" << Yout[0] << endl;
//    fout << Xin[0] <<"\t" << Yin[0] << endl;

    for(int i = 1; i < n+1; i++)
    {
        Xout[i] = RK1[0][i];
        Yout[i] = RK1[1][i];
        //gout << Xout[i] <<"\t" << Yout[i] << endl;
    }

    for(int i = 1; i < n+1; i++)
    {
        Xin[i] = RK2[0][i];
        Yin[i] = RK2[1][i];
        //fout << Xin[i] <<"\t" << Yin[i] << endl;
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
    //определ€ем какой по номеру в массиве внешнего решени€ по пор€дку идет элемент = R2
    for(int i = 0; i < l_out1; i++)
    {
        if(fabs(Xout[i] - R2) < 0.001)
            l_out2 = i;
    }

    //определ€ем какой по номеру в массив внутреннего решени€ по пор€дку идет элемент = R1
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

//    ofstream jout;
//    jout.open("file_INT_in.txt");

//    jout << R[0] << "\t" << YinINT[0] << endl;
    YinINT = Interpol_in(m,l_in1,R1,R2,Xin,Yin);
//    for(int i = 1; i < m-1; i ++)
//        jout << R[i] << "\t" << YinINT[i] << endl;
//    jout.close();

//    ofstream hout;
//    hout.open("file_INT_out.txt");

    YoutINT = Interpol_out(m,l_out2,R1,R2,Xout,Yout);
//    for(int i = 1; i < m-1; i++)
//        hout << R[i] << "\t" << YoutINT[i] << endl;
//    hout.close();

    R_jump = Jump(m,R,YinINT,YoutINT, Bo);

    if (fabs(R_jump - R1) < eps)
    {
//        mout << 1./Bo << "\t" << 2 << endl;
        mout << 1./Bo << "\t" << R_jump << endl;
        mout.close();
//        mout << 1./Bo << "\t" << 0 << endl;
//        nout << 1./Bo << "\t" << R_jump << endl;
        return 1;
    }
//    ofstream kout;
//    kout.open("file.txt");  //файл со скачком


/*    for(int i = 1; i < m-1; i++)
    {
        if(R[i] <= R_jump)
            kout << R[i] << "\t" << YinINT[i] << endl;
        else
            kout << R[i] << "\t" << YoutINT[i] << endl;
    }
    kout.close();*/

//    nout << 1./Bo << "\t" << R_jump << endl;

}
return 0;
}
