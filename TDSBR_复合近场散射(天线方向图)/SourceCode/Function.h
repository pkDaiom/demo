#ifndef Function_H
#define Function_H 
#include"MyClass.h"

int Sgn(double d);

void CFFT2(Complex *X, int N1, int N2, int kind);

void GSRN(double U, double G, double R, int N, double *A);

double sinc(double x);

//��˹���弰��ص���
double GaossPulse(double t, int DOrder);

//�ź���
void GateFunction(double g1, double g2, double xishu, vector<TDiracFunc<double>> &Vc);

void ProgressBar(int TotalNUM, int TotalN, int ChildNUM, int ChildN);

#endif