#include "Variables.h"
#include"Function.h"
#include "mkl_service.h"
#include "mkl_dfti.h"
int Sgn(double d)
{
	return d < 0.0 ? -1 : 1;
}

void CFFT2(Complex *X, int N1, int N2, int kind)
{
	int index;

	/* Pointer to input/output data */
	MKL_Complex16 *x = (MKL_Complex16*)mkl_malloc(N1*N2 * sizeof(MKL_Complex16), 64);

	for (index = 0; index < N1 * N2; index++)
	{
		x[index].real = X[index].x;
		x[index].imag = X[index].y;
	}

	/* Execution status */
	MKL_LONG status = 0;

	DFTI_DESCRIPTOR_HANDLE hand = 0;

	//"Create DFTI descriptor\n"
	MKL_LONG N[2]; N[0] = N1; N[1] = N2;
	status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_COMPLEX, 2, N);

	// "Commit DFTI descriptor\n"
	status = DftiCommitDescriptor(hand);

	if (kind == 0)//"Compute forward transform
		status = DftiComputeForward(hand, x);

	if (kind == 1)//"Compute backward transform
		status = DftiComputeBackward(hand, x);

	//Free DFTI descriptor
	DftiFreeDescriptor(&hand);

	for (index = 0; index < N1 * N2; index++)
	{
		X[index].x = x[index].real;
		X[index].y = x[index].imag;
	}
	mkl_free(x);
}

void GSRN(double U, double G, double R, int N, double *A)
{
	double S = 65536.0;
	double W = 2053.0;
	double V = 13849.0;
	double T;
	double M;
	int i, j;
	for (i = 0; i < N; i++)
	{
		T = 0.0;
		for (j = 0; j < 12; j++)
		{
			R = W * R + V;
			M = floor(R / S);
			R = R - M*S;
			T = T + R / S;
		}
		A[i] = U + G*(T - 6.0);
	}

}

double sinc(double x)
{
	if (x<EPS && x>-EPS)
	{
		return 1.0;
	}
	else
	{
		return sin(x) / x;
	}
}

//高斯脉冲及相关导数
double GaossPulse(double t, int DOrder)
{
	double t1, temp, temp1, temp2;
	double w0;
	int Order = 3;//激励脉冲为高斯调制脉冲的Order阶导数
	w0 = 2.*pi*f0;
	t1 = t - t0;
	temp = exp(-4.*pi*pow(t1, 2) / pow(tao, 2));

	switch (Order + DOrder)
	{
	case 0:
		return cos(w0*t1)*temp;
	case 1:
		return (-8 * pi*t1 / pow(tao, 2)*cos(w0*t1) - w0*sin(w0*t1))*temp;
	case 2:
		return ((64 * pow(pi, 2)* pow(t1, 2) / pow(tao, 4) - 8 * pi / pow(tao, 2) - pow(w0, 2))*cos(w0*t1) + 16 * pi*w0*t1 / pow(tao, 2)*sin(w0*t1))*temp;
	case 3:
		temp1 = -512 * pow(pi, 3) * pow(t1, 3) / pow(tao, 6) + 192 * pow(pi, 2) * t1 / pow(tao, 4) + 24 * pi*pow(w0, 2) * t1 / pow(tao, 2);
		temp2 = -192 * pow(pi, 2) * w0*pow(t1, 2) / pow(tao, 4) + 24 * pi*w0 / pow(tao, 2) + pow(w0, 3);
		return (temp1*cos(w0*t1) + temp2*sin(w0*t1))*temp;
	case 4:
		temp1 = 4096 * pow(pi, 4) * pow(t1, 4) / pow(tao, 8) - 3072 * pow(pi, 3) * pow(t1, 2) / pow(tao, 6)
			- 384 * pow(pi, 2) * pow(w0, 2) * pow(t1, 2) / pow(tao, 4) + 192 * pow(pi, 2) / pow(tao, 4)
			+ 48 * pi*pow(w0, 2) / pow(tao, 2) + pow(w0, 4);
		temp2 = 2048 * pow(pi, 3) * w0*pow(t1, 3) / pow(tao, 6) - 768 * pow(pi, 2) * w0 * t1 / pow(tao, 4)
			- 32 * pi*pow(w0, 3) * t1 / pow(tao, 2);
		return (temp1*cos(w0*t1) + temp2*sin(w0*t1))*temp;
	case 5:
		temp1 = -32768 * pow(pi, 5) * pow(t1, 5) / pow(tao, 10) + 40960 * pow(pi, 4) * pow(t1, 3) / pow(tao, 8)
			+ 5120 * pow(pi, 3) * pow(w0, 2) * pow(t1, 3) / pow(tao, 6) - 7680 * pow(pi, 3) * t1 / pow(tao, 6)
			- 1920 * pow(pi, 2) * pow(w0, 2) * t1 / pow(tao, 4) - 40 * pi*pow(w0, 4) * t1 / pow(tao, 2);
		temp2 = -20480 * pow(pi, 4) * w0*pow(t1, 4) / pow(tao, 8) + 15360 * pow(pi, 3) * w0*pow(t1, 2) / pow(tao, 6)
			+ 640 * pow(pi, 2) * pow(w0, 3) * pow(t1, 2) / pow(tao, 4) - 960 * pow(pi, 2) * w0 / pow(tao, 4)
			- 80 * pi*pow(w0, 3) / pow(tao, 2) - pow(w0, 5);
		return (temp1*cos(w0*t1) + temp2*sin(w0*t1))*temp;
	default:
		cout << "该高斯脉冲的" << DOrder << "阶导数不存在！无法进行后续计算";
		return 0;
	}
}

//门函数
void GateFunction(double g1, double g2, double xishu, vector<TDiracFunc<double>> &Vc)
{
	double temp;
	if (abs(g1 - g2)<EPS)
	{
		Vc.push_back(TDiracFunc<double>(g1, 0, xishu));
	}
	else
	{
		temp = 1.0 / (g1 - g2);
		Vc.push_back(TDiracFunc<double>(g1, -1, xishu*temp));
		Vc.push_back(TDiracFunc<double>(g2, -1, -xishu*temp));
	}
}

void ProgressBar(int TotalNUM, int TotalN, int ChildNUM, int ChildN)
{
	int N = 50;
	//int i, j;	
	if (ChildN >= ChildNUM)
		ChildN = ChildNUM;
	if (TotalN >= TotalNUM)
	{
		TotalN = TotalNUM;
		ChildN = ChildNUM;
	}
	//printf("\r%.1lf%%", a * 100.0 / NUM, a * 10.0 / NUM);
	cout << "\r" << fixed << setprecision(1) << "总进程：" << TotalN * 100.0 / TotalNUM << "%" << "  子进程：" << ChildN * 100.0 / ChildNUM << "%";

	/*cout << "[";
	for (i = 0; i<a*N / NUM; i++)
	{
	cout << "▇";
	}
	for (j = a*N / NUM; j<N; j++)
	cout << " ";
	cout << "]";*/
}
