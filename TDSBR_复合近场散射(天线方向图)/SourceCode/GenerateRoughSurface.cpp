#include "Variables.h"
#include"Function.h"
#include"GenerateRoughSurface.h"
//高斯谱
void Gauss_Spm(double *S, RoughSurface &roughsurface)
{
	int i, j;
	int Nx = roughsurface.LNx;
	int Ny = roughsurface.LNy;
	double delta = roughsurface.delta;
	double rlx = roughsurface.rlx;
	double rly = roughsurface.rly;
	double deltakx = 2.0*pi / roughsurface.Lx;
	double deltaky = 2.0*pi / roughsurface.Ly;
	double temp1 = 2.0*pi*sqrt(roughsurface.Lx*roughsurface.Ly);
	for (i = 0; i <= Nx / 2; i++)
	{
		for (j = -Ny / 2 + 1; j <= Ny / 2; j++)
		{
			S[i*Ny + (j + Ny / 2 - 1)] = temp1*sqrt(delta*delta * rlx*rly*exp(-(pow(deltakx*i*rlx, 2.0) + pow(deltaky*j*rly, 2.0)) / 4.) / (4.0 * pi));
		}
	}
}

void PM_Spm(double *S, RoughSurface &sea)
{
	int i, j;
	double alpha = 8.1e-3;
	double beta = 0.74;
	double g = 9.81;
	int Nx = sea.LNx;
	int Ny = sea.LNy;
	double deltakx = 2.0*pi / sea.Lx;
	double deltaky = 2.0*pi / sea.Ly;
	double k;
	double temp1 = 2.0*pi*sqrt(sea.Lx*sea.Ly);
	double Ss;
	double phi;
	for (i = 0; i <= Nx / 2; i++)
	{
		for (j = -Ny / 2 + 1; j <= Ny / 2; j++)
		{
			k = sqrt(pow(deltakx*i, 2) + pow(deltaky*j, 2));
			if (k == 0.0) k = EPS;
			Ss = alpha / (2 * pow(k, 3.0)) * exp(-beta*g*g / (k*k*pow(sea.WindSpeed, 4)));
			phi = atan((deltaky*j + 1.0e-8) / (deltakx*i + 1.0e-8));
			Ss = Ss * pow(cos((phi - sea.WindTheta*pi / 180.0) / 2.0), 2.0) / pi / k;
			S[i*Ny + (j + Ny / 2 - 1)] = temp1*sqrt(Ss);
		}
	}
}

double C(double k)
{
	return sqrt(9.8*(1 + (k*k) / (370 * 370)) / k);
}

void ELH_Spm(double *S, RoughSurface &sea)
{
	int i, j;
	double uf = 0.3286, g = 9.8, X0 = 2.2e4, km = 370.0, X = 3.0e4;
	double k, u, temp_X, omega_c, kp, omega, alfa_p, alfa_m, Lpm, gamma, sigma_j, temp_B, Jp, Fp, BL, Bs;
	int Nx = sea.LNx;
	int Ny = sea.LNy;
	double deltakx = 2.0*pi / sea.Lx;
	double deltaky = 2.0*pi / sea.Ly;
	double temp1 = 2.0*pi*sqrt(sea.Lx*sea.Ly);
	double Ss;
	double phi;
	u = sea.WindSpeed;
	temp_X = g*X / pow(u, 2);
	omega_c = 0.84*pow(tanh(pow(temp_X / X0, 0.4)), -0.75);
	kp = pow(omega_c, 2) * g / (u*u);
	omega = u / C(kp);
	alfa_p = 0.006*sqrt(omega);

	if (omega_c <= 1)
		gamma = 1.7;
	else if (omega_c >= 5)
		gamma = 2.7*pow(omega_c, 0.57);
	else
		gamma = 1.7 + 6 * log10(omega_c);
	if (omega_c < 5)
		sigma_j = 0.08*(1 + 4 * pow(omega_c, -3));
	else
		sigma_j = 0.16;
	for (i = 0; i <= Nx / 2; i++)
	{
		for (j = -Ny / 2 + 1; j <= Ny / 2; j++)
		{
			k = sqrt(pow(deltakx*i, 2) + pow(deltaky*j, 2));
			if (k == 0) k = EPS;
			Lpm = exp(-1.25*pow(k / kp, -2));
			temp_B = exp(-pow(pow(k / kp, 0.5) - 1, 2) / (2 * sigma_j*sigma_j));
			Jp = pow(gamma, temp_B);
			Fp = Lpm*Jp*exp(-omega / sqrt(10.0)*(pow(k / kp, 0.5) - 1));
			BL = 0.5*alfa_p*C(kp) / C(k)*Fp;
			if (uf <= 0.23)
				alfa_m = 0.01*(1 + log(uf / 0.23));
			else
				alfa_m = 0.01*(1 + 3 * log(uf / 0.23));
			Bs = 0.5*alfa_m*C(km) / C(k)*exp(-0.25*pow(k / km - 1.0, 2))*Lpm*Jp;
			Ss = (Bs + BL) / pow(k, 3.0);
			phi = atan((deltaky*j + 1.0e-8) / (deltakx*i + 1.0e-8));
			Ss = Ss * pow(cos((phi - sea.WindTheta*pi / 180.0) / 2.0), 2.0) / pi / k;
			S[i*Ny + (j + Ny / 2 - 1)] = temp1*sqrt(Ss);

		}
	}
}

// F为频率，T为温度，S_sw是盐度，世界大洋的平均盐度为35‰，这里S_sw取34.7。
void SeaDielec(double F, RoughSurface &sea)
{
	double T = sea.Temp;
	double S_sw = sea.Sault;
	double epsilon_sw_inf = 4.9;
	double epsilon_0 = 8.854e-12;

	// 首先计算离子电导率sigma_i_T_S

	double delta = 25.0 - T;
	double phi = delta *(2.033e-2 + 1.266e-4 * delta + 2.464e-6 * delta * delta - S_sw * (1.849e-5 - 2.551e-7 * delta + 2.551e-8 * delta * delta));
	double sigma_i_25_S = S_sw * (0.18252 - 1.4619e-3 * S_sw + 2.093e-5 * S_sw * S_sw - 1.282e-7 * S_sw * S_sw * S_sw);
	double sigma_i_T_S = sigma_i_25_S * exp(-phi);

	// 计算海水张弛时间tao_sw_T_S

	double tao_sw_T_0 = 1.1109e-10 - 3.824e-12 * T + 6.938e-14 * T * T - 5.096e-16 * T * T * T;
	double b_T_S = 1.0 + 2.282e-5 * T * S_sw - 7.638e-4 * S_sw - 7.760e-6 * S_sw * S_sw + 1.105e-8 * S_sw * S_sw * S_sw;
	double tao_sw_T_S = tao_sw_T_0 * b_T_S;

	// 计算静态介电常数epsilon_sw0_T_S
	double epsilon_sw0_T_0 = 87.134 - 1.949e-1 * T - 1.276e-2 * T * T + 2.491e-4 * T * T * T;
	double a_T_S = 1.0 + 1.613e-5 * T * S_sw - 3.656e-3 * S_sw + 3.210e-5 * S_sw * S_sw - 4.232e-7 * S_sw * S_sw * S_sw;
	double epsilon_sw0_T_S = epsilon_sw0_T_0 * a_T_S;

	double temp1 = 1.0 + (F * tao_sw_T_S) * (F * tao_sw_T_S);
	double temp2 = epsilon_sw0_T_S - epsilon_sw_inf;
	double temp = temp2 / temp1;

	sea.er.x = epsilon_sw_inf + temp;
	sea.er.y = tao_sw_T_S * temp * F + sigma_i_T_S / (2.0 * pi * epsilon_0 * F);
}

void GenerateRoughSurface(Target &Tar, RoughSurface &roughsurface)
{
	int i, j, m, n = 0;
	int Nx = roughsurface.LNx;
	int Ny = roughsurface.LNy;
	double Lx = roughsurface.Lx;
	double Ly = roughsurface.Ly;
	double High = roughsurface.High;

	double *A = new double[Nx*Ny];
	double *S = new double[(Nx / 2 + 1) * Ny];
	Complex *FF = new Complex[Nx*Ny];
	double *Z = new double[Nx*Ny];
	double x0 = -Lx / 2;
	double y0 = -Ly / 2;
	double dx = Lx / (Nx - 1);
	double	dy = Ly / (Ny - 1);

	Point3f P1, P2, P3, P4;
	Point3f tempv;

	GSRN(0.0, 1.0, roughsurface.Seed, Nx*Ny, A);//生成随机数
	switch (roughsurface.Spectrum)
	{
	case Gauss:
		Gauss_Spm(S, roughsurface);
		break;

	case PM:
		PM_Spm(S, roughsurface);
		break;

	case ELH:
		ELH_Spm(S, roughsurface);
		break;
	}

	FF[0 * Ny + 0] = S[0 * Ny + (Ny / 2 - 1)] * A[0];
	FF[0 * Ny + Ny / 2] = S[0 * Ny + (Ny - 1)] * A[1];
	FF[(Nx / 2) * Ny + 0] = S[Nx / 2 * Ny + (Ny / 2 - 1)] * A[2];
	FF[(Nx / 2)* Ny + Ny / 2] = S[Nx / 2 * Ny + (Ny - 1)] * A[3];

	for (m = 0; m <= Nx / 2; m = m + Nx / 2)
	{
		for (n = 1; n <= Ny / 2 - 1; n++)
		{
			FF[(m* Ny) + n] = S[m *Ny + (n + Ny / 2 - 1)] * (A[m*(Ny - 2) + n * 2 + 2] + J*A[m*(Ny - 2) + n * 2 + 3]);
			FF[(m* Ny) + (-n + Ny)] = FF[(m* Ny) + n].x - J * FF[(m* Ny) + n].y;
		}
	}

	for (m = 1; m <= Nx / 2 - 1; m++)
	{
		for (n = 1; n <= Ny / 2 - 1; n++)
		{
			FF[(m* Ny) + n] = S[m *Ny + (n + Ny / 2 - 1)] * (A[m*(Ny - 2) + n * 2 + 2] + J*A[m*(Ny - 2) + n * 2 + 3]) / sqrt(2.0);
			FF[(-m + Nx)*Ny + (-n + Ny)] = FF[(m* Ny) + n].x - J * FF[(m* Ny) + n].y;
		}
	}

	for (m = 1; m <= Nx / 2 - 1; m++)
	{
		for (n = 1; n <= Ny / 2 - 1; n++)
		{
			FF[m * Ny + (-n + Ny)] = S[m *Ny + (n + Ny / 2 - 1)] * (A[(m + Nx / 2)*(Ny - 2) + n * 2 + 2] + J*A[(m + Nx / 2)*(Ny - 2) + n * 2 + 3]) / sqrt(2.0);
			FF[(-m + Nx) * Ny + n] = FF[m * Ny + (-n + Ny)].x - J * FF[m * Ny + (-n + Ny)].y;
		}
	}

	for (n = 0; n <= Ny / 2; n = n + Ny / 2)
	{
		for (m = 1; m <= Nx / 2 - 1; m++)
		{
			FF[m * Ny + n] = S[m *Ny + (n + Ny / 2 - 1)] * (A[Nx*(Ny - 2) + (2 * n / Ny)*(Nx - 2) + m * 2 + 2] + J*A[Nx*(Ny - 2) + (2 * n / Ny)*(Nx - 2) + m * 2 + 3]);
			FF[(-m + Nx) * Ny + n] = FF[m * Ny + n].x - J * FF[m * Ny + n].y;
		}
	}
	delete[] A; A = nullptr;
	delete[] S; S = nullptr;

	CFFT2(FF, Nx, Ny, 1);

	for (m = 0; m <= Nx / 2; m++)
	{
		for (n = 0; n <= Ny / 2; n++)
		{
			Z[(m + Nx / 2 - 1)*Ny + (n + Ny / 2 - 1)] = FF[m * Ny + n].x / Lx / Ly + High;
		}
	}

	for (m = 0; m <= Nx / 2; m++)
	{
		for (n = Ny / 2 + 1; n <= Ny - 1; n++)
		{
			Z[(m + Nx / 2 - 1)*Ny + (n - Ny / 2 - 1)] = FF[m * Ny + n].x / Lx / Ly + High;
		}
	}
	for (m = Nx / 2 + 1; m <= Nx - 1; m++)
	{
		for (n = 0; n <= Ny / 2; n++)
		{
			Z[(m - Nx / 2 - 1)*Ny + (n + Ny / 2 - 1)] = FF[m * Ny + n].x / Lx / Ly + High;
		}
	}
	for (m = Nx / 2 + 1; m <= Nx - 1; m++)
	{
		for (n = Ny / 2 + 1; n <= Ny - 1; n++)
		{
			Z[(m - Nx / 2 - 1)*Ny + (n - Ny / 2 - 1)] = FF[m * Ny + n].x / Lx / Ly + High;
		}
	}

	delete[] FF; FF = nullptr;

	////////////////////////////////// 输出粗糙面/////////////////////////////////////////////
	fstream f1;
	f1.open(roughsurface.name, ios::out);
	for (m = 0; m < Nx; m++)
	{
		for (n = 0; n < Ny; n++)
		{
			f1 << x0 + dx*m << ' ' << y0 + dy*n << ' ' << Z[m*Ny + n] << endl;
		}
	}
	f1.close();
	/////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////粗糙面转换成目标类型////////////////////////////////////////////////////////
	if (roughsurface.Spectrum == Gauss)
	{
		Tar.er = roughsurface.er;
		Tar.mr = roughsurface.mr;
	}
	else
	{
		SeaDielec(f0, roughsurface);
		Tar.er = roughsurface.er;
		Tar.mr = roughsurface.mr;
	}

	Tar.TNum = (Nx - 1) * (Ny - 1) * 2;
	Tar.Create();
	m = 0;
	n = 0;
	for (i = 0; i<Nx - 1; ++i)
	{
		for (j = 0; j<Ny - 1; ++j)
		{
			P1.x = x0 + i * dx;		P1.y = y0 + j * dy;		P1.z = Z[i*Ny + j];
			P2.x = x0 + (i + 1) * dx;	P2.y = y0 + j * dy;		P2.z = Z[(i + 1)*Ny + j];
			P3.x = x0 + (i + 1) * dx; P3.y = y0 + (j + 1) * dy;	P3.z = Z[(i + 1)*Ny + (j + 1)];
			P4.x = x0 + i * dx;		P4.y = y0 + (j + 1) * dy;	P4.z = Z[i*Ny + (j + 1)];

			Tar.P[n] = P1;	++n;
			Tar.P[n] = P2;	++n;
			Tar.P[n] = P3;	++n;
			tempv = (P2 - P1) % (P3 - P1);
			Tar.N[m] = tempv.normlize();
			Tar.S[m] = 0.5*tempv.norm();
			Tar.C[m] = (P1 + P2 + P3) / 3.0;
			Tar.L[m] = (P2 - P1).norm();
			if ((P3 - P2).norm() > Tar.L[m])
				Tar.L[m] = (P3 - P2).norm();
			if ((P3 - P1).norm() >  Tar.L[m])
				Tar.L[m] = (P3 - P1).norm();

			if (j == 0)
				Tar.E[m].x = -1;
			else
				Tar.E[m].x = m - 1;
			if (i == Nx - 2)
				Tar.E[m].y = -1;
			else
				Tar.E[m].y = m + 2 * (Ny - 1) + 1;
			Tar.E[m].z = m + 1;
			++m;

			Tar.P[n] = P1;	++n;
			Tar.P[n] = P3;	++n;
			Tar.P[n] = P4;	++n;
			tempv = (P3 - P1) % (P4 - P1);
			Tar.N[m] = tempv.normlize();
			Tar.S[m] = 0.5*tempv.norm();
			Tar.C[m] = (P1 + P3 + P4) / 3.0;
			Tar.L[m] = (P3 - P1).norm();
			if ((P4 - P3).norm() > Tar.L[m])
				Tar.L[m] = (P4 - P3).norm();
			if ((P4 - P1).norm() >  Tar.L[m])
				Tar.L[m] = (P4 - P1).norm();

			Tar.E[m].x = m - 1;
			if (j == Ny - 2)
				Tar.E[m].y = -1;
			else
				Tar.E[m].y = m + 1;
			if (i == 0)
				Tar.E[m].z = -1;
			else
				Tar.E[m].z = m - 2 * (Ny - 1) - 1;
			++m;
		}
	}
	delete[] Z; Z = nullptr;
	Tar.mytree.CreateNode(Tar.P, Tar.TrianglesID, Tar.TNum);
	NTEMP = 0;
	Tar.CreateDevTree(&(Tar.mytree), Tar.TreePointer);
	cout << "粗糙面三角面片数：" << Tar.TNum << endl;
}
