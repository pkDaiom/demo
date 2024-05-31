#include "Function.h"
#include "Variables.h"
#include "SBR_RayTrace.h"
#include "CPUMethod.h"
#include "CPUSBRField.h"

void SwapVector(vector<int *> &V)
{
	vector<int *>::iterator it;
	for (it = V.begin(); it != V.end(); it++)
		if (nullptr != *it)
		{
			delete[] * it;
			*it = nullptr;
		}
	vector<int*>().swap(V);
}

//计算散射场
void CPUSBRField(Target &Tar, Target &TarSea, Radar &TRadar, Radar &RRadar, Result &EsResult)
{

	int i, j, k, n, m1, m2, m3;
	int iThNum, iPhNum, iRNum;
	int sThNum, sPhNum, sRNum;
	double thetai, phii, T_Radar_R;
	double thetas, phis, R_Radar_R;
	double t;
	Point3f T_Radar_Position, R_Radar_Position;

	Complex Eshh, Esvv;

	vector<int> Index1;           //记录到舰船的射线编号
	vector<int*> Index11;      //记录到舰船到舰船的射线编号 舰船考虑2次射线跟踪
	vector<int*> Index111;      //记录到舰船到舰船的射线编号 舰船考虑3次射线跟踪

	vector<int> Index2;				//记录到海面的射线编号

	vector<int*> Index12;     //记录到舰船到海面的射线编号
	vector<int*> Index112;     //记录到舰船到舰船到海面的射线编号
	vector<int*> Index1112;     //记录到舰船到舰船到舰船到海面的射线编号

	vector<int*> Index21;			//记录到海面到舰船的射线编号
	vector<int*> Index212;			//记录到海面到舰船到海面的射线编号

	TRadar.Creat(T_Radar_Theta1, T_Radar_Theta2, T_Radar_dTheta,
		T_Radar_Phi1, T_Radar_Phi2, T_Radar_dPhi,
		T_Radar_R1, T_Radar_R2, T_Radar_dR);

	if (SolutionType == "双站")
	{
		RRadar.Creat(R_Radar_Theta1, R_Radar_Theta2, R_Radar_dTheta,
			R_Radar_Phi1, R_Radar_Phi2, R_Radar_dPhi,
			R_Radar_R1, R_Radar_R2, R_Radar_dR);
		sThNum = RRadar.ThNum;
		sPhNum = RRadar.PhNum;
		sRNum = RRadar.RNum;
	}

	iThNum = TRadar.ThNum;
	iPhNum = TRadar.PhNum;
	iRNum = TRadar.RNum;

	if (SolutionType == "单站")
	{
		EsResult.Creat(iThNum, iPhNum, iRNum);
	}

	if (SolutionType == "双站")
	{
		EsResult.Creat(sThNum, sPhNum, sRNum);
	}

	int TotalProgressBarnum = 0;
	ProgressBar(100, TotalProgressBarnum, 100, 0);

	for (i = 0; i < iThNum; ++i)
	{

		thetai = TRadar.Theta[i];
		thetai = thetai* pi / 180.0;
		for (j = 0; j < iPhNum; ++j)
		{
			phii = TRadar.Phi[j];
			phii = phii* pi / 180.0;

			TBeamDirection = Point3f(-sin(thetai)*cos(phii), -sin(thetai)*sin(phii), -cos(thetai)).normlize();

			if (T_Radar_Polarization == "H")
				TPolarization = Point3f(-sin(phii), cos(phii), 0);//H极化
			else
				TPolarization = (TBeamDirection%Point3f(-sin(phii), cos(phii), 0)).normlize();//V极化

			for (k = 0; k < iRNum; ++k)
			{
				T_Radar_R = TRadar.R[k];
				T_Radar_Position = T_Radar_R*(-TBeamDirection);

				RefN_CPU1(Tar, TarSea, T_Radar_Position, Index1);
				RefN_CPU12(Tar, Tar, true, T_Radar_Position, Index1, Index11);
				RefN_CPU112(Tar, Tar, true, T_Radar_Position, Index11, Index111);

				if (computkind == combination)
				{
					RefN_CPU1(TarSea, Tar, T_Radar_Position, Index2);

					RefN_CPU12(Tar, TarSea, false, T_Radar_Position, Index1, Index12);//记录到舰船到海面的射线编号
					RefN_CPU112(Tar, TarSea, false, T_Radar_Position, Index11, Index112);//记录到舰船到舰船到海面的射线编号
					RefN_CPU1112(Tar, TarSea, false, T_Radar_Position, Index111, Index1112);//记录到舰船到舰船到舰船到海面的射线编号

					RefN_CPU12(TarSea, Tar, false, T_Radar_Position, Index2, Index21); //记录到海面到舰船的射线编号
					RefN_CPU121(TarSea, Tar, T_Radar_Position, Index21, Index212); //记录到海面到舰船到海面的射线编号
				}

				/*以上代码都是构建入射条件 */

				if (SolutionType == "单站")
				{
					RPolarization = TPolarization;
					R_Radar_Position = T_Radar_Position;

					GetSBRTimeSpane(Tar, TarSea, T_Radar_Position, R_Radar_Position,
						Index1, Index11, Index111, Index2, Index12, Index112, Index1112, Index21, Index212, tBegin, tEnd);

					EsResult.CreatTime(i, j, k, tBegin, tEnd, dt);

					for (n = 0; n < EsResult.tNum[i][j][k]; n++)
					{
						t = EsResult.t[i][j][k][n];
						EsResult.Pt[i][j][k][n] = GaossPulse(t - EsResult.t[i][j][k][0], 0);
						EsResult.Ei[i][j][k][n] = CPUEi(t - EsResult.t[i][j][k][0], T_Radar_Position);

						//计算舰船的Es
						EsResult.Es1[i][j][k][n] = (CPUSBREs1(Tar, Index1, t, T_Radar_Position, R_Radar_Position)
							+ CPUSBREs12(Tar, Tar, Index11, t, T_Radar_Position, R_Radar_Position)
							+ CPUSBREs112(Tar, Tar, Index111, t, T_Radar_Position, R_Radar_Position)
							)* RPolarization*eta0;

						//计算单独海面的Es
						EsResult.Es2[i][j][k][n] = CPUSBREs1(TarSea, Index2, t, T_Radar_Position, R_Radar_Position)* RPolarization*eta0;

						// 计算舰船到海面的Es
						EsResult.Es12[i][j][k][n] = (CPUSBREs12(Tar, TarSea, Index12, t, T_Radar_Position, R_Radar_Position)
							+ CPUSBREs112(Tar, TarSea, Index112, t, T_Radar_Position, R_Radar_Position)
							+ CPUSBREs1112(Tar, TarSea, Index1112, t, T_Radar_Position, R_Radar_Position)
							)* RPolarization*eta0;

						// 计算海面到舰船的Es
						EsResult.Es21[i][j][k][n] = CPUSBREs12(TarSea, Tar, Index21, t, T_Radar_Position, R_Radar_Position) * RPolarization*eta0;

						// 计算海面到舰船到海面的Es
						EsResult.Es212[i][j][k][n] = CPUSBREs121(TarSea, Tar, Index212, t, T_Radar_Position, R_Radar_Position) * RPolarization*eta0;

						//用Es12代表耦合场
						EsResult.Es12[i][j][k][n] = EsResult.Es12[i][j][k][n] + EsResult.Es21[i][j][k][n] + EsResult.Es212[i][j][k][n];

						EsResult.Es1[i][j][k][n] = EsResult.Es1[i][j][k][n] + EsResult.EsPTD[i][j][k][n];
						EsResult.Es[i][j][k][n] = EsResult.Es1[i][j][k][n] + EsResult.Es2[i][j][k][n] + EsResult.Es12[i][j][k][n];

						ProgressBar(iThNum*iPhNum*iRNum, TotalProgressBarnum, EsResult.tNum[i][j][k], n + 1);
					}

					ProgressBar(iThNum*iPhNum*iRNum, ++TotalProgressBarnum, 100, 0);
				}
				if (SolutionType == "双站")
				{
					for (m1 = 0; m1 < sThNum; ++m1)
					{
						thetas = RRadar.Theta[m1];
						thetas = thetas* pi / 180.0f;
						for (m2 = 0; m2 < sPhNum; ++m2)
						{
							phis = RRadar.Phi[m2];
							phis = phis* pi / 180.0f;

							RBeamDirection = Point3f(sin(thetas) * cos(phis), sin(thetas) * sin(phis), cos(thetas)).normlize();

							if (R_Radar_Polarization == "H")
								RPolarization = Point3f(-sin(phis), cos(phis), 0);
							else
								RPolarization = (RBeamDirection % Point3f(-sin(phis), cos(phis), 0)).normlize();

							for (m3 = 0; m3 < sRNum; m3++)
							{
								R_Radar_R = RRadar.R[m3];
								R_Radar_Position = R_Radar_R*RBeamDirection;

								GetSBRTimeSpane(Tar, TarSea, T_Radar_Position, R_Radar_Position,
									Index1, Index11, Index111, Index2, Index12, Index112, Index1112, Index21, Index212, tBegin, tEnd);

								EsResult.CreatTime(m1, m2, m3, tBegin, tEnd, dt);

								for (n = 0; n < EsResult.tNum[m1][m2][m3]; n++)
								{
									t = EsResult.t[m1][m2][m3][n];
									EsResult.Pt[m1][m2][m3][n] = GaossPulse(t - EsResult.t[m1][m2][m3][0], 0);

									EsResult.Es[m1][m2][m3][n] = 0.0;   EsResult.Es1[m1][m2][m3][n] = 0.0;    EsResult.EsPTD[m1][m2][m3][n] = 0.0;
									EsResult.Es2[m1][m2][m3][n] = 0.0;   EsResult.Es21[m1][m2][m3][n] = 0.0;    EsResult.Es12[m1][m2][m3][n] = 0.0; EsResult.Es212[m1][m2][m3][n] = 0.0;

									EsResult.Ei[m1][m2][m3][n] = CPUEi(t - EsResult.t[m1][m2][m3][0], T_Radar_Position);

									//计算舰船的Es
									EsResult.Es1[m1][m2][m3][n] = (CPUSBREs1(Tar, Index1, t, T_Radar_Position, R_Radar_Position)
										+ CPUSBREs12(Tar, Tar, Index11, t, T_Radar_Position, R_Radar_Position)
										+ CPUSBREs112(Tar, Tar, Index111, t, T_Radar_Position, R_Radar_Position)
										)* RPolarization*eta0;

									//计算单独海面的Es
									EsResult.Es2[m1][m2][m3][n] = CPUSBREs1(TarSea, Index2, t, T_Radar_Position, R_Radar_Position)* RPolarization*eta0;

									// 计算舰船到海面的Es
									EsResult.Es12[m1][m2][m3][n] = (CPUSBREs12(Tar, TarSea, Index12, t, T_Radar_Position, R_Radar_Position)
										+ CPUSBREs112(Tar, TarSea, Index112, t, T_Radar_Position, R_Radar_Position)
										+ CPUSBREs1112(Tar, TarSea, Index1112, t, T_Radar_Position, R_Radar_Position))* RPolarization*eta0;

									// 计算海面到舰船的Es
									EsResult.Es21[m1][m2][m3][n] = CPUSBREs12(TarSea, Tar, Index21, t, T_Radar_Position, R_Radar_Position)* RPolarization*eta0;

									// 计算海面到舰船再到海面的Es
									EsResult.Es212[m1][m2][m3][n] = CPUSBREs121(TarSea, Tar, Index212, t, T_Radar_Position, R_Radar_Position)* RPolarization*eta0;

									//用Es12代表耦合场
									EsResult.Es12[m1][m2][m3][n] = EsResult.Es12[m1][m2][m3][n] + EsResult.Es21[m1][m2][m3][n] + EsResult.Es212[m1][m2][m3][n];

									EsResult.Es1[m1][m2][m3][n] = EsResult.Es1[m1][m2][m3][n] + EsResult.EsPTD[m1][m2][m3][n];
									EsResult.Es[m1][m2][m3][n] = EsResult.Es1[m1][m2][m3][n] + EsResult.Es2[m1][m2][m3][n] + EsResult.Es12[m1][m2][m3][n];

									ProgressBar(iThNum*iPhNum*iRNum*sThNum*sPhNum*sRNum, TotalProgressBarnum, EsResult.tNum[m1][m2][m3], n + 1);
								}
								ProgressBar(iThNum*iPhNum*iRNum*sThNum*sPhNum*sRNum, ++TotalProgressBarnum, 100, 0);
							}
						}
					}
				}
				vector<int>().swap(Index1);
				SwapVector(Index11);
				SwapVector(Index111);
				vector<int>().swap(Index2);
				SwapVector(Index12);
				SwapVector(Index112);
				SwapVector(Index1112);
				SwapVector(Index21);
				SwapVector(Index212);
			}
		}
	}
}