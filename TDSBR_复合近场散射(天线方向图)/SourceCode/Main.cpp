#include "Function.h"
#include "GenerateRoughSurface.h"
#include "GenerateTarget.h"
#include "Variables.h"
#include "CPUSBRField.h"
using namespace std;
int main()
{
	cout.precision(8);
	f.precision(8);
	LARGE_INTEGER BegainTime;
	LARGE_INTEGER EndTime;
	LARGE_INTEGER Frequency;
	LONGLONG cputime;
	QueryPerformanceFrequency(&Frequency);
	QueryPerformanceCounter(&BegainTime);

	Target Tar;
	Target TarBackgrond;
	Result EsResult;
	Radar TRadar, RRadar;
	RoughSurface roughsurface;

	SetRoughSurface1Para(roughsurface);


	if (computkind == singletarget || computkind == combination)
		GenerateTarget(Tar, Tar_er, Tar_mr, targetfile);

	if (computkind == singlesea)
		GenerateRoughSurface(Tar, roughsurface);

	if (computkind == combination)
	{
		GenerateRoughSurface(TarBackgrond, roughsurface);
	}

	cout << "��ʼ����......" << endl;
	CPUSBRField(Tar, TarBackgrond, TRadar, RRadar, EsResult);
	cout << " ----�������" << endl << endl;
	QueryPerformanceCounter(&EndTime);
	cputime = (EndTime.QuadPart - BegainTime.QuadPart) * 1000 / Frequency.QuadPart;
	cout << "��ʱ��" << cputime << "����" << endl << endl;

	cout << "��ʼ������..." << endl;
	f.open(Resultfile, ios::out);
	f << cputime << "ms" << endl;
	f << "t " << "Pt " << "t " << "�ܳ� " << "Ŀ�� " << "���� " << "��� " << endl;
	double Es, Es1, Es2, Es12;
	for (int i = 0; i < EsResult.ThNum; i++)
	{
		for (int j = 0; j < EsResult.PhNum; j++)
		{
			for (int k = 0; k < EsResult.RNum; k++)
			{
				for (int n = 0; n < EsResult.tNum[i][j][k]; n++)
				{
					Es = EsResult.Es[i][j][k][n].abs()* Sgn(EsResult.Es[i][j][k][n].x);
					Es1 = EsResult.Es1[i][j][k][n].abs()* Sgn(EsResult.Es1[i][j][k][n].x);
					Es2 = EsResult.Es2[i][j][k][n].abs()* Sgn(EsResult.Es2[i][j][k][n].x);
					Es12 = EsResult.Es12[i][j][k][n].abs()* Sgn(EsResult.Es12[i][j][k][n].x);
					f << EsResult.t[i][j][k][n] - EsResult.t[i][j][k][0] << ' ' << EsResult.Pt[i][j][k][n] << ' ';
					f << EsResult.t[i][j][k][n] << ' ' << Es << ' ' << Es1 << ' ' << Es2 << ' ' << Es12 << endl;
				}
			}
		}
	}
	f.close();

	cout << "������" << endl;
	Tar.Release();
	TarBackgrond.Release();
	EsResult.Release();
	TRadar.Release();
	RRadar.Release();
	return 0;
}
