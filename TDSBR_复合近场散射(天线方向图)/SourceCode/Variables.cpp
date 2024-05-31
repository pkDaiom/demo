#include "Variables.h"
/*******************************����***********************************/
double pi = acos(-1.0);
double G = 9.8;
double EPS = 1E-16;
double EPB = 1E16;
Complex J(0.0, 1.0);
double eps0 = 8.85418781761e-12;
double mur0 = pi*4.0e-7;
double eta0 = sqrt(mur0 / eps0);
double c0 = 1.0 / sqrt(mur0*eps0)*1.0e-9;
Complex vacuum_er = Complex(1, 0);
Complex vacuum_mr = Complex(1, 0);
double K = 1;

/*****************************�˲���**********************************/
int MaxEachOrderTrianglesNum = 500;
int MaxOrderNum = 15;
int NTEMP = 0;
int CurrentNode = 0;

/*****************����***************************/
ComputKind computkind = singletarget;

string Resultfile = "TDResult.dat";
string SolutionType = "��վ";
fstream f;

/****************************����***********************************************/
double f0 = 8;
double fband = 4;
double tao = 4. / fband;
double t0 = 1.0*tao;
double dt = 1. / (4.*f0);
double tBegin = 0;
double tEnd = 0;

/****************************��������******************************************/
Point3f TBeamDirection;     //������������������
Point3f TPolarization;       //�������ߵ�����������
AntennaType antennatype = LinearAntenna;
string T_Radar_Polarization = "H"; //����������ʽ

double T_Radar_Theta1 = 0.0;
double T_Radar_Theta2 = 0.0;
double T_Radar_dTheta = 1;//����Ϊ0

double T_Radar_Phi1 = 30;
double T_Radar_Phi2 = 30;
double T_Radar_dPhi = 0.2;//����Ϊ0

double T_Radar_R1 = 30;
double T_Radar_R2 = 30;
double T_Radar_dR = 0.1;//����Ϊ0

						/****************************��������******************************************/
Point3f RBeamDirection;     //������������������
Point3f RPolarization; //�������߼�������

string R_Radar_Polarization = "H";//�������߼�����ʽ

double R_Radar_Theta1 = 0;
double R_Radar_Theta2 = 90;
double R_Radar_dTheta = 1;//����Ϊ0

double R_Radar_Phi1 = 0;
double R_Radar_Phi2 = 0;
double R_Radar_dPhi = 0.2; //����Ϊ0

double R_Radar_R1 = 1;
double R_Radar_R2 = 1;
double R_Radar_dR = 0.1;//����Ϊ0

						/****************************�ֲ���****************************/
string RoughS1file = "surface.nas";
Complex RoughS1_er = Complex(1, -1e9);
Complex RoughS1_mr = Complex(1, 0);
double RoughS1_High = 0.0;

void SetRoughSurface1Para(RoughSurface &roughsurface1)
{
	double lamda = c0 / f0;
	roughsurface1.name = "surface.dat";
	roughsurface1.Spectrum = Gauss;
	roughsurface1.Seed = 12356.0;
	roughsurface1.LNx = 320;
	roughsurface1.LNy = 200;
	roughsurface1.Lx = (roughsurface1.LNx - 1)*lamda / 8.0;
	roughsurface1.Ly = (roughsurface1.LNy - 1)*lamda / 8.0;
	roughsurface1.High = RoughS1_High;

	roughsurface1.delta = 0.5*lamda;
	roughsurface1.rlx = 1.5*lamda;
	roughsurface1.rly = 1.5*lamda;
	roughsurface1.er = RoughS1_er;
	roughsurface1.mr = RoughS1_mr;

	roughsurface1.WindSpeed = 10.0;
	roughsurface1.WindTheta = 0.0;
	roughsurface1.Temp = 20;
	roughsurface1.Sault = 34.7;
	roughsurface1.DeltT = 0;//�����Ų�����Ƶ����ʱȡ����FreqBegin����������ɢ��ʱ����Ҫ�޸ĳ���
}

/****************************Ŀ��***********************************************/
Complex Tar_er = Complex(1, -1e9);
Complex Tar_mr = Complex(1, 0);
string targetfile = "chuan15m_1.8m_3m.nas";