#ifndef VARIABLES_H
#define VARIABLES_H 1
#include "MyClass.h"

/*******************************����***********************************/
extern double pi;
extern double G;
extern Complex J;
extern double eps0;
extern double mur0;
extern double eta0;
extern double c0;
extern Complex vacuum_er;
extern Complex vacuum_mr;
extern double K;


/*****************����***************************/
extern ComputKind computkind;

extern string Resultfile;
extern string SolutionType;
extern fstream f;

/****************************����***********************************************/
extern double f0;
extern double fband;
extern double tao;
extern double t0;
extern double dt;
extern double tBegin;
extern double tEnd;

/****************************��������******************************************/
extern Point3f TBeamDirection;     //������������������
extern Point3f TPolarization;       //�������ߵ�����������
extern AntennaType antennatype;
extern string T_Radar_Polarization;

extern double T_Radar_Theta1;
extern double T_Radar_Theta2;
extern double T_Radar_dTheta;

extern double T_Radar_Phi1;
extern double T_Radar_Phi2;
extern double T_Radar_dPhi;

extern double T_Radar_R1;
extern double T_Radar_R2;
extern double T_Radar_dR;

/****************************��������******************************************/
extern Point3f RBeamDirection;     //������������������
extern Point3f RPolarization; //�������߼�������

extern string R_Radar_Polarization;

extern double R_Radar_Theta1;
extern double R_Radar_Theta2;
extern double R_Radar_dTheta;

extern double R_Radar_Phi1;
extern double R_Radar_Phi2;
extern double R_Radar_dPhi;

extern double R_Radar_R1;
extern double R_Radar_R2;
extern double R_Radar_dR;

/****************************�ֲ���****************************/
extern string RoughS1file;
extern Complex RoughS1_er;
extern Complex RoughS1_mr;
extern double RoughS1_High;

void SetRoughSurface1Para(RoughSurface &roughsurface1);

/****************************Ŀ��***********************************************/
extern Complex Tar_er;
extern Complex Tar_mr;
extern string targetfile;
#endif