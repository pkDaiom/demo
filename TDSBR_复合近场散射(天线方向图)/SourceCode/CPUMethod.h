#ifndef CPUMethod_H
#define CPUMethod_H 
#include"MyClass.h"
void AntennaPattern(Point3f ki, Point3c &EDir, vector<TDiracFunc<Complex>> &Amp);

Point3c CPUEi(double t, Point3f T_Radar_Position);

void GetSBRDistance1(Target &Tar, Point3f T_Radar_Position, Point3f R_Radar_Position, vector<int> &Index, double &dmin, double &dmax);

void GetSBRDistance12(Target &Tar1, Target &Tar2, Point3f T_Radar_Position, Point3f R_Radar_Position, vector<int*> &Index12, double &dmin, double &dmax);

void GetSBRDistance112(Target &Tar1, Target &Tar2, Point3f T_Radar_Position, Point3f R_Radar_Position, vector<int*> &Index112, double &dmin, double &dmax);

void GetSBRDistance1112(Target &Tar1, Target &Tar2, Point3f T_Radar_Position, Point3f R_Radar_Position, vector<int*> &Index1112, double &dmin, double &dmax);

void GetSBRDistance121(Target &Tar1, Target &Tar2, Point3f T_Radar_Position, Point3f R_Radar_Position, vector<int*> &Index121, double &dmin, double &dmax);

void GetSBRTimeSpane(Target &Tar1, Target &Tar2, Point3f T_Radar_Position, Point3f R_Radar_Position,
	vector<int> &Index1, vector<int*> &Index11, vector<int*> &Index111,
	vector<int> &Index2, vector<int*> &Index12, vector<int*> &Index112, vector<int*> &Index1112,
	vector<int*> &Index21, vector<int*> &Index212, double &timebegin, double &timeend);

void GOField(Complex er, Complex mr, Point3f N, Point3f Center1, Point3f Center2, AntennaCurrentMoment &TPolarization, Point3f T_Radar_Position, AntennaCurrentMoment &TPolarization1, Point3f &T_Radar_Position1);

Point3c MECAField(Target &Tar, int ID, double t, AntennaCurrentMoment &CurrentMoment, Point3f T_Radar_Position, Point3f R_Radar_Position);

//计算目标的单次散射
Point3c CPUSBREs1(Target &Tar, vector<int>&Index, double t, Point3f T_Radar_Position, Point3f R_Radar_Position);

//计算二次散射
Point3c CPUSBREs12(Target &Tar1, Target &Tar2, vector<int*>&Index, double t, Point3f T_Radar_Position, Point3f R_Radar_Position);

// 计算Tar1到Tar1再到Tar2三次散射
Point3c CPUSBREs112(Target &Tar1, Target &Tar2, vector<int*>&Index112, double t, Point3f T_Radar_Position, Point3f R_Radar_Position);

/// 计算Tar1到Tar11到Tar1再到Tar2四次散射
Point3c CPUSBREs1112(Target &Tar1, Target &Tar2, vector<int*>&Index1112, double t, Point3f T_Radar_Position, Point3f R_Radar_Position);

// 计算Tar1到Tar2再到Tar1三次散射
Point3c CPUSBREs121(Target &Tar1, Target &Tar2, vector<int*>&Index121, double t, Point3f T_Radar_Position, Point3f R_Radar_Position);

#endif
