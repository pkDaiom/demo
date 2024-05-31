#ifndef SBR_RayTrace_H
#define SBR_RayTrace_H 
#include"MyClass.h"
bool RayTriangle(Ray &ray, Target &Tar, int ID);

bool RayBox(Point3f orig, Point3f dir, Point3f center, Point3f cube, double * len);


int RayTree1(Target &Tar, Ray &ray);


int RayTree2(Target &Tar, Ray &ray);


void RefN_CPU1(Target &Tar, Target &Tar1, Point3f AntennaPosition, vector<int> &Index);

void RefN_CPU12(Target &Tar1, Target &Tar2, bool bl, Point3f AntennaPosition, vector<int> &Index1, vector<int*> &Index2);

//三次弹跳，计算雷达经 Tar1 到 Tar1再到Tar2 的射线. bl：Tar1和Tar2是否是同一目标，true：是。
void RefN_CPU112(Target &Tar1, Target &Tar2, bool bl, Point3f AntennaPosition, vector<int*> &Index11, vector<int*> &Index112);

//四次弹跳，计算雷达经 Tar1 到 Tar1到 Tar1再到Tar2 的射线. bl：Tar1和Tar2是否是同一目标，true：是。
void RefN_CPU1112(Target &Tar1, Target &Tar2, bool bl, Point3f AntennaPosition, vector<int*> &Index111, vector<int*> &Index1112);

//三次弹跳，计算雷达经 Tar1 到 Tar2再到Tar1 的射线 , Tar1 和 Tar2不能是同一个目标
void RefN_CPU121(Target &Tar1, Target &Tar2, Point3f AntennaPosition, vector<int*> &Index12, vector<int*> &Index121);
#endif