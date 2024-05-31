#ifndef SBR_RayTrace_H
#define SBR_RayTrace_H 
#include"MyClass.h"
bool RayTriangle(Ray &ray, Target &Tar, int ID);

bool RayBox(Point3f orig, Point3f dir, Point3f center, Point3f cube, double * len);


int RayTree1(Target &Tar, Ray &ray);


int RayTree2(Target &Tar, Ray &ray);


void RefN_CPU1(Target &Tar, Target &Tar1, Point3f AntennaPosition, vector<int> &Index);

void RefN_CPU12(Target &Tar1, Target &Tar2, bool bl, Point3f AntennaPosition, vector<int> &Index1, vector<int*> &Index2);

//���ε����������״ﾭ Tar1 �� Tar1�ٵ�Tar2 ������. bl��Tar1��Tar2�Ƿ���ͬһĿ�꣬true���ǡ�
void RefN_CPU112(Target &Tar1, Target &Tar2, bool bl, Point3f AntennaPosition, vector<int*> &Index11, vector<int*> &Index112);

//�Ĵε����������״ﾭ Tar1 �� Tar1�� Tar1�ٵ�Tar2 ������. bl��Tar1��Tar2�Ƿ���ͬһĿ�꣬true���ǡ�
void RefN_CPU1112(Target &Tar1, Target &Tar2, bool bl, Point3f AntennaPosition, vector<int*> &Index111, vector<int*> &Index1112);

//���ε����������״ﾭ Tar1 �� Tar2�ٵ�Tar1 ������ , Tar1 �� Tar2������ͬһ��Ŀ��
void RefN_CPU121(Target &Tar1, Target &Tar2, Point3f AntennaPosition, vector<int*> &Index12, vector<int*> &Index121);
#endif