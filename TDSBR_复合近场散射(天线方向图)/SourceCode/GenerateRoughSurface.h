#ifndef GenerateRoughSurface_H
#define GenerateRoughSurface_H 
#include"MyClass.h"
//��˹��
void Gauss_Spm(double *S, RoughSurface &roughsurface);

void PM_Spm(double *S, RoughSurface &sea);

double C(double k);

void ELH_Spm(double *S, RoughSurface &sea);

// FΪƵ�ʣ�TΪ�¶ȣ�S_sw���ζȣ���������ƽ���ζ�Ϊ35�룬����S_swȡ34.7��
void SeaDielec(double F, RoughSurface &sea);

void GenerateRoughSurface(Target &Tar, RoughSurface &roughsurface);
#endif