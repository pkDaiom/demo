#ifndef GenerateRoughSurface_H
#define GenerateRoughSurface_H 
#include"MyClass.h"
//高斯谱
void Gauss_Spm(double *S, RoughSurface &roughsurface);

void PM_Spm(double *S, RoughSurface &sea);

double C(double k);

void ELH_Spm(double *S, RoughSurface &sea);

// F为频率，T为温度，S_sw是盐度，世界大洋的平均盐度为35‰，这里S_sw取34.7。
void SeaDielec(double F, RoughSurface &sea);

void GenerateRoughSurface(Target &Tar, RoughSurface &roughsurface);
#endif