#ifndef GenerateTarget_H
#define GenerateTarget_H 
#include"MyClass.h"
void GetPTNum(string filename, int &TNum, int &PNum);

void ReadNas(Target &Tar, int PNum, string filename);

//��Ŀ���ļ�
void GenerateTarget(Target &Tar, Complex er, Complex mr, string filename);
#endif