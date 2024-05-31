#include "Variables.h"
#include "GenerateTarget.h"
void GetPTNum(string filename, int &TNum, int &PNum)
{
	fstream f2;
	int P = 0, T = 0;
	int i;
	string temp;
	char a;
	f2.open(filename, ios::in);
	for (i = 0; i<6; ++i)
	{
		getline(f2, temp);
	}
	f2 >> temp >> temp >> T;
	while (1 == 1)
	{
		f2.get(a);
		if (a == '$')
		{
			getline(f2, temp);
			continue;
		}
		if (a == 'G')
		{
			++P;
			getline(f2, temp);
			getline(f2, temp);
			continue;
		}
		if (a == 'C')
			break;
	}
	f2.close();
	TNum = T;
	PNum = P;
}

void ReadNas(Target &Tar, int PNum, string filename)
{
	int TNum = Tar.TNum;
	Point3f *Point = new Point3f[PNum];
	Point3i *Tri = new Point3i[TNum];//nas文件中每个三角形的三个顶点的编号

	fstream f2;
	string temp;
	char a;
	int tempi;
	f2.open(filename, ios::in);
	int i = 0, j = 0;
	while (1 == 1)
	{
		f2.get(a);
		if (a == '$')
		{
			getline(f2, temp);
			continue;
		}
		if (a == 'G')
		{
			f2 >> temp >> tempi >> Point[i].x >> Point[i].y >> tempi >> temp >> tempi >> Point[i].z; ++i;
			getline(f2, temp);
			continue;
		}
		if (a == 'C')
		{
			f2 >> temp >> tempi >> tempi >> Tri[j].x >> Tri[j].y >> Tri[j].z;

			getline(f2, temp);
			++j;
			continue;
		}
		if (a == 'E')
			break;
	}
	f2.close();


	for (i = 0; i<TNum; ++i)
	{
		Tri[i].x = Tri[i].x - 1;
		Tri[i].y = Tri[i].y - 1;
		Tri[i].z = Tri[i].z - 1;
	}

	for (i = 0; i<TNum; ++i)
	{
		//从0开始，每三个相邻点为同一个三角形的三个顶点
		Tar.P[3 * i] = Point[Tri[i].x];
		Tar.P[3 * i + 1] = Point[Tri[i].y];
		Tar.P[3 * i + 2] = Point[Tri[i].z];
	}

	double tempd;
	Point3f tempv;
	for (i = 0; i<TNum; ++i)
	{
		tempv = (Point[Tri[i].y] - Point[Tri[i].x]) % (Point[Tri[i].z] - Point[Tri[i].x]);
		tempd = tempv.norm();
		Tar.N[i] = tempv / tempd;
		Tar.S[i] = 0.5f * tempd;
		Tar.C[i] = (Point[Tri[i].x] + Point[Tri[i].y] + Point[Tri[i].z]) / 3.0f;
		Tar.L[i] = (Point[Tri[i].y] - Point[Tri[i].x]).norm();
		if ((Point[Tri[i].z] - Point[Tri[i].y]).norm() > Tar.L[i])
			Tar.L[i] = (Point[Tri[i].z] - Point[Tri[i].y]).norm();
		if ((Point[Tri[i].x] - Point[Tri[i].z]).norm() >  Tar.L[i])
			Tar.L[i] = (Point[Tri[i].x] - Point[Tri[i].z]).norm();
	}

	for (i = 0; i<TNum; ++i)
	{
		Tar.E[i].x = -1;
		Tar.E[i].y = -1;
		Tar.E[i].z = -1;
	}
	//找每一个三角形的相邻三个三角形
	int wedgeid = 0;
	for (i = 0; i<TNum; ++i)
	{
		for (j = i + 1; j<TNum; ++j)
		{
			//第i个三角形的第一条边
			if (Tar.E[i].x == -1)
			{
				if (Tri[i].x == Tri[j].y && Tri[i].y == Tri[j].x)
				{
					Tar.E[i].x = j;
					Tar.E[j].x = i;

				}
				else if (Tri[i].x == Tri[j].z && Tri[i].y == Tri[j].y)
				{
					Tar.E[i].x = j;
					Tar.E[j].y = i;

				}
				else if (Tri[i].x == Tri[j].x && Tri[i].y == Tri[j].z)
				{
					Tar.E[i].x = j;
					Tar.E[j].z = i;

				}

			}

			//第i个三角形的第二条边
			if (Tar.E[i].y == -1)
			{
				if (Tri[i].y == Tri[j].y && Tri[i].z == Tri[j].x)
				{
					Tar.E[i].y = j;
					Tar.E[j].x = i;

				}
				else if (Tri[i].y == Tri[j].z && Tri[i].z == Tri[j].y)
				{
					Tar.E[i].y = j;
					Tar.E[j].y = i;

				}
				else if (Tri[i].y == Tri[j].x && Tri[i].z == Tri[j].z)
				{
					Tar.E[i].y = j;
					Tar.E[j].z = i;

				}

			}
			//第i个三角形的第三条边
			if (Tar.E[i].z == -1)
			{
				if (Tri[i].z == Tri[j].y && Tri[i].x == Tri[j].x)
				{
					Tar.E[i].z = j;
					Tar.E[j].x = i;

				}
				else if (Tri[i].z == Tri[j].z && Tri[i].x == Tri[j].y)
				{
					Tar.E[i].z = j;
					Tar.E[j].y = i;

				}
				else if (Tri[i].z == Tri[j].x && Tri[i].x == Tri[j].z)
				{
					Tar.E[i].z = j;
					Tar.E[j].z = i;

				}

			}
		}

	}

	delete[] Point; Point = nullptr;
	delete[] Tri;	Tri = nullptr;
	return;
}

//读目标文件
void GenerateTarget(Target &Tar, Complex er, Complex mr, string filename)
{
	int PNum;
	Tar.er = Tar_er;
	Tar.mr = Tar_mr;
	GetPTNum(filename, Tar.TNum, PNum);
	Tar.Create();
	ReadNas(Tar, PNum, filename);
	cout << "目标三角面片数：" << Tar.TNum << endl;
	Tar.mytree.CreateNode(Tar.P, Tar.TrianglesID, Tar.TNum);
	NTEMP = 0;
	Tar.CreateDevTree(&(Tar.mytree), Tar.TreePointer);
}
