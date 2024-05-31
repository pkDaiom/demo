#ifndef MyClass_H
#define MyClass_H 
#include <math.h>
#include <memory.h>
#include <set>
#include<vector>
#include <windows.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
using namespace std;

extern int MaxEachOrderTrianglesNum;
extern int MaxOrderNum;
extern int CurrentNode;
extern int NTEMP;
extern double EPS;
extern double EPB;


#ifndef _COMPLEX_
#define _COMPLEX_
class Complex{
public:
	double x;
	double y;
	Complex(){x = 0.0; y = 0.0;}
	Complex(double _x, double _y):x(_x),y(_y) {}
	// 一元操作符重载	
	inline Complex operator +(){return *this;}
	inline Complex operator -(){return Complex(-x,-y);}
	inline Complex operator = ( Complex p ){x=p.x;y=p.y;return *this;}
	inline Complex operator = ( double p ){x=p;y=0.0f;return *this;}
	inline bool	operator   ==( Complex p ){return x==p.x&&y==p.y;}
	inline bool	operator   !=( Complex p ){return x!=p.x||y!=p.y;}
	// 求模
	 inline double abs() {return sqrt(x*x+y*y);}
	 //求共轭
	 inline Complex Conjg() { return Complex(x,-y); }
};
#endif

#ifndef _POINT3F_
#define _POINT3F_
class Point3f{
public:
	double x;
	double y;
	double z;
	Point3f(){x = 0.0;y=0.0;z=0.0;}
	Point3f(double _x, double _y, double _z):x(_x),y(_y),z(_z) {}
	// 一元操作符重载	
	inline Point3f operator +(){return *this;}
	inline Point3f operator -(){return Point3f(-x,-y,-z);}
	inline Point3f operator = ( Point3f p ){x=p.x;y=p.y;z=p.z;return *this;}
	inline bool	  operator ==( Point3f p ){return x==p.x&&y==p.y&&z==p.z;}
	inline bool	  operator !=( Point3f p ){return x!=p.x||y!=p.y||z!=p.z;}
	// 求模
	inline double norm() {return sqrt(x*x+y*y+z*z);}
	// 归一化
	inline Point3f normlize()
	{
		if (sqrt(x*x + y*y + z*z) < EPS)
			return Point3f(0.0, 0.0, 0.0);
		else
			return Point3f(x / sqrt(x*x + y*y + z*z), y / sqrt(x*x + y*y + z*z), z / sqrt(x*x + y*y + z*z));
	}
};
#endif

#ifndef _POINT3C_
#define _POINT3C_
class Point3c{
public:
	Complex x;
	Complex y;
	Complex z;

	// 构造函数
	 Point3c(){x.x =0.0;x.y=0.0;y.x=0.0;y.y=0.0;z.x=0.0;z.y=0.0;}
	 Point3c(Complex _x, Complex _y, Complex _z):x(_x),y(_y),z(_z) {}

	// 操作符重载
	 inline Point3c operator + (){return *this;}
	 inline Point3c operator - (){return Point3c(-x,-y,-z);}
	 inline Point3c	operator = ( Point3c p ){x=p.x;y=p.y;z=p.z;return *this;}
	 inline Point3c	operator = (Point3f p){ x = p.x; y = p.y; z = p.z; return *this; }
	 inline bool	operator ==( Point3c p) {return x==p.x&&y==p.y&&z==p.z;}
	 inline bool	operator !=( Point3c p) {return x!=p.x||y!=p.y||z!=p.z;}
	// 求模
	 double norm() { return sqrt(x.x * x.x + x.y*x.y + y.x*y.x+y.y*y.y + z.x*z.x+z.y*z.y);}
	 Point3f normlize() 
	{
		double rx = x.x*x.x + x.y*x.y;
		double ry = y.x*y.x + y.y*y.y;
		double rz = z.x*z.x + z.y*z.y;
		double r=sqrt(rx + ry + rz);
		if (r<1e-7)
			return Point3f(0,0,0);
		return Point3f(sqrt(rx / (rx + ry + rz)), sqrt(ry / (rx + ry + rz)), sqrt(rz / (rx + ry + rz)));
	}
};
#endif

#ifndef _MYOPERATORFUN_
#define _MYOPERATORFUN_

//Complex类的符号重载函数
inline Complex operator + (Complex p1,Complex p2)
{
	return Complex(p1.x+p2.x,p1.y+p2.y);
}
inline Complex operator + (Complex p1,double p2)
{
	return Complex(p1.x+p2,p1.y);
}
inline Complex operator + (double p1,Complex p2)
{
	return Complex(p1+p2.x,p2.y);
}

inline Complex operator - (Complex p1,Complex p2)
{
	return Complex(p1.x-p2.x,p1.y-p2.y);
}
inline Complex operator - (Complex p1,double p2)
{
	return Complex(p1.x-p2,p1.y);
}
inline Complex operator - (double p1,Complex p2)
{
	return Complex(p1-p2.x,-p2.y);
}

inline Complex operator * (Complex p1,Complex p2)
{
	return Complex(p1.x*p2.x - p1.y*p2.y, p1.y*p2.x + p1.x*p2.y);
}
inline Complex operator * (Complex p1,double p2)
{
	return Complex(p1.x*p2,p1.y*p2);
}
inline Complex operator * (double p1,Complex p2)
{
	return Complex(p1*p2.x,p1*p2.y);
}

inline Complex operator / (Complex p1,Complex p2)
{
	return Complex((p1.x*p2.x + p1.y*p2.y) / (p2.x*p2.x + p2.y*p2.y), (p1.y*p2.x - p1.x*p2.y) / (p2.x*p2.x + p2.y*p2.y));
}

inline Complex operator / (Complex p1,double p2)
{
	return Complex(p1.x/p2,p1.y/p2);
}
inline Complex operator / (double p1,Complex p2)
{
	return Complex(p1*p2.x / (p2.x*p2.x + p2.y*p2.y), -p1*p2.y / (p2.x*p2.x + p2.y*p2.y));
}

inline Complex powC(Complex p, int n)
{
	Complex p2;
	double q = atan2(p.y,p.x);
	double r = sqrt(p.x*p.x + p.y*p.y);
	if (r + 1.0 != 1.0)
	{
		r = n*log(r);
		r = exp(r);
	}
	p2.x = r*cos(n*q);
	p2.y = r*sin(n*q);
	return p2;
}

inline Complex powC(Complex p, double n)
{
	Complex p2;
	double q = atan2(p.y,p.x);
	double r = sqrt(p.x*p.x + p.y*p.y);
	if (r + 1.0 != 1.0)
	{
		r = n*log(r);
		r = exp(r);
	}
	p2.x = r*cos(n*q);
	p2.y = r*sin(n*q);
	return p2;
}

inline Complex sqrtC(Complex p, int n)
{
	Complex p2;
	double q = atan2(p.y,p.x)/n;
	double r = sqrt(p.x*p.x + p.y*p.y);
	if (r + 1.0 != 1.0)
	{
		r = (1.0/n)*log(r);
		r = exp(r);
	}
	p2.x = r*cos(q);
	p2.y = r*sin(q);
	return p2;
}

inline Complex expC(Complex p1)
{
	Complex p2;
	double p = exp(p1.x);
	p2.x = p * cos(p1.y);
	p2.y = p * sin(p1.y);
	return p2;
}

inline Complex logC(Complex p1)
{
	return Complex(log(sqrt(p1.x*p1.x+p1.y*p1.y)),atan2(p1.y,p1.x));
}

inline Complex sinC(Complex p1)
{
	Complex p2;
	double p = exp(p1.y),q = exp(-p1.y);
	p2.x = sin(p1.x) * (p+q) * 0.5;
	p2.y = cos(p1.x) * (p-q) * 0.5;
	return p2;
}

inline Complex cosC(Complex p1)
{
	Complex p2;
	double p = exp(p1.y),q = exp(-p1.y);
	p2.x = cos(p1.x) * (p+q) * 0.5;
	p2.y = -sin(p1.x) * (p-q) * 0.5;
	return p2;
}
//Point3f类、Point3c的符号重载函数

inline Point3f operator + (Point3f p1,Point3f p2)
{
	return Point3f(p1.x+p2.x,p1.y+p2.y,p1.z+p2.z);
}

inline Point3c operator + (Point3c p1,Point3c p2)
{
	return Point3c(p1.x+p2.x,p1.y+p2.y,p1.z+p2.z);
}

inline Point3c operator + (Point3f p1,Point3c p2)
{
	return Point3c(p1.x+p2.x,p1.y+p2.y,p1.z+p2.z);
}

inline Point3c operator + (Point3c p1,Point3f p2)
{
	return Point3c(p1.x+p2.x,p1.y+p2.y,p1.z+p2.z);
}

inline Point3f operator - (Point3f p1,Point3f p2)
{
	return Point3f(p1.x-p2.x,p1.y-p2.y,p1.z-p2.z);
}

inline Point3c operator - (Point3c p1,Point3c p2)
{
	return Point3c(p1.x-p2.x,p1.y-p2.y,p1.z-p2.z);
}

inline Point3c operator - (Point3f p1,Point3c p2)
{
	return Point3c(p1.x-p2.x,p1.y-p2.y,p1.z-p2.z);
}

inline Point3c operator - (Point3c p1,Point3f p2)
{
	return Point3c(p1.x-p2.x,p1.y-p2.y,p1.z-p2.z);
}

inline Point3f operator * (Point3f p,double f)
{
	return Point3f(f*p.x,f*p.y,f*p.z);
}
inline Point3f operator * (double f,Point3f p)
{
	return Point3f(f*p.x,f*p.y,f*p.z);
}
inline Point3c operator * (Complex f,Point3f p)
{
	return Point3c(f*p.x,f*p.y,f*p.z);
}
inline Point3c operator *(Point3f p, Complex f)
{
	return Point3c(f*p.x,f*p.y,f*p.z);
}
inline double  operator * (Point3f p1,Point3f p2)
{
	return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z;
}
inline Point3c operator * (Point3c p,double f)
{
	return Point3c(f*p.x,f*p.y,f*p.z);
}
inline Point3c operator * (double f,Point3c p)
{
	return Point3c(f*p.x,f*p.y,f*p.z);
}
inline Point3c operator * (Complex f,Point3c p)
{
	return Point3c(f*p.x,f*p.y,f*p.z);
}
inline Point3c operator *(Point3c p, Complex f)
{
	return Point3c(f*p.x,f*p.y,f*p.z);
}
inline Complex operator * (Point3c p1,Point3c p2)
{
	return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z;
}
inline Complex  operator * (Point3c p1,Point3f p2)
{
	return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z;
}
inline Complex  operator * (Point3f p1,Point3c p2)
{
	return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z;
}

inline Point3f operator / (Point3f p, double f)
{
	return Point3f(p.x/f,p.y/f,p.z/f);
}

inline Point3c operator / (Point3f p, Complex f)
{
	return Point3c(p.x / f, p.y / f, p.z / f);
}

inline Point3c operator / (Point3c p, double f)
{
	return Point3c(p.x/f,p.y/f,p.z/f);
}

inline Point3c operator / (Point3c p, Complex f)
{
	return Point3c(p.x / f, p.y / f, p.z / f);
}

inline Point3f operator % (Point3f p1,Point3f p2)
{
	return Point3f(p1.y*p2.z-p1.z*p2.y, p1.z*p2.x-p1.x*p2.z, p1.x*p2.y-p1.y*p2.x);
}
inline Point3c operator % (Point3c p1,Point3c p2)
{
	return Point3c(p1.y*p2.z-p1.z*p2.y, p1.z*p2.x-p1.x*p2.z, p1.x*p2.y-p1.y*p2.x);
}
inline Point3c operator % (Point3f p1,Point3c p2)
{
	return Point3c(p1.y*p2.z-p1.z*p2.y, p1.z*p2.x-p1.x*p2.z, p1.x*p2.y-p1.y*p2.x);
}
inline Point3c operator % (Point3c p1,Point3f p2)
{
	return Point3c(p1.y*p2.z-p1.z*p2.y, p1.z*p2.x-p1.x*p2.z, p1.x*p2.y-p1.y*p2.x);
}

#endif

#ifndef _POINT3I_
#define _POINT3I_
class Point3i{
public:
	int x,y,z;
};
#endif

#ifndef _POINT2I_
#define _POINT2I_
class Point2i{
public:
	int x,y;
};
#endif

#ifndef _MyENUM_
#define _MyENUM_
enum PowerSpectrum
{
	Gauss,
	PM,
	ELH
};
enum ComputKind
{
	singletarget,
	singlesea,
	combination
};
enum AntennaType
{
	DipoleAntenna,
	DipoleArray,
	LinearAntenna
};
#endif

#ifndef _AXIS_
#define _AXIS_
class Axis{
public:
	Point3f orig;
	Point3f dir;
};
#endif

#ifndef _TDiracFunc_
#define _TDiracFunc_
template <class TM>
class TDiracFunc
{
public:
	double Var0;
	int DOrder;
	TM Amp;
	TDiracFunc(double _Var0, int _DOrder, TM _Amp)
	{
		Var0 = _Var0;
		DOrder = _DOrder;
		Amp = _Amp;
	}
};
#endif

#ifndef _AntennaCurrentMoment_
#define _AntennaCurrentMoment_
struct AntennaCurrentMoment
{
	vector<TDiracFunc<Complex>> Amp;
	Point3c EDir;
};
#endif

#ifndef _KDTREE_
#define _KDTREE_
class Kdtree
{
public:
	int			m_NumT;		// Octree中的坐标点数目

	int			m_ifEndNode;			//	1 表示 最底层子节点，用于判断是否是最底层节点

	int			* m_TID;			// 一个数组，保存最底层树节点中三角面片的编号

	Point3f		m_Cube;				// m_vCube三维所标分别表示长宽高

	Point3f		m_Center;				// 立方体的中心

	Kdtree		*m_Left;

	Kdtree		*m_Right;

	Kdtree()
	{
		m_NumT = 0;
		m_ifEndNode = 0;
		m_TID = 0;
		m_Cube = Point3f(0, 0, 0);
		m_Center = Point3f(0, 0, 0);
		m_Left = 0;
		m_Right = 0;
	}

	void CreateNode(Point3f *Vertices, int *TrianglesID, int iNumTriangles)
	{
		int i;
		m_NumT = iNumTriangles;
		Point3f temp3f, vMinPoint = Vertices[3 * TrianglesID[0]], vMaxPoint = Vertices[3 * TrianglesID[0]];
		///vMinPoint和vMaxPoint是目标外接长方体对角线上的两个顶点
		for (i = 0; i<iNumTriangles; ++i)
		{
			for (int j = 0; j < 3; j++)
			{
				temp3f = Vertices[3 * TrianglesID[i] + j];
				if (temp3f.x < vMinPoint.x)	vMinPoint.x = temp3f.x;
				if (temp3f.y < vMinPoint.y)	vMinPoint.y = temp3f.y;
				if (temp3f.z < vMinPoint.z)	vMinPoint.z = temp3f.z;
				if (temp3f.x > vMaxPoint.x)	vMaxPoint.x = temp3f.x;
				if (temp3f.y > vMaxPoint.y)	vMaxPoint.y = temp3f.y;
				if (temp3f.z > vMaxPoint.z)	vMaxPoint.z = temp3f.z;
			}
		}
		m_Center = 0.5  * (vMinPoint + vMaxPoint);
		m_Cube = (vMaxPoint - vMinPoint)*1.01;
		//if (m_Cube.x <EPS)	m_Cube.x = 1;
		//if (m_Cube.y <EPS)	m_Cube.y = 1;
		//if (m_Cube.z <EPS)	m_Cube.z = 1;
		if ((iNumTriangles < MaxEachOrderTrianglesNum) || (CurrentNode >= MaxOrderNum))
		{
			m_TID = new int[iNumTriangles];
			memcpy(m_TID, TrianglesID, sizeof(int)*iNumTriangles);
			m_ifEndNode = 1;
			return;
		}
		else
		{
			double maxlen = m_Cube.x;
			int axis = 0;
			if (m_Cube.y>maxlen)
			{
				maxlen = m_Cube.y;
				axis = 1;
			}
			if (m_Cube.z>maxlen)
			{
				maxlen = m_Cube.z;
				axis = 2;
			}

			bool* pBoolLeft = new bool[iNumTriangles];
			bool* pBoolRight = new bool[iNumTriangles];
			//将pBoolLeft和pBoolRight初始化为0
			memset(pBoolLeft, 0, sizeof(bool)*iNumTriangles);
			memset(pBoolRight, 0, sizeof(bool)*iNumTriangles);

			int CountLeft = 0, CountRight = 0;
			if (axis == 0)
			{
				for (i = 0; i<iNumTriangles; ++i)
				{
					if (Vertices[3 * TrianglesID[i]].x <= m_Center.x || Vertices[3 * TrianglesID[i] + 1].x <= m_Center.x || Vertices[3 * TrianglesID[i] + 2].x <= m_Center.x)
					{
						pBoolLeft[i] = true;
						CountLeft++;
					}
					if (Vertices[3 * TrianglesID[i]].x >= m_Center.x || Vertices[3 * TrianglesID[i] + 1].x >= m_Center.x || Vertices[3 * TrianglesID[i] + 2].x >= m_Center.x)
					{
						pBoolRight[i] = true;
						CountRight++;
					}
				}

				int *LeftID = new int[CountLeft];//左边区域三角形的编号
				int *RightID = new int[CountRight];//右边区域三角形的编号

				int Lefti = 0, Righti = 0;
				for (i = 0; i<iNumTriangles; ++i)
				{
					if (pBoolLeft[i])
					{
						LeftID[Lefti] = TrianglesID[i];
						Lefti++;
					}
					if (pBoolRight[i])
					{
						RightID[Righti] = TrianglesID[i];
						Righti++;
					}
				}
				m_Left = new Kdtree;
				m_Right = new Kdtree;
				CurrentNode++;
				m_Left->CreateNode(Vertices, LeftID, CountLeft);
				m_Right->CreateNode(Vertices, RightID, CountRight);
				CurrentNode--;
				delete[]  LeftID, RightID;
				LeftID = nullptr; RightID = nullptr;
			}

			if (axis == 1)
			{
				for (i = 0; i<iNumTriangles; ++i)
				{
					if (Vertices[3 * TrianglesID[i]].y <= m_Center.y || Vertices[3 * TrianglesID[i] + 1].y <= m_Center.y || Vertices[3 * TrianglesID[i] + 2].y <= m_Center.y)
					{
						pBoolLeft[i] = true;
						CountLeft++;
					}
					if (Vertices[3 * TrianglesID[i]].y >= m_Center.y || Vertices[3 * TrianglesID[i] + 1].y >= m_Center.y || Vertices[3 * TrianglesID[i] + 2].y >= m_Center.y)
					{
						pBoolRight[i] = true;
						CountRight++;
					}
				}

				int *LeftID = new int[CountLeft];
				int *RightID = new int[CountRight];

				int Lefti = 0, Righti = 0;
				for (i = 0; i<iNumTriangles; ++i)
				{
					if (pBoolLeft[i])
					{
						LeftID[Lefti] = TrianglesID[i];
						Lefti++;
					}
					if (pBoolRight[i])
					{
						RightID[Righti] = TrianglesID[i];
						Righti++;
					}
				}
				m_Left = new Kdtree;
				m_Right = new Kdtree;
				CurrentNode++;
				m_Left->CreateNode(Vertices, LeftID, CountLeft);
				m_Right->CreateNode(Vertices, RightID, CountRight);
				CurrentNode--;
				delete[]  LeftID, RightID;
				LeftID = nullptr; RightID = nullptr;
			}

			if (axis == 2)
			{
				for (i = 0; i<iNumTriangles; ++i)
				{
					if (Vertices[3 * TrianglesID[i]].z <= m_Center.z || Vertices[3 * TrianglesID[i] + 1].z <= m_Center.z || Vertices[3 * TrianglesID[i] + 2].z <= m_Center.z)
					{
						pBoolLeft[i] = true;
						CountLeft++;
					}
					if (Vertices[3 * TrianglesID[i]].z >= m_Center.z || Vertices[3 * TrianglesID[i] + 1].z >= m_Center.z || Vertices[3 * TrianglesID[i] + 2].z >= m_Center.z)
					{
						pBoolRight[i] = true;
						CountRight++;
					}
				}

				int *LeftID = new int[CountLeft];
				int *RightID = new int[CountRight];

				int Lefti = 0, Righti = 0;
				for (i = 0; i<iNumTriangles; ++i)
				{
					if (pBoolLeft[i])
					{
						LeftID[Lefti] = TrianglesID[i];
						Lefti++;
					}
					if (pBoolRight[i])
					{
						RightID[Righti] = TrianglesID[i];
						Righti++;
					}
				}
				m_Left = new Kdtree;
				m_Right = new Kdtree;
				CurrentNode++;
				m_Left->CreateNode(Vertices, LeftID, CountLeft);
				m_Right->CreateNode(Vertices, RightID, CountRight);
				CurrentNode--;
				delete[]  LeftID, RightID;
				LeftID = nullptr; RightID = nullptr;
			}
			delete[] pBoolLeft; pBoolLeft = nullptr;
			delete[] pBoolRight; pBoolRight = nullptr;
		}
	}

	void Release()
	{
		if (m_ifEndNode == 1)
		{
			delete[] m_TID;
			m_TID = nullptr;
		}
		if (m_Left)
		{
			m_Left->Release();
			delete[] m_Left;
			m_Left = nullptr;
		}
		if (m_Right)
		{
			m_Right->Release();
			delete[] m_Right;
			m_Right = nullptr;
		}
	}
};
#endif

// 将Octee中的内存变量的地址连续存储在一个数组变量中，便于遍历。
#ifndef _KDTREE_INT_
#define _KDTREE_INT_
class Kdtree_int
{ 
public:
	Kdtree *Node;
	int		Num = 0;//此节点下的子节点的个数(包括节点本身)
};
#endif

#ifndef _TARGET_
#define _TARGET_
class Target{
public:

	Complex	er;
	Complex	mr;	
	
	int TNum;//三角形的个数
	Point3f	*P;	//point
	Point3f	*N;	//法向
	Point3f	*C;	//中心
	double	*S;	// 面积
	double	*L;	//最长边长度
	Point3i	*E;	// three near triangles' ID，每一个三角形的相邻三个三角形的编号

	int *TrianglesID; //每一个三角形的编号
	Kdtree_int *TreePointer;
	Kdtree mytree;

	Target()
	{
		TNum = 0;
		P = 0;
		N = 0;
		C = 0;	
		S = 0;
		L = 0;	

		E = 0;
		TrianglesID = 0;
		TreePointer = new Kdtree_int();
	}

	void Create()
	{
		if (TNum>0)
		{		
			P = new Point3f[TNum*3];
			N = new Point3f[TNum];
			C = new Point3f[TNum];
			S = new double[TNum];
			L = new double[TNum];
			E = new Point3i[TNum];
		
			TrianglesID = new int[TNum];
			for (int i = 0; i < TNum; ++i)
			{
				TrianglesID[i] = i;
				
			}
			TreePointer = new Kdtree_int[int((pow(2.0,MaxOrderNum+1)-1))];
			
		}
	}

	int CreateDevTree(Kdtree *pNode, Kdtree_int *TreePointer)
	{
		int n = 1;
		int n0 = 0, n1 = 0;
		int M = NTEMP; NTEMP++;
		TreePointer[M].Node = pNode;
		if (pNode->m_Left != nullptr)
		{
			n0 = CreateDevTree(pNode->m_Left, TreePointer);
		}
		if (pNode->m_Right != nullptr)
		{
			n1 = CreateDevTree(pNode->m_Right, TreePointer);
		}
		n = n + n0 + n1;
		TreePointer[M].Num = n;
		return n;
	}

	void Release()
	{
		if (TNum>0)
		{
			delete[] P; P = nullptr;
			delete[] N; N = nullptr;
			delete[] N; N = nullptr;
			delete[] S; S = nullptr;
			delete[] L; L = nullptr;
			delete[] E; E = nullptr;
			delete[]TrianglesID; TrianglesID = nullptr;
			delete[] TreePointer; TreePointer = nullptr;
			mytree.Release();	
		}
	}
};
#endif

#ifndef _ROUGHSURFACE_
#define _ROUGHSURFACE_
class RoughSurface{
public:
	string name;
	PowerSpectrum Spectrum;
	double Seed;
	double Lx;
	double Ly;					// 粗糙面的大小 单位m
	int	   LNx;
	int    LNy;					// x方向分了LNx个点，y方向分了LNy个点 
	double High;  //粗糙面平移高度

	double delta;       //均方根高度
	double rlx;
	double rly;         //相关长度

	double WindSpeed;			// 风速
	double WindTheta;			// 风向
	double DeltT;				// 时域海面的时间间隔
	double Temp;				// 海水温度
	double Sault;				// 盐度
	double Wave_H;				// 浪高

	Complex er;
	Complex mr;					// 相对介电常数
};
#endif

#ifndef _RADAR_
#define _RADAR_
class Radar{
public:

	int ThNum, PhNum, RNum, tNum;
	double *Theta;
	double *Phi;
	double *R;
	Radar()
	{
		ThNum = 0; PhNum = 0; RNum = 0; tNum = 0;
		Theta = 0; Phi = 0; R = 0;
	}
	void Creat(double ThetaBegin, double ThetaEnd, double dTheta, double PhiBegin, double PhiEnd, double dPhi, double RBegin, double REnd, double dR)
	{
		ThNum = dTheta == 0 ? 1 : int((ThetaEnd - ThetaBegin) / dTheta) + 1;
		PhNum = dPhi == 0 ? 1 : int((PhiEnd - PhiBegin) / dPhi) + 1;
		RNum = dR == 0 ? 1 : int((REnd - RBegin) / dR) + 1;

		Theta = new double[ThNum];
		Phi = new double[PhNum];
		R = new double[RNum];

		int i;
		for (i = 0; i < ThNum; ++i)
			Theta[i] = ThetaBegin + i * dTheta;

		for (i = 0; i < PhNum; ++i)
			Phi[i] = PhiBegin + i * dPhi;

		for (i = 0; i < RNum; ++i)
			R[i] = RBegin + i * dR;
	}

	void Release()
	{
		delete[] Theta, Phi, R;
		Theta = nullptr; Phi = nullptr; R = nullptr;
	}
};
#endif

#ifndef _RAY_
#define _RAY_
class Ray{
public:
	Point3f Dir;		//方向
	Point3f orig;
};
#endif

#ifndef _RESULT_
#define _RESULT_
class Result{
public:

	int ThNum, PhNum, RNum ;

	int ***tNum;
	double ****t;
	double ****Pt;
	Point3c   ****Ei;
	Complex   ****Es;
	Complex   ****EsDF;
	Complex   ****EsPTD;
	Complex   ****Es1;
	Complex   ****Es2;
	Complex   ****Es12;
	Complex   ****Es21;
	Complex   ****Es212;

	void Creat(int thNum, int phNum, int rNum)
	{
		ThNum = thNum;
		PhNum = phNum;
		RNum = rNum;
		tNum = new   int**[ThNum];
		t = new   double***[ThNum];
		Pt = new   double***[ThNum];
		Ei = new   Point3c***[ThNum];
		Es = new   Complex***[ThNum];
		EsDF = new   Complex***[ThNum];
		EsPTD = new   Complex***[ThNum];
		Es1 = new   Complex***[ThNum];
		Es2 = new   Complex***[ThNum];
		Es12 = new   Complex***[ThNum];
		Es21 = new   Complex***[ThNum];
		Es212 = new   Complex***[ThNum];
		for (int j = 0; j<ThNum; ++j)
		{
			tNum[j] = new   int*[PhNum];
			t[j] = new   double**[PhNum];
			Pt[j] = new   double**[PhNum];
			Ei[j] = new Point3c**[PhNum];
			Es[j] = new Complex**[PhNum];
			EsDF[j] = new Complex**[PhNum];
			EsPTD[j] = new Complex**[PhNum];
			Es1[j] = new Complex**[PhNum];
			Es2[j] = new Complex**[PhNum];
			Es12[j] = new Complex**[PhNum];
			Es21[j] = new Complex**[PhNum];
			Es212[j] = new Complex**[PhNum];
			for (int k = 0; k < PhNum; ++k)
			{
				tNum[j][k] = new   int[RNum];
				t[j][k] = new   double*[RNum];
				Pt[j][k] = new   double*[RNum];
				Ei[j][k] = new Point3c*[RNum];
				Es[j][k] = new Complex*[RNum];
				EsDF[j][k] = new Complex*[RNum];
				EsPTD[j][k] = new Complex*[RNum];
				Es1[j][k] = new Complex*[RNum];
				Es2[j][k] = new Complex*[RNum];
				Es12[j][k] = new Complex*[RNum];
				Es21[j][k] = new Complex*[RNum];
				Es212[j][k] = new Complex*[RNum];
			}
		}
	}

	void CreatTime(int i, int j, int k, double tBegin, double tEnd, double dt)
	{
		int tnum;
		tnum = dt == 0 ? 1 : int((tEnd - tBegin) / dt) + 1;
		
		tNum[i][j][k] = tnum;
		t[i][j][k] = new double[tnum];
		Pt[i][j][k] = new double[tnum];
		Ei[i][j][k] = new   Point3c[tnum];
		Es[i][j][k] = new   Complex[tnum];
		EsDF[i][j][k] = new   Complex[tnum];
		EsPTD[i][j][k] = new   Complex[tnum];
		Es1[i][j][k] = new   Complex[tnum];
		Es2[i][j][k] = new   Complex[tnum];
		Es12[i][j][k] = new   Complex[tnum];
		Es21[i][j][k] = new   Complex[tnum];
		Es212[i][j][k] = new   Complex[tnum];

		for (int m = 0; m < tnum; ++m)
			t[i][j][k][m] = tBegin + m* dt;
	}

	void Release()
	{
		int i, j, k;

		for (i = 0; i<ThNum; ++i)
		{
			for (j = 0; j<PhNum; ++j)
			{
				for (k = 0; k<RNum; ++k)
				{
					delete[] t[i][j][k], Pt[i][j][k], Ei[i][j][k], Es[i][j][k], EsDF[i][j][k], EsPTD[i][j][k], Es1[i][j][k], Es2[i][j][k], Es12[i][j][k], Es21[i][j][k], Es212[i][j][k];
				}
				delete[] tNum[i][j], t[i][j], Pt[i][j], Ei[i][j], Es[i][j], EsDF[i][j], EsPTD[i][j], Es1[i][j], Es2[i][j], Es12[i][j], Es21[i][j], Es212[i][j];
			}
			delete[]  tNum[i], t[i], Pt[i], Ei[i], Es[i], EsDF[i], EsPTD[i], Es1[i], Es2[i], Es12[i], Es21[i], Es212[i];
		}
		delete[] tNum, t, Pt, Ei, Es, EsDF, EsPTD, Es1, Es2, Es12, Es21, Es212;
		Ei = nullptr; Es = nullptr; EsDF = nullptr; EsPTD = nullptr; Es1 = nullptr; Es2 = nullptr; Es12 = nullptr; Es21 = nullptr; Es212 = nullptr;
	}

};
#endif

#endif