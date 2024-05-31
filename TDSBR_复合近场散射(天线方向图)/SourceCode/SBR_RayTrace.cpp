#include "SBR_RayTrace.h"
bool RayTriangle(Ray &ray, Target &Tar, int ID)
{
	Point3f P1 = Tar.P[3 * ID];
	Point3f P2 = Tar.P[3 * ID + 1];
	Point3f P3 = Tar.P[3 * ID + 2];

	Point3f center = ray.orig;

	double u, v;
	Point3f pvec = ray.Dir % (P3 - P1);

	double det = (P2 - P1) * pvec;
	Point3f tvec;
	if (det > 0)
		tvec = center - P1;
	else
	{
		tvec = P1 - center;
		det = -det;
	}
	if (det < -EPS)
		return false;

	u = tvec * pvec;
	if (u < -EPS || u > det)
		return false;

	pvec = tvec % (P2 - P1);
	v = ray.Dir * pvec;
	if (v < -EPS || u + v > det)
		return false;

	if ((P3 - P1) * pvec <= 0)
		return false;

	return true;
}

bool RayBox(Point3f orig, Point3f dir, Point3f center, Point3f cube, double * len)
{
	Point3f vmin = center - 0.5  * cube;
	Point3f vmax = center + 0.5  * cube;

	if (abs(dir.x)<EPS)
	{
		dir.x = EPS;
	}
	if (abs(dir.y)<EPS)
	{
		dir.y = EPS;
	}
	if (abs(dir.z)<EPS)
	{
		dir.z = EPS;
	}

	double tmin, tmax, tymin, tymax, tzmin, tzmax;
	if (dir.x >= 0)
	{
		tmin = (vmin.x - orig.x) / dir.x;
		tmax = (vmax.x - orig.x) / dir.x;
	}
	else
	{
		tmin = (vmax.x - orig.x) / dir.x;
		tmax = (vmin.x - orig.x) / dir.x;
	}

	if (dir.y >= 0)
	{
		tymin = (vmin.y - orig.y) / dir.y;
		tymax = (vmax.y - orig.y) / dir.y;
	}
	else
	{
		tymin = (vmax.y - orig.y) / dir.y;
		tymax = (vmin.y - orig.y) / dir.y;
	}

	if (tmin>tymax || tymin>tmax)
	{
		tmin = -1;
		*len = tmin;
		return false;
	}

	if (tymin > tmin)
	{
		tmin = tymin;
	}

	if (tymax < tmax)
	{
		tmax = tymax;
	}

	if (dir.z >= 0)
	{
		tzmin = (vmin.z - orig.z) / dir.z;
		tzmax = (vmax.z - orig.z) / dir.z;
	}
	else
	{
		tzmin = (vmax.z - orig.z) / dir.z;
		tzmax = (vmin.z - orig.z) / dir.z;
	}

	if (tmin>tzmax || tzmin>tmax)
	{
		tmin = -1;
		*len = tmin;
		return false;
	}

	if (tzmin > tmin)
	{
		tmin = tzmin;
	}

	if (tzmax < tmax)
	{
		tmax = tzmax;
	}

	*len = tmin;
	return true;
}


int RayTree1(Target &Tar, Ray &ray)
{
	double len = 0;
	int i = 0, j;
	Kdtree *inode;
	Point3f center = ray.orig;
	int TID;
	while (i < Tar.TreePointer[0].Num)
	{
		inode = Tar.TreePointer[i].Node;
		if (RayBox(center, ray.Dir, inode->m_Center, inode->m_Cube, &len))
		{
			if (inode->m_ifEndNode == 1)
			{
				for (j = 0; j < inode->m_NumT; ++j)
				{
					TID = inode->m_TID[j];
					if ((center - Tar.C[TID]).norm() > EPS && RayTriangle(ray, Tar, TID))
					{
						return TID;
					}
				}
			}
			++i;
			continue;
		}
		else
		{
			i = i + Tar.TreePointer[i].Num;
			continue;
		}
	}
	return -1;
}


int RayTree2(Target &Tar, Ray &ray)
{
	double len = 0;
	int i = 0, j;
	double aaa = EPB, bbb;
	int re = -1;
	Kdtree *inode;
	Point3f center = ray.orig;
	int TID;
	while (i < Tar.TreePointer[0].Num)
	{
		inode = Tar.TreePointer[i].Node;
		if (RayBox(center, ray.Dir, inode->m_Center, inode->m_Cube, &len))
		{
			if (len > aaa)
			{
				i = i + Tar.TreePointer[i].Num;
				continue;
			}
			if (inode->m_ifEndNode == 1)
			{
				for (j = 0; j<inode->m_NumT; ++j)
				{
					TID = inode->m_TID[j];
					if ((center - Tar.C[TID]).norm() > EPS && RayTriangle(ray, Tar, TID))
					{
						bbb = (center - Tar.C[TID]).norm();
						if (bbb<aaa)
						{
							aaa = bbb;
							re = TID;
						}
					}
				}
			}
			++i;
			continue;
		}
		else
		{
			i = i + Tar.TreePointer[i].Num;
			continue;
		}
	}
	return re;
}


void RefN_CPU1(Target &Tar, Target &Tar1, Point3f AntennaPosition, vector<int> &Index)
{
	int i = 0;
	Ray ray;

	for (i = 0; i<Tar.TNum; ++i)
	{
		ray.orig = Tar.C[i];
		ray.Dir = (AntennaPosition - ray.orig).normlize();
		if (ray.Dir * Tar.N[i] >0 && RayTree1(Tar, ray) == -1 && RayTree1(Tar1, ray) == -1)
		{
			Index.push_back(i);
		}
	}
	return;
}

void RefN_CPU12(Target &Tar1, Target &Tar2, bool bl, Point3f AntennaPosition, vector<int> &Index1, vector<int*> &Index2)
{
	int i = 0, j, p, m, m1, num, n;
	Ray ray;
	Point3f AntennaPositionMirror1;
	set<int> wedgeID;
	set<int>::iterator it;
	int *TreeE = new int[Tar2.TNum];
	for (n = 0; n<(int)Index1.size(); ++n)
	{
		i = Index1[n];
		AntennaPositionMirror1 = AntennaPosition - 2.0f * (AntennaPosition - Tar1.C[i])* Tar1.N[i] * Tar1.N[i];
		ray.orig = Tar1.C[i];
		ray.Dir = (ray.orig - AntennaPositionMirror1).normlize();
		p = RayTree2(Tar2, ray);
		if (!bl && RayTree1(Tar1, ray) != -1)
			p = -1;
		m1 = 0, num = 0;
		if (p != -1)
		{
			TreeE[num] = p; num++;
			if (Tar2.E[p].x != -1) { TreeE[num] = Tar2.E[p].x; num++; }
			if (Tar2.E[p].y != -1) { TreeE[num] = Tar2.E[p].y; num++; }
			if (Tar2.E[p].z != -1) { TreeE[num] = Tar2.E[p].z; num++; }
		}
		while (m1<num)
		{
			m = TreeE[m1];
			ray.orig = Tar2.C[m];
			ray.Dir = (AntennaPositionMirror1 - ray.orig).normlize();
			if (Tar2.N[m] * ray.Dir > 0 &&
				(Tar1.C[i] - Tar2.C[m] - (Tar1.C[i] - Tar2.C[m]) * ray.Dir * ray.Dir).norm() < Tar1.L[i]
				&& (bl || RayTree1(Tar2, ray) == -1)
				&& RayTree2(Tar1, ray) == i)
			{
				Index2.push_back(new int[2]{ i, m });

				j = 0;
				if (Tar2.E[m].x != -1)
				{
					for (j = 0; j<num; ++j)
					{
						if (TreeE[j] == Tar2.E[m].x)
							break;
					}
					if (j == num)
					{
						TreeE[num] = Tar2.E[m].x;
						num++;
					}
				}

				j = 0;
				if (Tar2.E[m].y != -1)
				{
					for (j = 0; j<num; ++j)
					{
						if (TreeE[j] == Tar2.E[m].y)
							break;
					}
					if (j == num)
					{
						TreeE[num] = Tar2.E[m].y;
						num++;
					}
				}

				j = 0;
				if (Tar2.E[m].z != -1)
				{
					for (j = 0; j<num; ++j)
					{
						if (TreeE[j] == Tar2.E[m].z)
							break;
					}
					if (j == num)
					{
						TreeE[num] = Tar2.E[m].z;
						num++;
					}
				}
			}
			m1++;
		}
	}

	delete[] TreeE;
	return;
}

//三次弹跳，计算雷达经 Tar1 到 Tar1再到Tar2 的射线. bl：Tar1和Tar2是否是同一目标，true：是。
void RefN_CPU112(Target &Tar1, Target &Tar2, bool bl, Point3f AntennaPosition, vector<int*> &Index11, vector<int*> &Index112)
{
	int ID1, ID2, ID3;
	int  j, p, m1, num, n;
	Ray ray;
	Point3f AntennaPositionMirror1, AntennaPositionMirror2;
	set<int> wedgeID;
	set<int>::iterator it;
	int *TreeE = new int[Tar2.TNum];
	for (n = 0; n<(int)Index11.size(); ++n)
	{
		ID1 = Index11[n][0];
		ID2 = Index11[n][1];
		AntennaPositionMirror1 = AntennaPosition - 2.0f * (AntennaPosition - Tar1.C[ID1])* Tar1.N[ID1] * Tar1.N[ID1];
		AntennaPositionMirror2 = AntennaPositionMirror1 - 2.0f * (AntennaPositionMirror1 - Tar1.C[ID2])* Tar1.N[ID2] * Tar1.N[ID2];
		ray.orig = Tar1.C[ID2];
		ray.Dir = (ray.orig - AntennaPositionMirror2).normlize();
		p = RayTree2(Tar2, ray);
		if (!bl && RayTree1(Tar1, ray) != -1)
			p = -1;
		m1 = 0, num = 0;
		if (p != -1)//找相交三角形的邻边三角形
		{
			TreeE[num] = p; num++;
			if (Tar2.E[p].x != -1) { TreeE[num] = Tar2.E[p].x; num++; }
			if (Tar2.E[p].y != -1) { TreeE[num] = Tar2.E[p].y; num++; }
			if (Tar2.E[p].z != -1) { TreeE[num] = Tar2.E[p].z; num++; }
		}
		while (m1<num)
		{
			ID3 = TreeE[m1];
			ray.orig = Tar2.C[ID3];
			ray.Dir = (AntennaPositionMirror2 - ray.orig).normlize();
			if (Tar2.N[ID3] * ray.Dir > 0 &&
				(Tar1.C[ID2] - Tar2.C[ID3] - (Tar1.C[ID2] - Tar2.C[ID3]) * ray.Dir * ray.Dir).norm() <Tar1.L[ID2]
				&& (bl || RayTree1(Tar2, ray) == -1)
				&& RayTree2(Tar1, ray) == ID2)
			{
				Index112.push_back(new int[3]{ ID1, ID2, ID3 });

				j = 0;
				if (Tar2.E[ID3].x != -1)//找被照射三角形的邻边三角形
				{
					for (j = 0; j<num; ++j)
					{
						if (TreeE[j] == Tar2.E[ID3].x)
							break;
					}
					if (j == num)
					{
						TreeE[num] = Tar2.E[ID3].x;
						num++;
					}
				}

				j = 0;
				if (Tar2.E[ID3].y != -1)
				{
					for (j = 0; j<num; ++j)
					{
						if (TreeE[j] == Tar2.E[ID3].y)
							break;
					}
					if (j == num)
					{
						TreeE[num] = Tar2.E[ID3].y;
						num++;
					}
				}

				j = 0;
				if (Tar2.E[ID3].z != -1)
				{
					for (j = 0; j<num; ++j)
					{
						if (TreeE[j] == Tar2.E[ID3].z)
							break;
					}
					if (j == num)
					{
						TreeE[num] = Tar2.E[ID3].z;
						num++;
					}
				}

			}
			m1++;

		}

	}

	delete[] TreeE;
	return;
}

//四次弹跳，计算雷达经 Tar1 到 Tar1到 Tar1再到Tar2 的射线. bl：Tar1和Tar2是否是同一目标，true：是。
void RefN_CPU1112(Target &Tar1, Target &Tar2, bool bl, Point3f AntennaPosition, vector<int*> &Index111, vector<int*> &Index1112)
{
	int ID1, ID2, ID3, ID4;
	int  j, p, m1, num, n;
	Ray ray;
	Point3f AntennaPositionMirror1, AntennaPositionMirror2, AntennaPositionMirror3;
	set<int> wedgeID;
	set<int>::iterator it;
	int *TreeE = new int[Tar2.TNum];
	for (n = 0; n<(int)Index111.size(); ++n)
	{
		ID1 = Index111[n][0];
		ID2 = Index111[n][1];
		ID3 = Index111[n][2];
		AntennaPositionMirror1 = AntennaPosition - 2.0f * (AntennaPosition - Tar1.C[ID1])* Tar1.N[ID1] * Tar1.N[ID1];
		AntennaPositionMirror2 = AntennaPositionMirror1 - 2.0f * (AntennaPositionMirror1 - Tar1.C[ID2])* Tar1.N[ID2] * Tar1.N[ID2];
		AntennaPositionMirror3 = AntennaPositionMirror2 - 2.0f * (AntennaPositionMirror2 - Tar1.C[ID3])* Tar1.N[ID3] * Tar1.N[ID3];
		ray.orig = Tar1.C[ID3];
		ray.Dir = (ray.orig - AntennaPositionMirror3).normlize();
		p = RayTree2(Tar2, ray);
		if (!bl && RayTree1(Tar1, ray) != -1)
			p = -1;
		m1 = 0, num = 0;
		if (p != -1)//找相交三角形的邻边三角形
		{
			TreeE[num] = p; num++;
			if (Tar2.E[p].x != -1) { TreeE[num] = Tar2.E[p].x; num++; }
			if (Tar2.E[p].y != -1) { TreeE[num] = Tar2.E[p].y; num++; }
			if (Tar2.E[p].z != -1) { TreeE[num] = Tar2.E[p].z; num++; }
		}
		while (m1<num)
		{
			ID4 = TreeE[m1];
			ray.orig = Tar2.C[ID4];
			ray.Dir = (AntennaPositionMirror2 - ray.orig).normlize();
			if (Tar2.N[ID4] * ray.Dir > 0 &&
				(Tar1.C[ID3] - Tar2.C[ID4] - (Tar1.C[ID3] - Tar2.C[ID4]) * ray.Dir * ray.Dir).norm() <Tar1.L[ID3]
				&& (bl || RayTree1(Tar2, ray) == -1)
				&& RayTree2(Tar1, ray) == ID3)
			{
				Index1112.push_back(new int[4]{ ID1, ID2, ID3, ID4 });

				j = 0;
				if (Tar2.E[ID4].x != -1)//找被照射三角形的邻边三角形
				{
					for (j = 0; j<num; ++j)
					{
						if (TreeE[j] == Tar2.E[ID4].x)
							break;
					}
					if (j == num)
					{
						TreeE[num] = Tar2.E[ID4].x;
						num++;
					}
				}

				j = 0;
				if (Tar2.E[ID4].y != -1)
				{
					for (j = 0; j<num; ++j)
					{
						if (TreeE[j] == Tar2.E[ID4].y)
							break;
					}
					if (j == num)
					{
						TreeE[num] = Tar2.E[ID4].y;
						num++;
					}
				}

				j = 0;
				if (Tar2.E[ID4].z != -1)
				{
					for (j = 0; j<num; ++j)
					{
						if (TreeE[j] == Tar2.E[ID4].z)
							break;
					}
					if (j == num)
					{
						TreeE[num] = Tar2.E[ID4].z;
						num++;
					}
				}

			}
			m1++;

		}
	}

	delete[] TreeE;
	return;
}

//三次弹跳，计算雷达经 Tar1 到 Tar2再到Tar1 的射线 , Tar1 和 Tar2不能是同一个目标
void RefN_CPU121(Target &Tar1, Target &Tar2, Point3f AntennaPosition, vector<int*> &Index12, vector<int*> &Index121)
{
	int ID1, ID2, ID3;
	int  j, p, m1, num, n;
	Ray ray;
	Point3f AntennaPositionMirror1, AntennaPositionMirror2;
	set<int> wedgeID;
	set<int>::iterator it;
	int *TreeE = new int[Tar2.TNum];
	for (n = 0; n<(int)Index12.size(); ++n)
	{
		ID1 = Index12[n][0];
		ID2 = Index12[n][1];
		AntennaPositionMirror1 = AntennaPosition - 2.0f * (AntennaPosition - Tar1.C[ID1])* Tar1.N[ID1] * Tar1.N[ID1];
		AntennaPositionMirror2 = AntennaPositionMirror1 - 2.0f * (AntennaPositionMirror1 - Tar2.C[ID2])* Tar2.N[ID2] * Tar2.N[ID2];
		ray.orig = Tar2.C[ID2];
		ray.Dir = (ray.orig - AntennaPositionMirror2).normlize();
		p = RayTree2(Tar1, ray);
		if (RayTree1(Tar2, ray) != -1)
			p = -1;
		m1 = 0, num = 0;
		if (p != -1)//找相交三角形的邻边三角形
		{
			TreeE[num] = p; num++;
			if (Tar1.E[p].x != -1) { TreeE[num] = Tar1.E[p].x; num++; }
			if (Tar1.E[p].y != -1) { TreeE[num] = Tar1.E[p].y; num++; }
			if (Tar1.E[p].z != -1) { TreeE[num] = Tar1.E[p].z; num++; }
		}
		while (m1<num)
		{
			ID3 = TreeE[m1];
			ray.orig = Tar1.C[ID3];
			ray.Dir = (AntennaPositionMirror2 - ray.orig).normlize();
			if (Tar1.N[ID3] * ray.Dir > 0 &&
				(Tar2.C[ID2] - Tar1.C[ID3] - (Tar2.C[ID2] - Tar1.C[ID3]) * ray.Dir * ray.Dir).norm() <Tar2.L[ID2]
				&& RayTree1(Tar1, ray) == -1
				&& RayTree2(Tar2, ray) == ID2)
			{
				Index121.push_back(new int[3]{ ID1, ID2, ID3 });

				j = 0;
				if (Tar1.E[ID3].x != -1)//找被照射三角形的邻边三角形
				{
					for (j = 0; j<num; ++j)
					{
						if (TreeE[j] == Tar1.E[ID3].x)
							break;
					}
					if (j == num)
					{
						TreeE[num] = Tar1.E[ID3].x;
						num++;
					}
				}

				j = 0;
				if (Tar1.E[ID3].y != -1)
				{
					for (j = 0; j<num; ++j)
					{
						if (TreeE[j] == Tar1.E[ID3].y)
							break;
					}
					if (j == num)
					{
						TreeE[num] = Tar1.E[ID3].y;
						num++;
					}
				}

				j = 0;
				if (Tar1.E[ID3].z != -1)
				{
					for (j = 0; j<num; ++j)
					{
						if (TreeE[j] == Tar1.E[ID3].z)
							break;
					}
					if (j == num)
					{
						TreeE[num] = Tar1.E[ID3].z;
						num++;
					}
				}

			}
			m1++;

		}
	}

	delete[] TreeE;
	return;
}
