#include <stdbool.h>

typedef double real;

real DotProduct(real* v1, real* v2, int ndims);
void CrossProduct(real v1[3], real v2[3], real v3[3]);
real PntProject(real c[3], real n[3], real p[3], real pp[3]);
int PntShadow(real c[3], real n[3], real t[3], real p[3], real ps[3]);
int Isinner(real p[][3], real p0[3], int nps);
int InterSect(real* p0, real* p1, real* p2, real* p3, real* pi);
real PolyArea(real p[][3], int nps);
void PolyProp(real p[][3], real c[3], real n[3], int nps);
void InsertSort(int* in, int* id, int nps);
void PolySlice(int slicebegin, int slicend, int* curid, int nps, real poly[][3], real polyslice[][3], bool startin, bool outloop);
bool OnPolyNodes(real p[][3], real p0[3], int* nodeid, int nps);
void PolyIntersect(real poly1[][3], real poly2[][3], real polyi[][3], int nps1, int nps2, int* npsi);
void PolyMerge(real poly1[][3], real poly2[][3], real polym[][3], int nps1, int nps2, int* npsm);
bool EraseSamePointsInPoly(real poly[][3], real polyn[][3], int nps, int* npsn);
