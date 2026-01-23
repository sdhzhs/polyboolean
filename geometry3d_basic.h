#pragma once

#include <math.h>
#include <stdbool.h>

#define TOLERANCE 1e-8 //geometric tolerence
typedef double real;

//comparison function of float pointing number based on tolerence
static inline int FltGT(real a, real b)
{
  return a-b > TOLERANCE;
}

static inline int FltGE(real a, real b)
{
  return a-b >= -TOLERANCE;
}

static inline int FltLT(real a, real b)
{
  return a-b < -TOLERANCE;
}

static inline int FltLE(real a, real b)
{
  return a-b <= TOLERANCE;
}

static inline int FltEQ(real a, real b)
{
  return fabs(a-b) <= TOLERANCE;
}

real DotProduct(real* v1, real* v2, int ndims);
void CrossProduct(real v1[3], real v2[3], real v3[3]);
real PntProject(real c[3], real n[3], real p[3], real pp[3]);
int PntShadow(real c[3], real n[3], real t[3], real p[3], real ps[3]);
int Isinner(real p[][3], real p0[3], int nps);
int InterSect(real* p0, real* p1, real* p2, real* p3, real* pi);
real PolyArea(real p[][3], int nps);
void PolyProp(real p[][3], real c[3], real n[3], int nps);
void InsertSort(int* in, int* id, int nps);
void PolyBound(real poly[][3], real boundbox[2][3], int nps);
bool IsAxisIntersect(real axis1[2], real axis2[2]);
bool IsBBoxIntersect(real boundbox1[2][3], real boundbox2[2][3]);