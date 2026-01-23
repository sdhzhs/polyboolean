#include <float.h>
#include "geometry3d_basic.h"

//basic geometric calculation
real DotProduct(real* v1, real* v2, int ndims) //dot product of two n-dimentional vectors
{
  int j;
  real dp;

  dp = 0.0;
  for(j=0; j<ndims; j++)
  {
    dp += v1[j]*v2[j];
  }

  return dp;
}

void CrossProduct(real v1[3], real v2[3], real v3[3]) //cross product of two 3-dimentional vectors
{
  v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
  v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
  v3[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

real PntProject(real c[3], real n[3], real p[3], real pp[3])  //normal projection of one point on one plane
{
  int j;
  real tp;

  tp = DotProduct(n,c,3) - DotProduct(n,p,3);
  for(j=0; j<3; j++)
  {
    pp[j] = p[j]+n[j]*tp;
  }

  return fabs(tp);
}

int PntShadow(real c[3], real n[3], real t[3], real p[3], real ps[3]) //ray projection of one point on one plane
{
  int j;
  real tp, tol;
  real v[3];

  tol = TOLERANCE*TOLERANCE;

  tp = DotProduct(n,t,3);

  if(fabs(tp)<tol)
  {
    for(j=0; j<3; j++)
    {
      ps[j] = 0;
    }
    return 0;
  }
  else
  {
    for(j=0; j<3; j++)
    {
      v[j] = c[j]-p[j];
    }
    tp = DotProduct(n,v,3)/tp;
    for(j=0; j<3; j++)
    {
      ps[j] = p[j]+t[j]*tp;
    }
    return 1;
  }
}

int Isinner(real p[][3], real p0[3], int nps) //predicate a point is inside a 3D convex polygon, including point on vertices, excluding point on edges
{
  int i, j, inext, id;
  real ip, tol, epsi, epsi0;
  real v0[3], v1[3], v2[3], v3[3];

  tol = TOLERANCE;

  for(j=0; j<3; j++)
    v0[j] = 0.0;

  epsi0 = 0.0;

  for(i=0; i<nps-1; i++)
  {
    for(j=0; j<3; j++)
    {
      v1[j] = p[i][j] - p0[j];
      v2[j] = p[i+1][j] - p0[j];
    }
    CrossProduct(v1, v2, v0);
    epsi0 = sqrt(DotProduct(v1,v1,3)*DotProduct(v2,v2,3));
    if(sqrt(DotProduct(v0,v0,3))>tol*epsi0) break;
  }

  id = 1;
  for(i=0; i<nps; i++)
  {
    inext = (i+1)%nps;
    for(j=0; j<3; j++)
    {
      v1[j] = p[i][j] - p0[j];
      v2[j] = p[inext][j] - p0[j];
    }
    CrossProduct(v1, v2, v3);
    epsi = tol*epsi0*tol*sqrt(DotProduct(v1,v1,3)*DotProduct(v2,v2,3));
    ip = DotProduct(v0,v3,3);
    if(ip<epsi)
    {
      id = 0;
      break;
    }
  }
  return id;
}

int InterSect(real* p0, real* p1, real* p2, real* p3, real* pi) //predicate intersecting status of two segments, excluding collinear and common point intersections
{
  int j, ndims;
  real ip, tp, sp, tol, epsi1, epsi2;
  real v0[3], v1[3], v2[3], v3[3], v4[3];

  ndims = 3;
  tol = TOLERANCE;

  for(j=0; j<ndims; j++)
  {
    v0[j] = p0[j] - p2[j];
    v1[j] = p1[j] - p2[j];
  }
  CrossProduct(v0, v1, v2);
  epsi1 = sqrt(DotProduct(v0,v0,ndims)*DotProduct(v1,v1,ndims));
  for(j=0; j<ndims; j++)
  {
    v0[j] = p0[j] - p3[j];
    v1[j] = p1[j] - p3[j];
  }
  CrossProduct(v0, v1, v3);
  epsi1 *= tol*sqrt(DotProduct(v0,v0,ndims)*DotProduct(v1,v1,ndims));
  ip = DotProduct(v2,v3,ndims);
  for(j=0; j<ndims; j++)
  {
    v0[j] = p2[j] - p0[j];
    v1[j] = p3[j] - p0[j];
  }
  CrossProduct(v0, v1, v2);
  epsi2 = sqrt(DotProduct(v0,v0,ndims)*DotProduct(v1,v1,ndims));
  for(j=0; j<ndims; j++)
  {
    v0[j] = p2[j] - p1[j];
    v1[j] = p3[j] - p1[j];
  }
  CrossProduct(v0, v1, v3);
  epsi2 *= tol*sqrt(DotProduct(v0,v0,ndims)*DotProduct(v1,v1,ndims));
  tp = DotProduct(v2,v3,ndims);
  if(ip < epsi1 && tp < epsi2 && (fabs(ip) > epsi1 || fabs(tp) > epsi2))
  //if(ip <= epsi1 && tp <= epsi2)
  {
    for(j=0; j<ndims; j++)
    {
      v0[j] = p1[j] - p0[j];
      v1[j] = p3[j] - p2[j];
      v2[j] = p2[j] - p0[j];
    }
    CrossProduct(v0, v1, v3);
    CrossProduct(v1, v2, v4);
    sp = sqrt(DotProduct(v4,v4,ndims)/DotProduct(v3,v3,ndims));
    for(j=0; j<ndims; j++)
    {
      pi[j] = p0[j] + v0[j]*sp;
    }
    return 1;
  }
  else
  {
    for(j=0; j<ndims; j++)
    {
      pi[j] = 0;
    }
    return 0;
  }
}

real PolyArea(real p[][3], int nps) //calculate area of a 3D polygon
{
  int i, j;
  real dA, A;
  real v0[3], v1[3], v2[3];

  A = 0;

  for(i=1; i<nps-1; i++)
  {
    for(j=0; j<3; j++)
    {
      v0[j] = p[i][j] - p[0][j];
      v1[j] = p[i+1][j] - p[0][j];
    }
    CrossProduct(v0, v1, v2);
    dA = 0.5 * sqrt(DotProduct(v2,v2,3));
    A += dA;
  }

  return A;
}

void PolyProp(real p[][3], real c[3], real n[3], int nps) //calculate normal vector and geometric center of a 3D polygon
{
  int i, j;
  real dn;
  real v0[3], v1[3], v2[3];

  for(j=0; j<3; j++)
  {
    c[j] = 0;
    n[j] = 0;
  }
  for(i=0; i<nps; i++)
  {
    for(j=0; j<3; j++)
    {
      c[j] += p[i][j];
    }
  }
  for(i=1; i<nps-1; i++)
  {
    for(j=0; j<3; j++)
    {
      v0[j] = p[i][j] - p[0][j];
      v1[j] = p[i+1][j] - p[0][j];
    }
    CrossProduct(v0, v1, v2);
    for(j=0; j<3; j++)
    {
      n[j] += v2[j];
    }
  }

  if(nps == 0) return;
  for(j=0; j<3; j++)
  {
    c[j] /= nps;
  }

  dn = sqrt(DotProduct(n,n,3));
  if(dn == 0.0) return;
  for(j=0; j<3; j++)
  {
    n[j] /= dn;
  }
}

void InsertSort(int* in, int* id, int nps)  //insert sort
{
  int i, j, key;

  for(j=1; j<nps; j++)
  {
    key = id[j];
    i = j-1;
    while(i>=0 && id[i]>key)
    {
      id[i+1] = id[i];
      in[i+1] = in[i];
      i--;
    }
    id[i+1] = key;
    in[i+1] = j;
  }
}

void PolyBound(real poly[][3], real boundbox[2][3], int nps) //get a 3D boundary box of a 3D polygon
{
  for(int j=0; j<3; j++)
  {
    boundbox[0][j] = DBL_MAX;
    boundbox[1][j] = -DBL_MAX;
  }

  for(int i=0; i<nps; i++)
  {
    for(int j=0; j<3; j++)
    {
      if(poly[i][j]<boundbox[0][j]) boundbox[0][j] = poly[i][j];
      if(poly[i][j]>boundbox[1][j]) boundbox[1][j] = poly[i][j];
    }
  }
}

bool IsAxisIntersect(real axis1[2], real axis2[2]) //check whether two ranges on 1D axis are intersected
{
  if(axis1[0] > axis2[1] || axis1[1] < axis2[0])
  {
    return false;
  }
  
  return true;
}

bool IsBBoxIntersect(real boundbox1[2][3], real boundbox2[2][3]) //check whether two 3D boundary boxes are intersected
{
  bool IsInsec = true;
  real axis1[2], axis2[2];

  for(int j=0; j<3; j++)
  {
    axis1[0] = boundbox1[0][j];
    axis1[1] = boundbox1[1][j];
    axis2[0] = boundbox2[0][j];
    axis2[1] = boundbox2[1][j];

    IsInsec = IsInsec && IsAxisIntersect(axis1, axis2);
  }

  return IsInsec;
}

// -------------------------------------------------------------