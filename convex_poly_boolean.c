#include <stdio.h>
#include <math.h>
#include <float.h>
#include "convex_poly_boolean.h"

#define ND_ND 3
#define MAX_INTERSECTIONS 4
#define TOLERANCE 1e-6

//comparison function of float pointing number based on tolerence
inline int FltGT(real a, real b)
{
  return a-b > TOLERANCE;
}

inline int FltGE(real a, real b)
{
  return a-b >= -TOLERANCE;
}

inline int FltLT(real a, real b)
{
  return a-b < -TOLERANCE;
}

inline int FltLE(real a, real b)
{
  return a-b <= TOLERANCE;
}

inline int FltEQ(real a, real b)
{
  return fabs(a-b) <= TOLERANCE;
}

//basic geometry calculation
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
  real tp;
  real v[3];

  tp = DotProduct(n,t,3);

  if(FltEQ(tp, 0.0))
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

  tol = TOLERANCE; //geometric tolerence

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
  tol = TOLERANCE; //geometric tolerence

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
  for(j=0; j<3; j++)
  {
    c[j] /= nps;
    n[j] /= nps-2;
  }
  dn = sqrt(DotProduct(n,n,3));
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

void PolySlice(int slicebegin, int slicend, int* curid, int nps, real poly[][3], real polyslice[][3], bool startin, bool outloop) //slice of polygon using begin and end indices, controlled by 2 flags
{
  int i, j, startid, endid;

  if(startin)
    startid = slicebegin;
  else
    startid = slicebegin+1;

  endid = slicend;

  if(slicend > slicebegin)
  {
    for(i=startid; i<=endid; i++)
    {
      for(j=0; j<3; j++)
      {
        polyslice[*curid][j] = poly[i][j];
      }
      (*curid)++;
    }
  }
  else if(slicend < slicebegin || (slicend == slicebegin && outloop))
  {
    for(i=startid; i<nps; i++)
    {
      for(j=0; j<3; j++)
      {
        polyslice[*curid][j] = poly[i][j];
      }
      (*curid)++;
    }
    for(i=0; i<=endid; i++)
    {
      for(j=0; j<3; j++)
      {
        polyslice[*curid][j] = poly[i][j];
      }
      (*curid)++;
    }
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

bool OnPolyNodes(real p[][3], real p0[3], int* nodeid, int nps) //predicate a point is on the vertices of a 3D polygon, get the index if on
{
  int i, j;
  real tp, v0[3];
  bool isOnNodes = false;

  for(i=0; i<nps; i++)
  {
    for(j=0; j<3; j++)
    {
      v0[j] = p0[j] - p[i][j];
    }
    tp = sqrt(DotProduct(v0,v0,3));
    if(FltEQ(tp, 0.0))
    {
      isOnNodes = true;
      *nodeid = i;
      break;
    }
  }

  return isOnNodes;
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

//boolean operations between two convex polygons
void PolyIntersect(real poly1[][3], real poly2[][3], real polyi[][3], int nps1, int nps2, int* npsi) //intersect operation of two 3D convex polygons
{
  int i, j, n, l, inext, iprev = 0, ninners, ninsecs, npoints, innerflag, insecflag, ibinsec[4];
  int innerflags[nps1], insecflags[nps1], ibinsecs[2*nps1];
  real tp, sp;
  real pa[3], pb[3], v0[3], v1[3], boundbox1[2][3], boundbox2[2][3];
  real pinsec[4][3], pinsecs[2*nps1][3], polygon[2*nps1+nps2][3];

  if (!poly1 || !poly2 || !polyi || !npsi)
  {
      printf("Unallocated memory for input or ouput polygons\n");

      return;
  }

  PolyBound(poly1, boundbox1, nps1);
  PolyBound(poly2, boundbox2, nps2);
  if (!IsBBoxIntersect(boundbox1, boundbox2))
  {
    (*npsi) = 0;
    return;
  }

  ninners = 0;
  ninsecs = 0;
  for(i=0; i<nps1; i++)
  {
    innerflags[i] = 0;
    insecflags[i] = 0;
  }

  for(n=0; n<nps1; n++)
  {
    pa[0] = poly1[n][0];
    pa[1] = poly1[n][1];
    pa[2] = poly1[n][2];
    inext = (n+1)%nps1;
    pb[0] = poly1[inext][0];
    pb[1] = poly1[inext][1];
    pb[2] = poly1[inext][2];
    innerflag = Isinner(poly2, pa, nps2);
    if(innerflag == 1)
    {
      ninners++;
      innerflags[n] = 1;
    }

    l = 0;
    for(i=0; i<nps2; i++)
    {
      inext = (i+1)%nps2;
      insecflag = InterSect(pa, pb, poly2[i], poly2[inext], pinsec[l]);
      if(insecflag == 1)
      {
        ibinsec[l] = i;
        l++;
        ninsecs++;
      }
      if(l == MAX_INTERSECTIONS) break;
    }
    if(l < 2)
      insecflags[n] = l;
    else
      insecflags[n] = 2;

    if(l == MAX_INTERSECTIONS-1)
    {
      for(j=0; j<3; j++)
      {
        v0[j] = pinsec[0][j] - pinsec[1][j];
      }
      tp = sqrt(DotProduct(v0,v0,3));
      if(FltEQ(tp, 0.0))
      {
        for(j=0; j<3; j++)
          pinsec[1][j] = pinsec[2][j];
        ibinsec[1] = ibinsec[2];
      }
    }
    else if(l == MAX_INTERSECTIONS)
    {
      for(j=0; j<3; j++)
        pinsec[1][j] = pinsec[2][j];
      ibinsec[1] = ibinsec[2];
    }

    if(l == 1)
    {
      for(j=0; j<3; j++)
      {
        pinsecs[2*n][j] = pinsec[0][j];
      }
      ibinsecs[2*n] = ibinsec[0];
      ibinsecs[2*n+1] = ibinsec[0];
    }
    else if(l > 0)
    {
      for(j=0; j<3; j++)
      {
        v0[j] = pinsec[0][j] - pa[j];
        v1[j] = pinsec[1][j] - pa[j];
      }
      tp = DotProduct(v0,v0,3);
      sp = DotProduct(v1,v1,3);
      if(tp > sp)
      {
        for(j=0; j<3; j++)
        {
          pinsecs[2*n][j] = pinsec[1][j];
          pinsecs[2*n+1][j] = pinsec[0][j];
        }
        ibinsecs[2*n] = ibinsec[1];
        ibinsecs[2*n+1] = ibinsec[0];
      }
      else
      {
        for(j=0; j<3; j++)
        {
          pinsecs[2*n][j] = pinsec[0][j];
          pinsecs[2*n+1][j] = pinsec[1][j];
        }
        ibinsecs[2*n] = ibinsec[0];
        ibinsecs[2*n+1] = ibinsec[1];
      }
    }
  }

  printf("Points of inner and intersect: %d, %d\n", ninners, ninsecs);

  npoints = 0;
  insecflag = -1;
  innerflag = -1;
  for(n=0; n<nps1; n++)
  {
    inext = (n+1)%nps1;
    if(innerflags[n] == 0 && innerflags[inext] == 0)
    {
      if(insecflags[n] > 0)
      {
        if(insecflag == -1)
        {
          for(i=1; i<nps1; i++)
          {
            iprev = (nps1+n-i)%nps1;
            if(insecflags[iprev] > 0)
            {
              insecflag = ibinsecs[2*iprev+1];
              break;
            }
          }
        }

        if(insecflag >= 0)
        {
          PolySlice(insecflag, ibinsecs[2*n], &npoints, nps2, poly2, polygon, false, false);
        }
        else
        {
          for(j=0; j<3; j++)
          {
            v0[j] = pinsecs[2*n][j] - pinsecs[2*n+1][j];
          }
          tp = sqrt(DotProduct(v0,v0,3));

          if(FltGT(tp, 0.0)) PolySlice(ibinsecs[2*n+1], ibinsecs[2*n], &npoints, nps2, poly2, polygon, false, false);
        }
      }

      if(insecflags[n] > 0) insecflag = ibinsecs[2*n+1];

      for(i=0; i<insecflags[n]; i++)
      {
        for(j=0; j<3; j++)
        {
          polygon[npoints][j] = pinsecs[2*n+i][j];
        }
        npoints++;
      }
    }
    else if(innerflags[n] == 0 && innerflags[inext] == 1)
    {
      if(insecflags[n] > 0)
      {
        if(insecflag == -1)
        {
          for(i=1; i<nps1; i++)
          {
            iprev = (nps1+n-i)%nps1;
            if(insecflags[iprev] > 0)
            {
              insecflag = ibinsecs[2*iprev+1];
              break;
            }
          }
        }

        if(insecflag >= 0) PolySlice(insecflag, ibinsecs[2*n], &npoints, nps2, poly2, polygon, false, false);
      }
      else
      {
        if(insecflag == -1)
        {
          for(i=1; i<nps1; i++)
          {
            iprev = (nps1+n-i)%nps1;
            if(insecflags[iprev] == 0 && innerflags[iprev] == 1 && innerflags[(iprev+1)%nps1] == 0)
            {
              break;
            }
          }

          OnPolyNodes(poly2, poly1[iprev], &insecflag, nps2);
        }

        OnPolyNodes(poly2, poly1[inext], &innerflag, nps2);

        if(insecflag >= 0 && innerflag >= 0) PolySlice(insecflag, innerflag, &npoints, nps2, poly2, polygon, false, false);
      }
        
      for(i=0; i<insecflags[n]; i++)
      {
        for(j=0; j<3; j++)
        {
          polygon[npoints][j] = pinsecs[2*n+i][j];
        }
        npoints++;
      }
    }
    else if(innerflags[n] == 1 && innerflags[inext] == 0)
    {
      for(j=0; j<3; j++)
      {
        polygon[npoints][j] = poly1[n][j];
      }
      npoints++;
      for(i=0; i<insecflags[n]; i++)
      {
        for(j=0; j<3; j++)
        {
          polygon[npoints][j] = pinsecs[2*n+i][j];
        }
        npoints++;
      }

      if(insecflags[n] > 0)
      {
        insecflag = ibinsecs[2*n+1];
      }
      else
      {
        OnPolyNodes(poly2, poly1[n], &insecflag, nps2);
      }
    }
    else if(innerflags[n] == 1 && innerflags[inext] == 1)
    {
      for(j=0; j<3; j++)
      {
        polygon[npoints][j] = poly1[n][j];
      }
      npoints++;
    }
  }
 
  EraseSamePointsInPoly(polygon, polyi, npoints, npsi);

  /**npsi = npoints;
  for(i=0;i<npoints;i++)
  {
    for(j=0; j<ND_ND; j++)
    {
      polyi[i][j] = polygon[i][j];
    }
  }*/
}

void PolyMerge(real poly1[][3], real poly2[][3], real polym[][3], int nps1, int nps2, int* npsm) //union operation of two 3D convex polygons
{
  int i, j, n, l, inext, iprev = 0, ninners, ninsecs, npoints, innerflag, insecflag, ibinsec[4];
  int innerflags[nps1], insecflags[nps1], ibinsecs[2*nps1];
  real tp, sp;
  real pa[3], pb[3], v0[3], v1[3], pinsec[4][3], boundbox1[2][3], boundbox2[2][3];
  real pinsecs[2*nps1][3], polygon[3*nps1+nps2][3];

  if (!poly1 || !poly2 || !polym || !npsm)
  {
      printf("Unallocated memory for input or ouput polygons\n");

      return;
  }

  PolyBound(poly1, boundbox1, nps1);
  PolyBound(poly2, boundbox2, nps2);
  if (!IsBBoxIntersect(boundbox1, boundbox2))
  {
    (*npsm) = nps1;
    for(n=0; n<nps1; n++)
    {
      for(j=0; j<3; j++)
      {
        polym[n][j] = poly1[n][j];
      }
    }
    return;
  }

  ninners = 0;
  ninsecs = 0;
  for(i=0; i<nps1; i++)
  {
    innerflags[i] = 0;
    insecflags[i] = 0;
  }

  for(n=0; n<nps1; n++)
  {
    pa[0] = poly1[n][0];
    pa[1] = poly1[n][1];
    pa[2] = poly1[n][2];
    inext = (n+1)%nps1;
    pb[0] = poly1[inext][0];
    pb[1] = poly1[inext][1];
    pb[2] = poly1[inext][2];
    innerflag = Isinner(poly2, pa, nps2);
    if(innerflag == 1)
    {
      ninners++;
      innerflags[n] = 1;
    }

    l = 0;
    for(i=0; i<nps2; i++)
    {
      inext = (i+1)%nps2;
      insecflag = InterSect(pa, pb, poly2[i], poly2[inext], pinsec[l]);
      if(insecflag == 1)
      {
        ibinsec[l] = i;
        l++;
        ninsecs++;
      }
      if(l == MAX_INTERSECTIONS) break;
    }
    if(l < 2)
      insecflags[n] = l;
    else
      insecflags[n] = 2;

    if(l == MAX_INTERSECTIONS-1)
    {
      for(j=0; j<3; j++)
      {
        v0[j] = pinsec[0][j] - pinsec[1][j];
      }
      tp = sqrt(DotProduct(v0,v0,3));
      if(FltEQ(tp, 0.0))
      {
        for(j=0; j<3; j++)
          pinsec[1][j] = pinsec[2][j];
        ibinsec[1] = ibinsec[2];
      }
    }
    else if(l == MAX_INTERSECTIONS)
    {
      for(j=0; j<3; j++)
        pinsec[1][j] = pinsec[2][j];
      ibinsec[1] = ibinsec[2];
    }

    if(l == 1)
    {
      for(j=0; j<3; j++)
      {
        pinsecs[2*n][j] = pinsec[0][j];
      }
      ibinsecs[2*n] = ibinsec[0];
      ibinsecs[2*n+1] = ibinsec[0];
    }
    else if(l > 0)
    {
      for(j=0; j<3; j++)
      {
        v0[j] = pinsec[0][j] - pa[j];
        v1[j] = pinsec[1][j] - pa[j];
      }
      tp = DotProduct(v0,v0,3);
      sp = DotProduct(v1,v1,3);
      if(tp > sp)
      {
        for(j=0; j<3; j++)
        {
          pinsecs[2*n][j] = pinsec[1][j];
          pinsecs[2*n+1][j] = pinsec[0][j];
        }
        ibinsecs[2*n] = ibinsec[1];
        ibinsecs[2*n+1] = ibinsec[0];
      }
      else
      {
        for(j=0; j<3; j++)
        {
          pinsecs[2*n][j] = pinsec[0][j];
          pinsecs[2*n+1][j] = pinsec[1][j];
        }
        ibinsecs[2*n] = ibinsec[0];
        ibinsecs[2*n+1] = ibinsec[1];
      }
    }
  }

  printf("Points of inner and intersect: %d, %d\n", ninners, ninsecs);

  if(ninners == nps1 && ninsecs == 0)
  {
    *npsm = nps2;
    for(n=0; n<nps2; n++)
    {
      for(j=0; j<3; j++)
      {
        polym[n][j] = poly2[n][j];
      }
    }
    return;
  }

  npoints = 0;
  insecflag = -1;
  innerflag = -1;
  for(n=0; n<nps1; n++)
  {
    inext = (n+1)%nps1;
    if(innerflags[n] == 0 && innerflags[inext] == 0)
    {
      if(insecflags[n] == 1 && ninners == 0)
      {
        for(i=1; i<nps1; i++)
        {
          iprev = (nps1+n-i)%nps1;
          if(insecflags[iprev] == 1 && innerflags[iprev] == 0 && innerflags[(iprev+1)%nps1] == 0)
          {
            break;
          }
        }

        for(j=0; j<3; j++)
        {
          v0[j] = pinsecs[2*n][j] - poly1[n][j];
        }
        tp = sqrt(DotProduct(v0,v0,3));

        if(FltEQ(tp, 0.0)) PolySlice(ibinsecs[2*iprev], ibinsecs[2*n], &npoints, nps2, poly2, polygon, false, true);
      }

      for(j=0; j<3; j++)
      {
        polygon[npoints][j] = poly1[n][j];
      }
      npoints++;

      if(insecflags[n] == 1 && ninners == 0)
      {
        for(j=0; j<3; j++)
        {
          v0[j] = pinsecs[2*n][j] - poly1[inext][j];
        }
        tp = sqrt(DotProduct(v0,v0,3));

        if(FltEQ(tp, 0.0))
        {
          for(j=0; j<3; j++)
          {
            polygon[npoints][j] = pinsecs[2*n][j];
          }
          npoints++;
        }
      }
      
      if(insecflags[n] == 2)
      {
        for(j=0; j<3; j++)
        {
          polygon[npoints][j] = pinsecs[2*n][j];
        }
        npoints++;

        for(j=0; j<3; j++)
        {
          v0[j] = pinsecs[2*n][j] - pinsecs[2*n+1][j];
        }
        tp = sqrt(DotProduct(v0,v0,3));

        if(FltGT(tp, 0.0))
        {
          PolySlice(ibinsecs[2*n], ibinsecs[2*n+1], &npoints, nps2, poly2, polygon, false, false);
        }
        else if(ninners==0 && ninsecs == 2)
        {
          if(ibinsecs[2*n+1] > ibinsecs[2*n])
            PolySlice(ibinsecs[2*n+1], ibinsecs[2*n], &npoints, nps2, poly2, polygon, false, false);
          else
            PolySlice(ibinsecs[2*n], ibinsecs[2*n+1], &npoints, nps2, poly2, polygon, false, false);
        }

        for(j=0; j<3; j++)
        {
          polygon[npoints][j] = pinsecs[2*n+1][j];
        }
        npoints++;
      }
    }
    else if(innerflags[n] == 0 && innerflags[inext] == 1)
    {
      for(j=0; j<3; j++)
      {
        polygon[npoints][j] = poly1[n][j];
      }
      npoints++;

      for(i=0; i<insecflags[n]; i++)
      {
        for(j=0; j<3; j++)
        {
          polygon[npoints][j] = pinsecs[2*n+i][j];
        }
        npoints++;
      }

      if(insecflags[n] > 0)
      {
        insecflag = ibinsecs[2*n+1];
      }
      else
      {
        OnPolyNodes(poly2, poly1[inext], &insecflag, nps2);
      }
    }
    else if(innerflags[n] == 1 && innerflags[inext] == 0)
    {
      if(insecflags[n] > 0)
      {
        if(insecflag == -1)
        {
          for(i=1; i<nps1; i++)
          {
            iprev = (nps1+n-i)%nps1;
            if(insecflags[iprev] > 0)
            {
              insecflag = ibinsecs[2*iprev+1];
              break;
            }
          }
        }

        if(insecflag >= 0) PolySlice(insecflag, ibinsecs[2*n], &npoints, nps2, poly2, polygon, false, true);
      }
      else
      {
        if(insecflag == -1)
        {
          for(i=1; i<nps1; i++)
          {
            iprev = (nps1+n-i)%nps1;
            if(insecflags[iprev] == 0 && innerflags[iprev] == 0 && innerflags[(iprev+1)%nps1] == 1)
            {
              break;
            }
          }

          OnPolyNodes(poly2, poly1[(iprev+1)%nps1], &insecflag, nps2);
        }

        OnPolyNodes(poly2, poly1[n], &innerflag, nps2);

        if(insecflag >= 0 && innerflag >= 0) PolySlice(insecflag, innerflag, &npoints, nps2, poly2, polygon, true, true);
      }

      for(i=0; i<insecflags[n]; i++)
      {
        for(j=0; j<3; j++)
        {
          polygon[npoints][j] = pinsecs[2*n+i][j];
        }
        npoints++;
      }
    }
  }
 
  EraseSamePointsInPoly(polygon, polym, npoints, npsm);

  /**npsm = npoints;
  for(i=0;i<npoints;i++)
  {
    for(j=0; j<ND_ND; j++)
    {
      polym[i][j] = polygon[i][j];
    }
  }*/
}

bool EraseSamePointsInPoly(real poly[][ND_ND], real polyn[][ND_ND], int nps, int* npsn) //erase repeated points in one polygon
{
  int inext;
  bool hasSamePoints, allSamePoints;
  real dis;
  real v[ND_ND];

  *npsn = 0;
  hasSamePoints = false;
  allSamePoints = true;
  
  for(int i=0; i<nps; i++)
  {
    inext = (i+1)%nps;
    for(int j=0; j<ND_ND; j++)
    {   
      v[j] = poly[inext][j] - poly[i][j];
    }
    dis = sqrt(DotProduct(v,v,ND_ND));
    if(FltEQ(dis, 0.0))
    {
      hasSamePoints = true;
    }
    else
    {
      allSamePoints = false;
      for(int j=0; j<ND_ND; j++)
        polyn[*npsn][j] = poly[i][j];
      (*npsn)++;
    }
  }

  if(nps > 0 && allSamePoints)
  {
    for(int j=0; j<ND_ND; j++)
      polyn[*npsn][j] = poly[0][j];
    (*npsn)++;
  }

  return hasSamePoints;
}

// -------------------------------------------------------------
