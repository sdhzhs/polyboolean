#include <stdio.h>
#include "convex_poly_boolean.h"

#define ND_ND 3
#define MAX_INTERSECTIONS 4

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

void BuildBooleanTopo(real poly1[][3], real poly2[][3], int nps1, int nps2, int* innerflags, int* insecflags, int* ibinsecs, real pinsecs[][3], int* innersnum, int* insecsnum)
{
  int i, j, n, l, inext, ninners, ninsecs, innerflag, insecflag, ibinsec[MAX_INTERSECTIONS];
  real tp, sp;
  real pa[3], pb[3], v0[3], v1[3], pinsec[MAX_INTERSECTIONS][3];

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
        pinsecs[2*n+1][j] = pinsec[0][j];
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
      if(FltGT(tp, sp) || (FltEQ(tp, sp) && ibinsec[0] == 0 && ibinsec[1] == nps2-1))
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

  (*innersnum) = ninners;
  (*insecsnum) = ninsecs;
}

//boolean operations between two convex polygons
void PolyIntersect(real poly1[][3], real poly2[][3], real polyi[][3], int nps1, int nps2, int* npsi) //intersect operation of two 3D convex polygons
{
  int i, j, n, inext, iprev = 0, ninners, ninsecs, npoints, innerflag, insecflag;
  int innerflags[nps1], insecflags[nps1], ibinsecs[2*nps1];
  real tp, pa[3], pb[3], v0[3], v1[3], boundbox1[2][3], boundbox2[2][3];
  real pinsecs[2*nps1][3], polygon[2*nps1+nps2][3];

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

  PolyProp(poly1, pa, v0, nps1);
  PolyProp(poly2, pb, v1, nps2);
  if(FltLT(DotProduct(v0,v1,3), 0.0))
  {
    for(i=0; i<nps2/2; i++)
    {
      for(j=0; j<3; j++)
      {
        pa[j] = poly2[nps2-1-i][j];
        poly2[nps2-1-i][j] = poly2[i][j];
        poly2[i][j] = pa[j];
      }
    }
  }

  BuildBooleanTopo(poly1, poly2, nps1, nps2, innerflags, insecflags, ibinsecs, pinsecs, &ninners, &ninsecs);

  if(ninners == 0 && ninsecs == 0 && Isinner(poly1, poly2[0], nps1))
  {
    *npsi = nps2;
    for(n=0; n<nps2; n++)
    {
      for(j=0; j<3; j++)
      {
        polyi[n][j] = poly2[n][j];
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
      if(insecflags[n] == 2)
      {
        if(insecflag == -1)
        {
          for(i=1; i<nps1; i++)
          {
            iprev = (nps1+n-i)%nps1;
            if(insecflags[iprev] == 2 || (insecflags[iprev] == 1 && innerflags[iprev] == 1 && innerflags[(iprev+1)%nps1] == 0))
            {
              insecflag = ibinsecs[2*iprev+1];
              break;
            }
            else if(insecflags[iprev] == 0 && innerflags[iprev] == 1 && innerflags[(iprev+1)%nps1] == 0)
            {
              OnPolyNodes(poly2, poly1[iprev], &insecflag, nps2);
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

          if(FltGT(tp, 0.0))
            PolySlice(ibinsecs[2*n+1], ibinsecs[2*n], &npoints, nps2, poly2, polygon, false, false);
          else
          {
            for(j=0; j<3; j++)
            {
              v0[j] = poly1[inext][j] - poly1[n][j];
              v1[j] = pinsecs[2*n][j] - poly2[ibinsecs[2*n]][j];
            }
            if(FltGT(DotProduct(v0,v1,3),0.0)) PolySlice(ibinsecs[2*n+1], ibinsecs[2*n], &npoints, nps2, poly2, polygon, false, false);
          }
        }
      }

      if(insecflags[n] == 2) insecflag = ibinsecs[2*n+1];

      if(insecflags[n] == 2 || (insecflags[n] == 1 && ((ninners == 0 && ninsecs == 2) || ninsecs == 1)))
      {
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
            else if(innerflags[iprev] == 1 && innerflags[(iprev+1)%nps1] == 0)
            {
              OnPolyNodes(poly2, poly1[iprev], &insecflag, nps2);
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
            if(insecflags[iprev] > 0)
            {
              insecflag = ibinsecs[2*iprev+1];
              break;
            }
            else if(innerflags[iprev] == 1 && innerflags[(iprev+1)%nps1] == 0)
            {
              OnPolyNodes(poly2, poly1[iprev], &insecflag, nps2);
              break;
            }
          }
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
  int i, j, n, inext, iprev = 0, ninners, ninsecs, npoints, innerflag, insecflag;
  int innerflags[nps1], insecflags[nps1], ibinsecs[2*nps1];
  real tp, sp, pa[3], pb[3], v0[3], v1[3], boundbox1[2][3], boundbox2[2][3];
  real pinsecs[2*nps1][3], polygon[3*nps1+nps2][3];
  bool found;

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

  PolyProp(poly1, pa, v0, nps1);
  PolyProp(poly2, pb, v1, nps2);
  if(FltLT(DotProduct(v0,v1,3), 0.0))
  {
    for(i=0; i<nps2/2; i++)
    {
      for(j=0; j<3; j++)
      {
        pa[j] = poly2[nps2-1-i][j];
        poly2[nps2-1-i][j] = poly2[i][j];
        poly2[i][j] = pa[j];
      }
    }
  }

  BuildBooleanTopo(poly1, poly2, nps1, nps2, innerflags, insecflags, ibinsecs, pinsecs, &ninners, &ninsecs);

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
      found = false;
      if(insecflags[n] == 1 && ninners == 0 && ninsecs == 2)
      {
        for(j=0; j<3; j++)
        {
          v0[j] = pinsecs[2*n][j] - poly1[n][j];
        }
        tp = sqrt(DotProduct(v0,v0,3));
        for(j=0; j<3; j++)
        {
          v0[j] = pinsecs[2*n][j] - poly2[(ibinsecs[2*n]+1)%nps2][j];
        }
        sp = sqrt(DotProduct(v0,v0,3));

        if(FltEQ(sp, 0.0)) found = true;

        if(FltEQ(tp, 0.0) || FltEQ(sp, 0.0))
        {
          for(i=1; i<nps1; i++)
          {
            iprev = (nps1+n-i)%nps1;
            if(insecflags[iprev] == 1 && innerflags[iprev] == 0 && innerflags[(iprev+1)%nps1] == 0)
            {
              break;
            }
          }
          PolySlice(ibinsecs[2*iprev], ibinsecs[2*n], &npoints, nps2, poly2, polygon, false, true);
        }
      }
      
      if(!found)
      {
        for(j=0; j<3; j++)
        {
          polygon[npoints][j] = poly1[n][j];
        }
        npoints++;
      }

      if(insecflags[n] == 1 && ninners == 0 && ninsecs == 2)
      {
        for(j=0; j<3; j++)
        {
          polygon[npoints][j] = pinsecs[2*n][j];
        }
        npoints++;
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
        else if(ninners == 0 && ninsecs == 2)
        {
          for(j=0; j<3; j++)
          {
            v0[j] = poly1[inext][j] - poly1[n][j];
            v1[j] = pinsecs[2*n][j] - poly2[ibinsecs[2*n]][j];
          }
          if(FltLT(DotProduct(v0,v1,3),0.0)) PolySlice(ibinsecs[2*n+1], ibinsecs[2*n], &npoints, nps2, poly2, polygon, false, false);
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
        for(j=0; j<3; j++)
        {
          polygon[npoints][j] = poly1[inext][j];
        }
        npoints++;
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
            else if(innerflags[iprev] == 0 && innerflags[(iprev+1)%nps1] == 1)
            {
              OnPolyNodes(poly2, poly1[(iprev+1)%nps1], &insecflag, nps2);
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
            if(insecflags[iprev] > 0)
            {
              insecflag = ibinsecs[2*iprev+1];
              break;
            }
            else if(innerflags[iprev] == 0 && innerflags[(iprev+1)%nps1] == 1)
            {
              OnPolyNodes(poly2, poly1[(iprev+1)%nps1], &insecflag, nps2);
              break;
            }
          }
        }

        OnPolyNodes(poly2, poly1[n], &innerflag, nps2);

        for(j=0; j<3; j++)
        {
          v0[j] = poly1[inext][j] - poly1[n][j];
          v1[j] = poly2[innerflag][j] - poly2[(innerflag-1+nps2)%nps2][j];
        }

        CrossProduct(v0,v1,pa);
        tp = DotProduct(pa,pa,3);
        sp = DotProduct(v0,v1,3);
        if(FltEQ(tp, 0.0) && FltLT(sp, 0.0)) innerflag = (innerflag-1+nps2)%nps2;

        if(insecflag >= 0 && innerflag >= 0) PolySlice(insecflag, innerflag, &npoints, nps2, poly2, polygon, false, true);
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