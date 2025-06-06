#pragma once

#include "geometry3d_basic.h"

void PolySlice(int slicebegin, int slicend, int* curid, int nps, real poly[][3], real polyslice[][3], bool startin, bool outloop);
bool OnPolyNodes(real p[][3], real p0[3], int* nodeid, int nps);
void BuildBooleanTopo(real poly1[][3], real poly2[][3], int nps1, int nps2, int* innerflags, int* insecflags, int* ibinsecs, real pinsecs[][3], int* innersnum, int* insecsnum);
void PolyIntersect(real poly1[][3], real poly2[][3], real polyi[][3], int nps1, int nps2, int* npsi);
void PolyMerge(real poly1[][3], real poly2[][3], real polym[][3], int nps1, int nps2, int* npsm);
bool EraseSamePointsInPoly(real poly[][3], real polyn[][3], int nps, int* npsn);