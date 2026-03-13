#include <stdio.h>
#include "convex_poly_boolean.h"

int main (int argc, char *argv[])
{
   int i, j, k, npoints;
   real poly1[6][3], poly2[6][3], polyi[18][3];

   for(i=0; i<3; i++)
   {
     poly1[i][2] = 0.0;
     poly2[i][2] = 0.0;
   }

   for(i=0; i<9; i++)
   {
     polyi[i][0] = 0.0;
     polyi[i][1] = 0.0;
     polyi[i][2] = 0.0;
   }

   poly1[0][0] = -1.0;
   poly1[0][1] = 0.0;
   poly1[1][0] = 1.0;
   poly1[1][1] = 0.0;
   poly1[2][0] = 0.0;
   poly1[2][1] = 1.0;

   //intersection and union of two congruent triangles, horizontal movement
   printf("Case: 1\n");
   for(j=-1; j<4; j++)
   {
     printf("Move: %d\n", j);

     poly2[0][0] = -3.0+j*1.0;
     poly2[0][1] = 0.0;
     poly2[1][0] = -1.0+j*1.0;
     poly2[1][1] = 0.0;
     poly2[2][0] = -2.0+j*1.0;
     poly2[2][1] = 1.0;

     PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }

   printf("-----------------------------\n");

   //intersection and union of two congruent triangles, vertical movement
   printf("Case: 2\n");
   for(j=0; j<5; j++)
   {
     printf("Move: %d\n", j);

     poly2[0][0] = -1.0;
     poly2[0][1] = -1.0+j*0.5;
     poly2[1][0] = 1.0;
     poly2[1][1] = -1.0+j*0.5;
     poly2[2][0] = 0.0;
     poly2[2][1] = j*0.5;

     PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }

   printf("-----------------------------\n");

   //intersection and union of two symmetric triangles, horizontal movement
   printf("Case: 3\n");
   for(j=0; j<3; j++)
   {
     printf("Move: %d\n", j);

     poly2[0][0] = -2.0+j*1.0;
     poly2[0][1] = 1.0;
     poly2[1][0] = -1.0+j*1.0;
     poly2[1][1] = 0.0;
     poly2[2][0] = j*1.0;
     poly2[2][1] = 1.0;

     PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }

   printf("-----------------------------\n");

   //intersection and union of two symmetric triangles, vertical movement
   printf("Case: 4\n");
   for(j=0; j<5; j++)
   {
     printf("Move: %d\n", j);

     poly2[0][0] = -1.0;
     poly2[0][1] = j*0.5;
     poly2[1][0] = 0.0;
     poly2[1][1] = -1.0+j*0.5;
     poly2[2][0] = 1.0;
     poly2[2][1] = j*0.5;

     PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }

   printf("-----------------------------\n");

   //intersection and union of two colinear different triangles, horizontal movement
   printf("Case: 5\n");
   for(j=0; j<5; j++)
   {
     printf("Move: %d\n", j);

     poly2[0][0] = -3.0+j*1.0;
     poly2[0][1] = 0.0;
     poly2[1][0] = -2.0+j*1.0;
     poly2[1][1] = -1.0;
     poly2[2][0] = -1.0+j*1.0;
     poly2[2][1] = 0.0;
  
     PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);
  
     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
  
     PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);
  
     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
  
     PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);
  
     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
  
     PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);
  
     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }

   printf("-----------------------------\n");

   //intersection and union of two colinear different triangles, vertical movement
   printf("Case: 6\n");
   for(j=0; j<4; j++)
   {
     printf("Move: %d\n", j);

     poly2[0][0] = -1.0;
     poly2[0][1] = j*0.5;
     poly2[1][0] = 1.0;
     poly2[1][1] = -1.5+j*0.5;
     poly2[2][0] = 1.0;
     poly2[2][1] = j*0.5;
  
     PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);
  
     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
  
     PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);
  
     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
  
     PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);
  
     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
  
     PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);
  
     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }

   printf("-----------------------------\n");

   //intersection and union of two node-on-edge different triangles, vertical movement
   printf("Case: 7\n");
   for(j=0; j<8; j++)
   {
     printf("Move: %d\n", j);

     poly2[0][0] = -0.4;
     poly2[0][1] = j*0.2;
     poly2[1][0] = 0.0;
     poly2[1][1] = -0.4+j*0.2;
     poly2[2][0] = 0.4;
     poly2[2][1] = j*0.2;
  
     PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);
  
     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
  
     PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);
  
     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
  
     PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);
  
     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
  
     PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);
  
     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }

   printf("-----------------------------\n");

   //intersection and union of two node-on-edge different triangles, horizontal movement
   printf("Case: 8\n");
   for(j=0; j<5; j++)
   {
     printf("Move: %d\n", j);

     poly2[0][0] = -1.5+j*0.5;
     poly2[0][1] = 0.5;
     poly2[1][0] = -0.5+j*0.5;
     poly2[1][1] = 0.5;
     poly2[2][0] = -1.0+j*0.5;
     poly2[2][1] = 1.5;
  
     PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);
  
     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
  
     PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);
  
     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
  
     PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);
  
     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
  
     PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);
  
     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }

   printf("-----------------------------\n");

   //intersection and union of two congruent triangles, horizontal movement, opposite points order
   printf("Case: 9\n");
   for(j=-1; j<4; j++)
   {
     printf("Move: %d\n", j);

     poly2[2][0] = -3.0+j*1.0;
     poly2[2][1] = 0.0;
     poly2[1][0] = -1.0+j*1.0;
     poly2[1][1] = 0.0;
     poly2[0][0] = -2.0+j*1.0;
     poly2[0][1] = 1.0;

     PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }

   printf("-----------------------------\n");

   //intersection and union of two similar triangles, small sized one horizontal movement
   printf("Case: 10\n");
   for(j=0; j<5; j++)
   {
     printf("Move: %d\n", j);

     poly2[0][0] = -2.0+j*0.8;
     poly2[0][1] = 0.3;
     poly2[1][0] = -1.2+j*0.8;
     poly2[1][1] = 0.3;
     poly2[2][0] = -1.6+j*0.8;
     poly2[2][1] = 0.7;

     PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }

   printf("-----------------------------\n");

   //intersection and union of two congruent triangles, horizontal movement, on yoz plane
   printf("Case: 11\n");
   real p0[3] = {-5.0, 0.0, 0.0}, n[3] = {1.0, 0.0, 0.0}, t[3] = {-1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0)};
   
   for(i=0; i<3; i++)
   {
     PntShadow(p0, n, t, poly1[i], polyi[i]);

     for(j=0; j<3; j++)
       poly1[i][j] = polyi[i][j];
   }

   for(j=-1; j<4; j++)
   {
     printf("Move: %d\n", j);

     poly2[0][0] = -3.0+j*1.0;
     poly2[0][1] = 0.0;
     poly2[1][0] = -1.0+j*1.0;
     poly2[1][1] = 0.0;
     poly2[2][0] = -2.0+j*1.0;
     poly2[2][1] = 1.0;

     for(i=0; i<3; i++)
     {
       poly2[i][2] = 0.0;
     }

     for(i=0; i<3; i++)
     {
       PntShadow(p0, n, t, poly2[i], polyi[i]);

       for(k=0; k<3; k++)
         poly2[i][k] = polyi[i][k];
     }

     PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }

   printf("-----------------------------\n");

   for(i=0; i<6; i++)
   {
     poly1[i][2] = 0.0;
   }
   for(i=0; i<4; i++)
   {
     poly2[i][2] = 0.0;
   }

   for(i=0; i<16; i++)
   {
     polyi[i][0] = 0.0;
     polyi[i][1] = 0.0;
     polyi[i][2] = 0.0;
   }

   poly1[0][0] = -0.5;
   poly1[0][1] = 0.5;
   poly1[1][0] = -1.0;
   poly1[1][1] = 0.0;
   poly1[2][0] = -0.5;
   poly1[2][1] = -0.5;
   poly1[3][0] = 0.5;
   poly1[3][1] = -0.5;
   poly1[4][0] = 1.0;
   poly1[4][1] = 0.0;
   poly1[5][0] = 0.5;
   poly1[5][1] = 0.5;

   //intersection and union between one static hexagon and one moving square, vertical movement
   printf("Case: 12\n");
   for(j=0; j<7; j++)
   {
     printf("Move: %d\n", j);

     poly2[0][0] = -1.0;
     poly2[0][1] = -0.5+j*0.5;
     poly2[1][0] = -1.0;
     poly2[1][1] = -2.5+j*0.5;
     poly2[2][0] = 1.0;
     poly2[2][1] = -2.5+j*0.5;
     poly2[3][0] = 1.0;
     poly2[3][1] = -0.5+j*0.5;

     PolyIntersect(poly1, poly2, polyi, 6, 4, &npoints);
     
     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyIntersect(poly2, poly1, polyi, 4, 6, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly1, poly2, polyi, 6, 4, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly2, poly1, polyi, 4, 6, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }

   printf("-----------------------------\n");

   //intersection and union between one static hexagon and one moving square, horizontal movement
   printf("Case: 13\n");
   for(j=0; j<9; j++)
   {
     printf("Move: %d\n", j);

     poly2[0][0] = -3.0+j*0.5;
     poly2[0][1] = 1.0;
     poly2[1][0] = -3.0+j*0.5;
     poly2[1][1] = -1.0;
     poly2[2][0] = -1.0+j*0.5;
     poly2[2][1] = -1.0;
     poly2[3][0] = -1.0+j*0.5;
     poly2[3][1] = 1.0;

     PolyIntersect(poly1, poly2, polyi, 6, 4, &npoints);
     
     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyIntersect(poly2, poly1, polyi, 4, 6, &npoints);

     printf("Points of intersected polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Intersected polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly1, poly2, polyi, 6, 4, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }

     PolyMerge(poly2, poly1, polyi, 4, 6, &npoints);

     printf("Points of merged polygon: %d\n", npoints);
     for(i=0; i<npoints; i++)
     {
       printf("Merged polygon: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
     }
   }


   //Case 14: Two polygons touching at a single vertex
   printf("Case: 14 - Touching at Single Vertex\n");
   
   poly1[0][0] = -1.0;
   poly1[0][1] = 0.0;
   poly1[0][2] = 0.0;
   poly1[1][0] = 1.0;
   poly1[1][1] = 0.0;
   poly1[1][2] = 0.0;
   poly1[2][0] = 0.0;
   poly1[2][1] = 1.5;
   poly1[2][2] = 0.0;

   poly2[0][0] = 0.0;
   poly2[0][1] = 1.5;
   poly2[0][2] = 0.0;
   poly2[1][0] = 1.0;
   poly2[1][1] = 1.5;
   poly2[1][2] = 0.0;
   poly2[2][0] = 0.5;
   poly2[2][1] = 2.5;
   poly2[2][2] = 0.0;

   PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Intersection (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Intersection (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   printf("-----------------------------\n");

   //Case 15: Two polygons sharing an edge
   printf("Case: 15 - Sharing an Edge\n");
   
   poly1[0][0] = -1.0;
   poly1[0][1] = -1.0;
   poly1[0][2] = 0.0;
   poly1[1][0] = 1.0;
   poly1[1][1] = -1.0;
   poly1[1][2] = 0.0;
   poly1[2][0] = 0.0;
   poly1[2][1] = 1.0;
   poly1[2][2] = 0.0;

   poly2[0][0] = -1.0;
   poly2[0][1] = -1.0;
   poly2[0][2] = 0.0;
   poly2[1][0] = 1.0;
   poly2[1][1] = -1.0;
   poly2[1][2] = 0.0;
   poly2[2][0] = 0.0;
   poly2[2][1] = -2.0;
   poly2[2][2] = 0.0;

   PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Intersection (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Intersection (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }


   //Case 16: Two identical polygons (coincident)
   printf("Case: 16 - Identical/Coincident Polygons\n");
   
   poly1[0][0] = -1.0;
   poly1[0][1] = -1.0;
   poly1[0][2] = 0.0;
   poly1[1][0] = 1.0;
   poly1[1][1] = -1.0;
   poly1[1][2] = 0.0;
   poly1[2][0] = 0.0;
   poly1[2][1] = 1.0;
   poly1[2][2] = 0.0;

   poly2[0][0] = -1.0;
   poly2[0][1] = -1.0;
   poly2[0][2] = 0.0;
   poly2[1][0] = 1.0;
   poly2[1][1] = -1.0;
   poly2[1][2] = 0.0;
   poly2[2][0] = 0.0;
   poly2[2][1] = 1.0;
   poly2[2][2] = 0.0;

   PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Intersection (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Intersection (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }


   //Case 17: Vertex of one polygon on edge interior of another
   printf("Case: 17 - Vertex on Edge Interior\n");
   
   poly1[0][0] = -3.0;
   poly1[0][1] = -1.0;
   poly1[0][2] = 0.0;
   poly1[1][0] = 3.0;
   poly1[1][1] = -1.0;
   poly1[1][2] = 0.0;
   poly1[2][0] = 0.0;
   poly1[2][1] = 2.0;
   poly1[2][2] = 0.0;

   poly2[0][0] = 0.0;
   poly2[0][1] = -1.0;
   poly2[0][2] = 0.0;
   poly2[1][0] = 1.5;
   poly2[1][1] = 0.5;
   poly2[1][2] = 0.0;
   poly2[2][0] = -1.5;
   poly2[2][1] = 0.5;
   poly2[2][2] = 0.0;

   PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Intersection (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Intersection (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   printf("-----------------------------\n");

   //Case 18: Partial collinear edge overlap
   printf("Case: 18 - Partial Collinear Edge Overlap\n");
   
   poly1[0][0] = -2.0;
   poly1[0][1] = -1.0;
   poly1[0][2] = 0.0;
   poly1[1][0] = 2.0;
   poly1[1][1] = -1.0;
   poly1[1][2] = 0.0;
   poly1[2][0] = 0.0;
   poly1[2][1] = 2.0;
   poly1[2][2] = 0.0;

   poly2[0][0] = -1.0;
   poly2[0][1] = -1.0;
   poly2[0][2] = 0.0;
   poly2[1][0] = 3.0;
   poly2[1][1] = -1.0;
   poly2[1][2] = 0.0;
   poly2[2][0] = 1.0;
   poly2[2][1] = 2.0;
   poly2[2][2] = 0.0;

   PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Intersection (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Intersection (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   printf("-----------------------------\n");

   //Case 19: Multiple vertices on edges - triangle and quadrilateral
   printf("Case: 19 - Multiple Vertices on Edges (Triangle-Quadrilateral)\n");
   
   poly1[0][0] = -1.5;
   poly1[0][1] = -1.0;
   poly1[0][2] = 0.0;
   poly1[1][0] = 1.5;
   poly1[1][1] = -1.0;
   poly1[1][2] = 0.0;
   poly1[2][0] = 0.0;
   poly1[2][1] = 1.5;
   poly1[2][2] = 0.0;

   poly2[0][0] = -1.0;
   poly2[0][1] = -0.5;
   poly2[0][2] = 0.0;
   poly2[1][0] = 1.0;
   poly2[1][1] = -0.5;
   poly2[1][2] = 0.0;
   poly2[2][0] = 1.0;
   poly2[2][1] = 0.8;
   poly2[2][2] = 0.0;
   poly2[3][0] = -1.0;
   poly2[3][1] = 0.8;
   poly2[3][2] = 0.0;

   PolyIntersect(poly1, poly2, polyi, 3, 4, &npoints);
   printf("Intersection (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyIntersect(poly2, poly1, polyi, 4, 3, &npoints);
   printf("Intersection (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly1, poly2, polyi, 3, 4, &npoints);
   printf("Merged polygon (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly2, poly1, polyi, 4, 3, &npoints);
   printf("Merged polygon (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }


   //Case 20: Two Quadrilaterals with spiral crossing pattern
   printf("Case: 20 - Two Quadrilaterals Spiral Pattern\n");
   
   poly1[0][0] = -1.0;
   poly1[0][1] = -1.5;
   poly1[0][2] = 0.0;
   poly1[1][0] = 1.5;
   poly1[1][1] = -0.5;
   poly1[1][2] = 0.0;
   poly1[2][0] = 1.0;
   poly1[2][1] = 1.5;
   poly1[2][2] = 0.0;
   poly1[3][0] = -1.5;
   poly1[3][1] = 0.5;
   poly1[3][2] = 0.0;

   poly2[0][0] = -1.5;
   poly2[0][1] = -0.5;
   poly2[0][2] = 0.0;
   poly2[1][0] = 1.0;
   poly2[1][1] = -1.0;
   poly2[1][2] = 0.0;
   poly2[2][0] = 1.5;
   poly2[2][1] = 1.0;
   poly2[2][2] = 0.0;
   poly2[3][0] = -0.5;
   poly2[3][1] = 1.5;
   poly2[3][2] = 0.0;

   PolyIntersect(poly1, poly2, polyi, 4, 4, &npoints);
   printf("Intersection (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyIntersect(poly2, poly1, polyi, 4, 4, &npoints);
   printf("Intersection (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly1, poly2, polyi, 4, 4, &npoints);
   printf("Merged polygon (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly2, poly1, polyi, 4, 4, &npoints);
   printf("Merged polygon (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   printf("-----------------------------\n");

   //Case 21: Hexagon-Triangle with vertex on edge and deep overlap
   printf("Case: 21 - Hexagon-Triangle Overlapping with Edge Contact\n");
   
   poly1[0][0] = -1.5;
   poly1[0][1] = -1.0;
   poly1[0][2] = 0.0;
   poly1[1][0] = 1.5;
   poly1[1][1] = -1.0;
   poly1[1][2] = 0.0;
   poly1[2][0] = 2.0;
   poly1[2][1] = 0.5;
   poly1[2][2] = 0.0;
   poly1[3][0] = 1.0;
   poly1[3][1] = 1.5;
   poly1[3][2] = 0.0;
   poly1[4][0] = -1.0;
   poly1[4][1] = 1.5;
   poly1[4][2] = 0.0;
   poly1[5][0] = -2.0;
   poly1[5][1] = 0.5;
   poly1[5][2] = 0.0;

   poly2[0][0] = -1.0;
   poly2[0][1] = -0.5;
   poly2[0][2] = 0.0;
   poly2[1][0] = 1.5;
   poly2[1][1] = 0.5;
   poly2[1][2] = 0.0;
   poly2[2][0] = -0.25;
   poly2[2][1] = 1.5;
   poly2[2][2] = 0.0;

   PolyIntersect(poly1, poly2, polyi, 6, 3, &npoints);
   printf("Intersection (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyIntersect(poly2, poly1, polyi, 3, 6, &npoints);
   printf("Intersection (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly1, poly2, polyi, 6, 3, &npoints);
   printf("Merged polygon (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly2, poly1, polyi, 3, 6, &npoints);
   printf("Merged polygon (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   printf("-----------------------------\n");

   //Case 22: Two edges crossing at interior points (edge-to-edge X pattern)
   printf("Case: 22 - Two Edges Crossing at Interior Points\n");
   
   poly1[0][0] = -1.0;
   poly1[0][1] = -1.0;
   poly1[0][2] = 0.0;
   poly1[1][0] = 2.0;
   poly1[1][1] = 1.0;
   poly1[1][2] = 0.0;
   poly1[2][0] = 0.0;
   poly1[2][1] = 2.0;
   poly1[2][2] = 0.0;

   poly2[0][0] = 0.0;
   poly2[0][1] = -1.0;
   poly2[0][2] = 0.0;
   poly2[1][0] = 1.5;
   poly2[1][1] = 0.5;
   poly2[1][2] = 0.0;
   poly2[2][0] = 0.0;
   poly2[2][1] = 1.5;
   poly2[2][2] = 0.0;

   PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Intersection (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Intersection (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   printf("-----------------------------\n");

   //Case 23: Inner polygon with one edge tangent to outer polygon boundary
   printf("Case: 23 - Inner Polygon with Edge Tangent to Boundary\n");
   
   poly1[0][0] = -2.0;
   poly1[0][1] = -2.0;
   poly1[0][2] = 0.0;
   poly1[1][0] = 2.0;
   poly1[1][1] = -2.0;
   poly1[1][2] = 0.0;
   poly1[2][0] = 2.0;
   poly1[2][1] = 2.0;
   poly1[2][2] = 0.0;
   poly1[3][0] = -2.0;
   poly1[3][1] = 2.0;
   poly1[3][2] = 0.0;

   poly2[0][0] = -1.0;
   poly2[0][1] = 0.0;
   poly2[0][2] = 0.0;
   poly2[1][0] = 1.0;
   poly2[1][1] = 0.0;
   poly2[1][2] = 0.0;
   poly2[2][0] = 0.0;
   poly2[2][1] = 1.0;
   poly2[2][2] = 0.0;

   PolyIntersect(poly1, poly2, polyi, 4, 3, &npoints);
   printf("Intersection (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyIntersect(poly2, poly1, polyi, 3, 4, &npoints);
   printf("Intersection (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly1, poly2, polyi, 4, 3, &npoints);
   printf("Merged polygon (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly2, poly1, polyi, 3, 4, &npoints);
   printf("Merged polygon (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   printf("-----------------------------\n");

   //Case 24: Partial collinear edge overlap with one vertex inside
   printf("Case: 24 - Collinear Edge Overlap with Vertex Inside\n");
   
   poly1[0][0] = -1.5;
   poly1[0][1] = 0.0;
   poly1[0][2] = 0.0;
   poly1[1][0] = 1.5;
   poly1[1][1] = 0.0;
   poly1[1][2] = 0.0;
   poly1[2][0] = 0.0;
   poly1[2][1] = 1.5;
   poly1[2][2] = 0.0;

   poly2[0][0] = -0.5;
   poly2[0][1] = 0.0;
   poly2[0][2] = 0.0;
   poly2[1][0] = 2.0;
   poly2[1][1] = 0.0;
   poly2[1][2] = 0.0;
   poly2[2][0] = 0.5;
   poly2[2][1] = 1.0;
   poly2[2][2] = 0.0;

   PolyIntersect(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Intersection (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyIntersect(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Intersection (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly1, poly2, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly2, poly1, polyi, 3, 3, &npoints);
   printf("Merged polygon (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   printf("-----------------------------\n");

   //Case 25: Quadrilateral-Triangle T-junction (vertex on edge creating T pattern)
   printf("Case: 25 - T-Junction with Vertex on Edge\n");
   
   poly1[0][0] = -1.5;
   poly1[0][1] = -1.0;
   poly1[0][2] = 0.0;
   poly1[1][0] = 1.5;
   poly1[1][1] = -1.0;
   poly1[1][2] = 0.0;
   poly1[2][0] = 1.5;
   poly1[2][1] = 1.5;
   poly1[2][2] = 0.0;
   poly1[3][0] = -1.5;
   poly1[3][1] = 1.5;
   poly1[3][2] = 0.0;

   poly2[0][0] = 0.0;
   poly2[0][1] = -1.0;
   poly2[0][2] = 0.0;
   poly2[1][0] = 1.0;
   poly2[1][1] = 0.5;
   poly2[1][2] = 0.0;
   poly2[2][0] = -1.0;
   poly2[2][1] = 0.5;
   poly2[2][2] = 0.0;

   PolyIntersect(poly1, poly2, polyi, 4, 3, &npoints);
   printf("Intersection (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyIntersect(poly2, poly1, polyi, 3, 4, &npoints);
   printf("Intersection (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly1, poly2, polyi, 4, 3, &npoints);
   printf("Merged polygon (poly1, poly2) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   PolyMerge(poly2, poly1, polyi, 3, 4, &npoints);
   printf("Merged polygon (poly2, poly1) points: %d\n", npoints);
   for(i=0; i<npoints; i++)
   {
     printf("Point: %le, %le, %le\n", polyi[i][0], polyi[i][1], polyi[i][2]);
   }

   return 0;
}
