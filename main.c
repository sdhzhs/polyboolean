#include <stdio.h>
#include <math.h>
#include "convex_poly_boolean.h"

int main (int argc, char *argv[])
{
   int i, j, k, npoints;
   real poly1[6][3], poly2[4][3], polyi[16][3];

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

   //intersection and union of two congruent triangles, horizontal movement, opposite points order
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

   return 0;
}
