#include <stdio.h>
#include "convex_poly_boolean.h"

int main (int argc, char *argv[])
{
   int i, j, npoints;
   double poly1[6][3], poly2[4][3], polyi[16][3];

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

   for(j=0; j<4; j++)
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

   return 0;
}
