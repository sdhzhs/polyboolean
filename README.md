# polyboolean
Background

The C codes in this repository are extracted and extended from one Computer Aided Engineering (CAE) project.  The main goal is to get intersections or unions of two 3D convex polygons.  The codes emphasize the implementation of a general algorithm for Boolean operations of polygons, in which the intersection and inner points are used to construct new polygon for each Boolean operation.  This algorithm is similar to the widely known Sutherland-Hodgman algorithm for convex polygon clipping, and implements general idea for both intersection and merge operations.

Data structure

The data structure for polygons is just array in built-in floating number which represents ordered points set in 3D space.  The order of points set in the first (poly1) and second (poly2) polygon should be identical (both clockwise or both counterclockwise if viewing along a certain direction).  Obviously the polygons should be planar (or approximately planar based on geometric tolerance) for valid Boolean operations.  The reason why the polygons are 3D is that the application scenarios in CAE are always three dimensional and some approximate cases should be treated.  Because the focus of functions in the codes is on Boolean operations of two polygons, no extra tree-like data structure is introduced to improve the efficiency of geometric operations, which is commonly used in Boolean operations of polygons (graphs) in large amount.  The application in CAE shows the codes are efficient if the number of points and Boolean operations are not larger than 1M.

Algorithm

Even though some corner scenarios for predication of inner and intersection relationships are considered by the codes, the corner cases may still exist in application and these bugs of corner type will be fixed as it is discovered.  Note that the intersection of two convex polygons is still convex, but the union of two convex polygons may not be convex.  Furthermore, the number of result polygons may be more than one (topological change) in other Boolean operations than intersection and union (e.g., difference set of two polygons).  Despite this algorithm is believed to be general for all Boolean operations, the implementation of other Boolean operations is non-trivial and may be extended in the future.

Usage

There is a Makefile in this repository, just type

make

in a Unix-like OS environment in which gcc has been deployed, then a shared library named "libpolyboolean.so" and an exec named "polyop" will be generated. You can also use the script in file cclib2exe to build this programme. 

The codes in source file main.c represent some examples which call Boolean operation functions in source file convex_poly_boolean.c.  The examples include
1.  intersection and union sets of two moving triangles.
2.  intersection and union sets between one static hexagon and one moving square.

There is a golden file output_golden which records all the ouputs of examples, after running exec file like

bin/polyop > output

the accuracy of the codes can be monitored using the command

diff output output_golden

Note that enough memory should be allocated for result polygons after Boolean operations. There are comments in source file convex_poly_boolean.c, which describe each Boolean operation and helper function. The functions in geometry3d_basic.c belong to low-level geometric operations and can be used for other targets of computational geometry.
