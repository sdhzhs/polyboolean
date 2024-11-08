Background

The C codes in this repository are extracted and extended from one Computer Aided Engineering (CAE) project.  The main goal is to get intersections or unions of two 3D convex polygons.  The codes emphasize the implementation of a general algorithm for boolean operations of polygons, in which the intersection and inner points are used to construct new polygon for each boolean operation.  This algorithm is similar to the widely known Sutherland-Hodgman algorithm for convex polygon clipping, and implement general idea for both intersection and merge operations.

Data structure

The data structure for polygons is just array in double precision floating number which represents ordered points set in 3D space.  The order of ponits set in the first (poly1) and second (poly2) polygon should be identical (both clockwise or both counterclockwise).  Obviously the polygons should be planar (or approximately planar based on geometric tolerence) for valid boolean operations.  The reason why the polygons are 3D is the application scenarios in CAE is always three dimensional and some approximate cases should be treated.  Because the focus of functions in the codes is on boolean operations of two polygons, no extra tree-like data structure is introduced to improve the efficiency of geometric calculation, which is commonly used in boolean operations of polygons (graphs) in large amount.  The application in CAE shows the codes are efficient if the number of points and boolean operations are not larger than 1M.

Algorithm

Even though some corner scenarios for predication of inner and intersection relationships are considered by the codes, the corner cases may still exist in application and these bugs of corner type will be fixed as it is discovered.  Note that the intersection of two convex polygons is still convex, but the union of two convex polygons may not be convex.  Furthermore, the number of result polygons may be more than one (topology change) in other boolean operations than intersection and union (e.g., difference set of two polygons).  Despite this algorithm is believed to be general for all boolean operations, the implementation of other boolean operations is non-trivial and may be extended in the future.

Usage

There is a Makefile in this repository, just type

make

in a Unix-like OS environment in which gcc has been deployed, then an exec named "polyop" will be generated.  You can also use the script in file cclib2exe to build this programme.

The codes in source file main.c represent some examples which call geometric functions in source file convex_poly_boolean.c.  The examples include
1.  intersection and union sets of two moving triangles.
2.  intersection and union sets between one static hexagon and one moving square.

Note that enough memory should be allocated for result polygon after boolean operations. There are comments in source file convex_poly_boolean.c, which describe each geometric function. The functions in section of basic geometry calculation are low-level and can be used for other target of computational geometry.
