#include <cstdlib>
//#include <iostream>
#include "InlineArguments.h"
#include "SegmentationParameters.h"

using namespace std;

void PrintUsage()
{
  printf("Usage: segmentation_complete -i <input text file>\n");
  printf("                             -o <output text file>\n");
  
  printf("\n         General optional parameters\n");
  printf("                   -rd <remove double points>\n");
  printf("                   -h  <number of header records, def: 0>\n");

  printf("\n         Parameters for specifying columns of attributes\n");
  printf("                   -x  <column number of X-coordinate, def: 1>\n");
  printf("                   -y  <column number of Y-coordinate, def: 2>\n");
  printf("                   -z  <column number of Z-coordinate, def: 3>\n");
  printf("Type \"segmentation_complete -parameters\" to get the segmentation parameters\n");
}

void PrintParameters()
{
  printf("Segmentation parameters of the segmentation_complete programme\n");
  printf("        [-tin OR -octree OR -knn (default)] : neighbourhood type\n");
  printf("        [-ocbinmax <size> (def: 100)        : maximum octree bin size\n");
  printf("        [-ocbinoverlap <size> (def: 1.0)    : maximum octree bin overlap\n");
  printf("        [-knn <number> (def: 20)            : number of nearest neighbours\n");
  printf("        [-dim <2 or 3> (def: 3)             : neighbour metric dimension\n");
  printf("        [-minsegsize <size> (def: 10)       : minimum segment size\n\n");
  printf("        [-seednbh <0: direct (def), 1: dist>: seed neighbourhood definition\n");
  printf("        [-seedradius <number> (def: 1.0)    : seed neighbourhood radius\n");
  printf("        [-maxslope <degrees> (def: 90.0)    : maximum plane slope\n");
  printf("        [-binsizeslope <degrees> (def: 3.0) : Hough bin size of slope\n");
  printf("        [-binsizedist <number> (def: 0.2)   : Hough bin size of distance\n");
  printf("        [-minsizeseed <number> (def: 10)    : minimum seed size (#pts)\n");
  printf("        [-maxdistseed <number> (def: 0.2)   : maximum distance of point to seed\n\n");
  printf("        [-smooth] OR [-plane (def)]         : surface type to be grown\n");
  printf("        [-grownbh <0: direct (def), 1: dist>: growing neighbourhood definition\n");
  printf("        [-growradius <number> (def: 0.3)    : growing search radius\n");
  printf("        [-maxdistgrow <number> (def: 0.3)   : maximum distance of point to surface\n");
  printf("        [-mindistrecompute <num> (def: 0.15): minimum distance to recompute surface\n");
  printf("        [-compete]                          : let planes compete for points on edge\n");
}

int main(int argc, char *argv[])
{
  InlineArguments *args = new InlineArguments(argc, argv);
  SegmentationParameters par;
  
  void segmentation_complete(char *, char *, int, int, int, int, int,
                    const SegmentationParameters &);
  
  // Check on required input files
  if (args->Contains("-usage")) {
    PrintUsage();
    exit(0);
  }
  if (args->Contains("-parameters")) {
    PrintParameters();
    exit(0);
  }                  
  if (!args->Contains("-i")) {
    printf("Error: -i is a required argument.\n");
    PrintUsage();
    exit(0);
  }
  if (!args->Contains("-o")) {
    printf("Error: -o is a required argument.\n");
    PrintUsage();
    exit(0);
  }
  
  // Get all individually set segmentation parameters
  // Neighbourhood storage model
  if (args->Contains("-tin")) par.NeighbourhoodStorageModel() = 0;
  else if (args->Contains("-octree")) {
    par.NeighbourhoodStorageModel() = 1;
    par.OctreeBinMaxNumberOfPoints() = args->Integer("-ocbinmax", 100);
    par.OctreeBinOverlap() = args->Double("-ocbinoverlap", 1.0);
  }
  else if (!args->Contains("-par")) {
    par.NeighbourhoodStorageModel() = 2; // knn, default
    par.NumberOfNeighbours() = args->Integer("-knn", 20);
  }
  if (args->Contains("-dim"))
    par.DistanceMetricDimension() = args->Integer("-dim", 3);
  
  // Minimum component size
  if (args->Contains("-minsegsize"))
    par.MinNumberOfPointsComponent() = args->Integer("-minsegsize", 10);

  // Seed selection parameters
  if (args->Contains("-seednbh"))
    par.SeedNeighbourhoodDefinition() = args->Integer("-seednbh", 0);
  if (args->Contains("-seedradius"))
    par.SeedNeighbourhoodRadius() = args->Double("-seedradius", 1.0);
  if (args->Contains("-maxslope"))
    par.MaxSlopeAngle() = args->Double("-maxslope", 90.0) * atan(1.0) / 45.0;
  if (args->Contains("-binsizeslope"))
    par.BinSizeSlopeAngle() = args->Double("-binsizeslope", 3.0) * atan(1.0) / 45.0;
  if (args->Contains("-binsizedist"))
    par.BinSizeDistance() = args->Double("-binsizedist", 0.2);
  if (args->Contains("-minsizeseed"))
    par.MinNumberOfPointsSeed() = args->Integer("-minsizeseed", 10);
  if (args->Contains("-maxdistseed"))
    par.MaxDistanceSeedPlane() = args->Double("-maxdistseed", 0.2);
  
  // Surface growing parameters
  if (args->Contains("-smooth")) par.SurfaceModel() = 1;
  else if (!args->Contains("-par")) par.SurfaceModel() = 0; // Plane
  if (args->Contains("-grownbh"))
    par.GrowingNeighbourhoodDefinition() = args->Integer("-grownbh", 0);
  if (args->Contains("-growradius"))
    par.GrowingRadius() = args->Double("-growradius", 0.3);
  if (args->Contains("-maxdistgrow"))
    par.MaxDistanceSurface() = args->Double("-maxdistgrow", 0.3);
  if (args->Contains("-mindistrecompute"))
    par.MinDistanceRecompute() = args->Double("-mindistrecompute", 0.15);
  if (args->Contains("-compete")) par.SurfacesCompete() = true;
  
  // Segment the data
  segmentation_complete(args->String("-i"),args->String("-o"),
                args->Integer("-x", 1),args->Integer("-y", 2),args->Integer("-z", 3),
                args->Contains("-rd"),args->Integer("-h", 0),par); 
                
  return EXIT_SUCCESS;
}
