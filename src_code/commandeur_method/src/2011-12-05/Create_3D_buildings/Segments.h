#include <vector>
#include <map>
#include "ogrsf_frmts.h"

struct Segment{
	int id;
	//int Normal_X, Normal_Y, Normal_Z;
	std::vector<int> normal;
	std::vector<int> cells;
	OGRMultiPoint* points;
	OGRPolygon* convexHull;
	double averageHeight;
	virtual ~Segment(){};
};

//Segments
std::map<int,Segment> findSegmentsInBuilding(OGRLayer* pointsLayer, OGRGeometry* buildingGeometry);
double calculateEavesHeight(std::map<int,Segment> segments,std::vector<int> segmentNumbers);
double calculateRidgeHeight(std::map<int,Segment> segments,std::vector<int> segmentNumbers);
void calculateSegmentHeights(std::map<int,Segment> &segments);
double calculateWeightedAverageCellHeight(std::map<int,Segment> &segments, std::vector<int> segmentNumbers, OGRMultiPolygon* roof);
double calculateHeightMultipleSegments(std::map<int,Segment> &segments, std::vector<int> segmentNumbers);