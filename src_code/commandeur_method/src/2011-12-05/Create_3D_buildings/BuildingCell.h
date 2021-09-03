#ifndef _SEGMENTS_H_INCLUDED
#define _SEGMENTS_H_INCLUDED
#include "Segments.h"
#endif

const double PI = 4*atan(1.0);

typedef enum RoofType{
	UnknownType = 0,
	Flat = 1,
	Shed = 2,
	Gabled = 3,
	NoSegments = 4
} RoofType;

struct BuildingCell{
	std::vector<int> segmentNumbers;
	RoofType roofType;
	std::vector<int> roofDirection;
	double roofAngle;
	double area;
	double roofHeight;
	double eavesHeight;
	double ridgeHeight;
	OGRMultiPolygon* roof;
	OGRMultiPolygon* walls;
	OGRPolygon* floor;
};

void processBuildingCell(BuildingCell &buildingCell, std::map<int,Segment> segments, OGRMultiPolygon* cellPolygon);
bool isShedRoof(BuildingCell cell, std::map<int,Segment> segments);
bool isGabledRoof(BuildingCell cell, std::map<int,Segment> segments);
double angleXYBetweenVectors(std::vector<int> vector1, std::vector<int> vector2);