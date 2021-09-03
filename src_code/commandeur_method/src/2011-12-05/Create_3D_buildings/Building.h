#ifndef _BUILDINGCELL_H_INCLUDED
#define _BUILDINGCELL_H_INCLUDED
#include "BuildingCell.h"
#endif

#ifndef _SEGMENTS_H_INCLUDED
#define _SEGMENTS_H_INCLUDED
#include "Segments.h"
#endif

struct Building{
	long id;
	std::vector<long> adjacentBuildings;
	OGRMultiPolygon* buildingGeom;
	std::vector<BuildingCell> buildingCells;
	std::map<int,Segment> segments;
};

//Building
Building createBuilding(long id, OGRMultiPolygon* buildingGeom, std::map<int,Segment> segments);
void processBuilding(Building &building, double averageHeight);
std::vector<BuildingCell> processBuilding(OGRMultiPolygon* buildingGeom, std::map<int,Segment> segments);
void createBuildingGeometryLOD1(std::vector<BuildingCell> &building);
void createBuildingGeometryLOD2(std::vector<BuildingCell> &building);
std::vector<BuildingCell> findSegmentsCoveringBuildingCells(OGRMultiPolygon* building, std::map<int,Segment> &segments);
