#include "Building.h"
#include "geos_c.h"
#include <set>

using namespace std;

//FillBuilding
Building createBuilding(long id, OGRMultiPolygon* buildingGeom, map<int,Segment> segments){
	Building building;
	building.id = id;
	building.buildingGeom = buildingGeom;
	building.segments = segments;
	return building;
}

vector<BuildingCell> findSegmentsCoveringBuildingCells(OGRMultiPolygon* building, map<int,Segment> &segments){
	vector<BuildingCell> buildingCells;
	for(map<int,Segment>::iterator segment = segments.begin();segment!=segments.end();segment++){
		segment->second.cells.clear();
	}

	for(int CellNr=0 ; CellNr < building->getNumGeometries() ; CellNr++){
		BuildingCell cell;
		cell.area = ((OGRPolygon*)building->getGeometryRef(CellNr))->get_Area();
		for(map<int,Segment>::iterator segment = segments.begin();segment!=segments.end();segment++){
			if(segment->second.convexHull->Intersects(building->getGeometryRef(CellNr))){
				OGRGeometry* intersectGeom = (OGRPolygon*)segment->second.convexHull->Intersection(building->getGeometryRef(CellNr));
				//printf("area perc %lf\n", ((OGRPolygon*)intersectGeom)->get_Area()/segment->convexHull->get_Area());
				//printf("cov perc %lf\n", ((OGRPolygon*)intersectGeom)->get_Area()/((OGRPolygon*)building->getGeometryRef(CellNr))->get_Area());
				double intersectArea = ((OGRPolygon*)intersectGeom)->get_Area();
				double convexHullArea = segment->second.convexHull->get_Area();
				double cellArea = ((OGRPolygon*)building->getGeometryRef(CellNr))->get_Area();
				//test if intersection is larger then 5 percent of convexHull or 10 percent of cell
				if(intersectGeom != NULL && (intersectArea/convexHullArea > 0.05 || intersectArea/cellArea > 0.1) ){
					cell.segmentNumbers.push_back(segment->second.id);
					segment->second.cells.push_back(CellNr);
				}
			}
		}
		buildingCells.push_back(cell);
	}
	return buildingCells;
}

//Merge Cells
int GetGEOSNumberOfGeometries(GEOSGeom pGEOSGeom){
 	if(NULL == pGEOSGeom){
 		return 0;
 	}
 		 
 	int geometryType = GEOSGeomTypeId(pGEOSGeom);
 	if(geometryType == GEOS_POINT || geometryType == GEOS_LINESTRING || geometryType == GEOS_LINEARRING || geometryType == GEOS_POLYGON){
 		return 1;
 	}
 		 
 	//calling GEOSGetNumGeometries is save for multi types and collections also in geos2
 	return GEOSGetNumGeometries(pGEOSGeom);
} 

OGRMultiPolygon* mergePolygons(OGRMultiPolygon* Polygon,OGRMultiLineString* MergeLine){
	if(MergeLine == NULL){
		CPLError( CE_Failure, CPLE_ObjectNull, "Merge Line is null" );
		return NULL;
	}

	//1. Check intersects
	if(!Polygon->Intersects(MergeLine)){
		MergeLine->dumpReadable(NULL);
		CPLError( CE_Warning, CPLE_AppDefined, "Merge Line does not intersect polygon" );
		return NULL;
	}
	//2. Boundary
	OGRGeometry *pPolyBoundary = Polygon->Boundary();
	if(NULL == pPolyBoundary){
 		CPLError( CE_Failure, CPLE_ObjectNull, "Error get Boundary" );
 		return NULL;
	}
	//3. Union
	OGRGeometry *pNewPoly = pPolyBoundary->Difference(MergeLine);
	if(NULL == pNewPoly){
 		CPLError( CE_Failure, CPLE_ObjectNull, "Error union Polygon & Cut Line" );
 		return NULL;
	}
	GEOSGeometry *nodedGeometry = pNewPoly->exportToGEOS();

 	GEOSGeometry *cutEdges = GEOSPolygonizer_getCutEdges(&nodedGeometry, 1);
 	if (cutEdges){
 		if(GetGEOSNumberOfGeometries(cutEdges) > 0){
 		    GEOSGeom_destroy(cutEdges);
 		    GEOSGeom_destroy(nodedGeometry);
 		    CPLError( CE_Failure, CPLE_AppDefined, "Unknown error" );
 		    return NULL;
 		}
 		GEOSGeom_destroy( cutEdges );
 	}

	GEOSGeometry *polygons = GEOSPolygonize(&nodedGeometry, 1);
 	if (!polygons || GetGEOSNumberOfGeometries(polygons) == 0){
 		if (polygons)
 		GEOSGeom_destroy(polygons);
 		GEOSGeom_destroy(nodedGeometry);
 		CPLError( CE_Failure, CPLE_AppDefined, "Cut error" );
 		return NULL;
 	}

	GEOSGeom_destroy(nodedGeometry);
	OGRMultiPolygon* pOutputGeom = (OGRMultiPolygon*)OGRGeometryFactory::createGeometry(wkbMultiPolygon);
 		 
	for ( int i = 0; i < GetGEOSNumberOfGeometries(polygons); i++ ){
 		OGRGeometry* pGeom = OGRGeometryFactory::createFromGEOS((GEOSGeom)GEOSGetGeometryN(polygons, i));
 		pOutputGeom->addGeometry(pGeom);
	}
	GEOSGeom_destroy(polygons);

	return pOutputGeom;
}

bool TouchesByLine(OGRGeometry* cell1, OGRGeometry* cell2){
	OGRGeometry* intersection = cell1->Intersection(cell2);
	if(intersection!=NULL && (wkbFlatten(intersection->getGeometryType()) == wkbLineString
		|| wkbFlatten(intersection->getGeometryType()) == wkbMultiLineString) ){
		return true;
	}
	return false;
}

void cleanPolygons(OGRMultiPolygon* &polygon){
	for(int cellNr=0 ; cellNr < polygon->getNumGeometries() ; cellNr++){
		OGRPolygon* polyPtr = (OGRPolygon*)polygon->getGeometryRef(cellNr);
		if(polyPtr->getExteriorRing()->getNumPoints() > 5){
			OGRLinearRing ring;
			int pointCount = polyPtr->getExteriorRing()->getNumPoints();
			for(int coordNr=0;coordNr<pointCount-1;coordNr++){
				int pointNr1 = coordNr;
				int pointNr2 = coordNr+1;
				int pointNr3 = coordNr+2;
				
				if(coordNr==pointCount-2){
					pointNr1 = pointCount-2;
					pointNr2 = 0;
					pointNr3 = 1;
				}
				//printf("testing %d,%d,%d \n",pointNr1,pointNr2,pointNr3);
				
				OGRPoint point1;
				OGRPoint point2;
				OGRPoint point3;
				polyPtr->getExteriorRing()->getPoint(pointNr1,&point1);
				polyPtr->getExteriorRing()->getPoint(pointNr2,&point2);
				polyPtr->getExteriorRing()->getPoint(pointNr3,&point3);
				OGRLineString line;
				line.addPoint(&point1); 
				line.addPoint(&point3);
				//printf("distance is %lf\n",line->Distance(point2));
				if(line.Distance(&point2) > 0.01){
					ring.addPoint(&point2);
				}
			}
			ring.closeRings();
			//temp->dumpReadable(NULL);
			polyPtr->empty();
			polyPtr->addRing(&ring);
			//temp->dumpReadable(NULL);
		}
	}
}

bool cellsHaveSingleGabledRoof(vector<BuildingCell> buildingCells, map<int,Segment> segments, int mainCell, int subCell){
	//Test is on two cells covered by two segments of opposite direction where one cell contains a part of one segment
	bool merge=true;
	if(buildingCells.at(mainCell).segmentNumbers.size()+buildingCells.at(subCell).segmentNumbers.size() == 3){
		//mainCell could be gable
		if(buildingCells.at(mainCell).segmentNumbers.size() == 2){
			if(buildingCells.at(mainCell).segmentNumbers.at(0) == buildingCells.at(subCell).segmentNumbers.at(0)
					|| buildingCells.at(mainCell).segmentNumbers.at(1) == buildingCells.at(subCell).segmentNumbers.at(0) ){
				if(isGabledRoof(buildingCells.at(mainCell), segments)){
					return true;
				}
			}
		}else{ //subCell could be gable
			if(buildingCells.at(subCell).segmentNumbers.at(0) == buildingCells.at(mainCell).segmentNumbers.at(0)
					|| buildingCells.at(subCell).segmentNumbers.at(1) == buildingCells.at(mainCell).segmentNumbers.at(0) ){
				if(isGabledRoof(buildingCells.at(subCell), segments)){
					return true;
				}
			}
		}
	}
	else if(buildingCells.at(mainCell).segmentNumbers.size()+buildingCells.at(subCell).segmentNumbers.size() == 2){
		if(isShedRoof(buildingCells.at(mainCell),segments) && isShedRoof(buildingCells.at(subCell),segments)){
			//segments.at(segmentNr1).normal,segments.at(segmentNr2).normal)
			double angle = angleXYBetweenVectors(segments.at(buildingCells.at(mainCell).segmentNumbers.at(0)).normal, 
												segments.at(buildingCells.at(subCell).segmentNumbers.at(0)).normal);
			if(angle  > 175 && angle < 185){
					return true;
			}
		}
	}
	return false;
}

bool mergeDecompositionCells(OGRMultiPolygon* &building, vector<BuildingCell> buildingCells, map<int,Segment> segments){
	//x - make an array of segment numbers covering which building cell numbers
	//x - make a list of touching building cells that are covered by a single segment
	//x - merge building cells
	OGRMultiLineString mergeLines;

	for(int mainCellNr=0 ; mainCellNr < building->getNumGeometries() ; mainCellNr++){
		//Check if cell is covered by one segment
		for(int subCellNr=mainCellNr+1 ; subCellNr < building->getNumGeometries() ; subCellNr++){
			bool merge=false;
			//Check if other cell is covered by one segment and cells are covered by the same segment
			if(!buildingCells.at(mainCellNr).segmentNumbers.empty() && !buildingCells.at(subCellNr).segmentNumbers.empty()){
				//Merge cells containing identical segments
				if(buildingCells.at(mainCellNr).segmentNumbers == buildingCells.at(subCellNr).segmentNumbers){
					merge=true;
				}
				//Merge cells containing the same gabled roof
				else if(cellsHaveSingleGabledRoof(buildingCells,segments,mainCellNr,subCellNr)){
					merge=true;
				}
			}
			if(merge==true){
				//Check if cells are touching
				if(TouchesByLine(building->getGeometryRef(mainCellNr),building->getGeometryRef(subCellNr))){
					//They touch
					//printf("Cells %d and %d can be merged\n",mainCellNr,subCellNr);
					//add linestring of cells to list
					mergeLines.addGeometry((OGRLineString*)building->getGeometryRef(mainCellNr)->Boundary()->Intersection(building->getGeometryRef(subCellNr)->Boundary()));
				}
			}
		}
	}
	if(!mergeLines.IsEmpty()){
		building = mergePolygons(building,&mergeLines);
		cleanPolygons(building);
		//building->dumpReadable(NULL);
		return true;
	}
	return false;
}

void averageCellHeights(vector<BuildingCell> &building, double averageHeight){
	if(building.size() > 1){
		vector<int> flatList;
		//for(vector<BuildingCell>::iterator cell1 = building.begin();cell1!=building.end();cell1++){
		for(int cell=0;cell<(int)building.size();cell++){
			if(building.at(cell).roofType == Flat || building.at(cell).roofType == UnknownType){
				flatList.push_back(cell);
			}
		}
		//for(vector<int>::iterator cell1 = flatList.begin();cell1!=flatList.end();cell1++){
		for(int cell1 = 0;cell1<(int)flatList.size();cell1++){
			vector<int> sameHeightList;
			double height = building.at(flatList.at(cell1)).roofHeight * building.at(flatList.at(cell1)).area;
			double size = building.at(flatList.at(cell1)).area;
			
			//for(vector<int>::iterator cell2 = flatList.begin();cell2!=flatList.end();cell2++){
			for(int cell2 = cell1+1;cell2<(int)flatList.size();cell2++){
				if(abs(building.at(flatList.at(cell1)).roofHeight - building.at(flatList.at(cell2)).roofHeight) < averageHeight){
					sameHeightList.push_back(flatList.at(cell2));
					height += building.at(flatList.at(cell2)).roofHeight * building.at(flatList.at(cell2)).area;
					size += building.at(flatList.at(cell2)).area;
					flatList.erase(flatList.begin()+cell2);
					cell2--;
				}
			}
			if(sameHeightList.size() > 0){
				double newHeight = ceil((height/size)*100)/100;
				building.at(flatList.at(cell1)).roofHeight = newHeight;
				for(vector<int>::iterator cell2 = sameHeightList.begin();cell2!=sameHeightList.end();cell2++){
					building.at(*cell2).roofHeight = newHeight;
				}
			}
			flatList.erase(flatList.begin()+cell1);
		}
	}
}

void mergeAveragedCells(vector<BuildingCell> &building){
	if(building.size() > 1){
		for(int cell1 = 0;cell1<(int)building.size();cell1++){
			if(building.at(cell1).roofType == UnknownType){
				for(int cell2 = cell1+1;cell2<(int)building.size();cell2++){
					if(building.at(cell2).roofType == UnknownType
							&& TouchesByLine(building.at(cell1).roof, building.at(cell2).roof)
							&& abs(building.at(cell1).roofHeight - building.at(cell2).roofHeight) < 0.001){	
						set<int> segmentNumbers;
						//add SegmentNumbers
						segmentNumbers.insert(building.at(cell1).segmentNumbers.begin(),building.at(cell1).segmentNumbers.end());
						segmentNumbers.insert(building.at(cell2).segmentNumbers.begin(),building.at(cell2).segmentNumbers.end());
						building.at(cell1).segmentNumbers.clear();
						building.at(cell1).segmentNumbers.insert(building.at(cell1).segmentNumbers.begin(),segmentNumbers.begin(),segmentNumbers.end());
						//merge roof polygons
						OGRMultiPolygon* roof = (OGRMultiPolygon*)OGRGeometryFactory::createGeometry(wkbMultiPolygon);
						roof->addGeometry(((OGRPolygon*)building.at(cell1).roof)->Union((OGRPolygon*)building.at(cell2).roof));
						cleanPolygons(roof);
						building.at(cell1).roof = (OGRMultiPolygon*)roof->getGeometryRef(0);
						building.at(cell1).area = ((OGRPolygon*)building.at(cell1).roof)->get_Area();
						building.erase(building.begin()+cell2);
						cell2=cell1;
					}
				}
			}else if(building.at(cell1).roofType == Flat){
				for(int cell2 = cell1+1;cell2<(int)building.size();cell2++){
					if(building.at(cell2).roofType == Flat
							&& TouchesByLine(building.at(cell1).roof, building.at(cell2).roof)
							&& abs(building.at(cell1).roofHeight - building.at(cell2).roofHeight) < 0.001){	
						set<int> segmentNumbers;
						//add SegmentNumbers
						segmentNumbers.insert(building.at(cell1).segmentNumbers.begin(),building.at(cell1).segmentNumbers.end());
						segmentNumbers.insert(building.at(cell2).segmentNumbers.begin(),building.at(cell2).segmentNumbers.end());
						building.at(cell1).segmentNumbers.clear();
						building.at(cell1).segmentNumbers.insert(building.at(cell1).segmentNumbers.begin(),segmentNumbers.begin(),segmentNumbers.end());
						//merge roof polygons
						OGRMultiPolygon* roof = (OGRMultiPolygon*)OGRGeometryFactory::createGeometry(wkbMultiPolygon);
						roof->addGeometry(((OGRPolygon*)building.at(cell1).roof)->Union((OGRPolygon*)building.at(cell2).roof));
						cleanPolygons(roof);
						building.at(cell1).roof = (OGRMultiPolygon*)roof->getGeometryRef(0);
						building.at(cell1).area = ((OGRPolygon*)building.at(cell1).roof)->get_Area();
						building.erase(building.begin()+cell2);
						cell2=cell1;
					}
				}
			}
		}
	}
}

//ProcessBuilding
void processBuilding(Building &building, double averageHeight){
	building.buildingCells = findSegmentsCoveringBuildingCells(building.buildingGeom,building.segments);
	if( mergeDecompositionCells(building.buildingGeom,building.buildingCells,building.segments) ){
		building.buildingCells = findSegmentsCoveringBuildingCells(building.buildingGeom,building.segments);
	}

	//Calculate average height of segments
	calculateSegmentHeights(building.segments);

	//Check dec-cells one by one
	for(int buildingPart=0 ; buildingPart<building.buildingGeom->getNumGeometries() ; buildingPart++){
		processBuildingCell(building.buildingCells.at(buildingPart),building.segments,(OGRMultiPolygon*)building.buildingGeom->getGeometryRef(buildingPart));
	}

	averageCellHeights(building.buildingCells, averageHeight);

	mergeAveragedCells(building.buildingCells);
}

//CreateBuilding
bool linePerpendicularToDirection(vector<int> roofDirection, OGRPoint* point1, OGRPoint* point2, double dirEpsilon){
	dirEpsilon = dirEpsilon / (180/PI);

	double A = point2->getY()-point1->getY();
	double B = point2->getX()-point1->getX();
	double norm = sqrt((A*A)+(B*B));
	double normX = A / norm;
	double normY = -B / norm;
	double direction = atan2(normX,normY)-(0.5*PI);

	double segNorm = sqrt((double)(roofDirection[0]*roofDirection[0]) + (roofDirection[1]*roofDirection[1]) + (roofDirection[2]*roofDirection[2]));
	double segNormX = roofDirection[0]/segNorm;
	double segNormY = roofDirection[1]/segNorm;
	double segDirection = atan2(segNormX,segNormY);
	
	if(abs(direction-segDirection) < dirEpsilon //Direction is perpendicular
			|| abs((direction-PI)-segDirection) < dirEpsilon //Direction -180deg is perpendicular
			|| abs((direction+PI)-segDirection) < dirEpsilon //Direction +180deg is perpendicular
			|| abs((direction-(2*PI))-segDirection) < dirEpsilon //Direction -360deg is perpendicular
			|| abs((direction+(2*PI))-segDirection) < dirEpsilon){ //Direction +360deg is perpendicular
		return true;
	}
	return false;
}

bool lineParallelToDirection(vector<int> roofDirection, OGRPoint* point1, OGRPoint* point2, double dirEpsilon){
	dirEpsilon = dirEpsilon / (180/PI);

	double A = point2->getY()-point1->getY();
	double B = point2->getX()-point1->getX();
	double norm = sqrt((A*A)+(B*B));
	double normX = A / norm;
	double normY = -B / norm;
	double direction = atan2(normX,normY);

	double segNorm = sqrt((double)(roofDirection[0]*roofDirection[0]) + (roofDirection[1]*roofDirection[1]) + (roofDirection[2]*roofDirection[2]));
	double segNormX = roofDirection[0]/segNorm;
	double segNormY = roofDirection[1]/segNorm;
	double segDirection = atan2(segNormX,segNormY);
	
	if(abs(direction-segDirection) < dirEpsilon
		|| abs(direction-segDirection) < dirEpsilon){
		return true;
	}
	return false;
}

void addHeightToCell(BuildingCell &buildingCell){
	OGRMultiPolygon* newCellPolygon = (OGRMultiPolygon*)OGRGeometryFactory::createGeometry(wkbMultiPolygon);
	OGRLinearRing ring;
	if(!buildingCell.roof->IsEmpty() && wkbFlatten(buildingCell.roof->getGeometryType()) == wkbPolygon){
		OGRPolygon* roofPoly = (OGRPolygon*)buildingCell.roof;
		for(int pointNr=0;pointNr<roofPoly->getExteriorRing()->getNumPoints();pointNr++){
			OGRPoint point;
			roofPoly->getExteriorRing()->getPoint(pointNr,&point);
			point.setZ(buildingCell.roofHeight);
			ring.addPoint(&point);
		}
	}else{
		printf("addHeightToCell: roof is not a polygon");
	}
	OGRPolygon newPolygon;
	newPolygon.addRing(&ring);
	newCellPolygon->addGeometry(&newPolygon);
	buildingCell.roof = newCellPolygon;
}

bool addHeightToShedRoof(BuildingCell &buildingCell){
	if(!buildingCell.roof->IsEmpty() && wkbFlatten(buildingCell.roof->getGeometryType()) == wkbPolygon){
		OGRPolygon* roofPoly = (OGRPolygon*)buildingCell.roof;
		roofPoly->setCoordinateDimension(3);
		OGRLinearRing* roofRing = (OGRLinearRing*)roofPoly->getExteriorRing();
		OGRLineString ridgeLine;

		bool foundStartPoint = false;
		int startPoint = 0;
		for(startPoint;startPoint<roofRing->getNumPoints()-1;startPoint++){
			OGRPoint point1;
			OGRPoint point2;
			roofRing->getPoint(0+startPoint,&point1);
			roofRing->getPoint(1+startPoint,&point2);

			if(lineParallelToDirection(buildingCell.roofDirection,&point1,&point2,10)){
				foundStartPoint = true;
				ridgeLine.addPoint(&point1);
				ridgeLine.addPoint(&point2);
				break;
			}
		}

		if(!foundStartPoint){
			return false;
			//throw "ERROR: addHeightToShedRoof StartPoint not found";
		}

		for(int pointNr=0;pointNr<roofRing->getNumPoints();pointNr++){
			OGRPoint point;
			roofPoly->getExteriorRing()->getPoint(pointNr,&point);
			if(point.Intersects(&ridgeLine)){
				point.setZ(buildingCell.ridgeHeight); //set ridgeHeight
			}else point.setZ(buildingCell.eavesHeight); //set eavesHeight
			roofPoly->getExteriorRing()->setPoint(pointNr,&point);
		}
	}else{
		printf("addHeightToShedRoof: roof is not a Polygon\n");
	}
	OGRMultiPolygon* multiPoly = (OGRMultiPolygon*)OGRGeometryFactory::createGeometry(wkbMultiPolygon);
    multiPoly->addGeometryDirectly(buildingCell.roof);
	buildingCell.roof = multiPoly;
	return true;
}

void addHeightToGabledRoof(BuildingCell &buildingCell){
	if(!buildingCell.roof->IsEmpty() && wkbFlatten(buildingCell.roof->getGeometryType()) == wkbMultiPolygon
			&& buildingCell.roof->getNumGeometries() == 2){
		OGRLineString* ridgeLine = (OGRLineString*)buildingCell.roof->getGeometryRef(0)->Intersection(buildingCell.roof->getGeometryRef(1));
		
		for(int roofPart=0;roofPart<buildingCell.roof->getNumGeometries();roofPart++){
			OGRPolygon* roofPoly = (OGRPolygon*)buildingCell.roof->getGeometryRef(roofPart);
			roofPoly->setCoordinateDimension(3);

			for(int pointNr=0;pointNr<roofPoly->getExteriorRing()->getNumPoints();pointNr++){
				OGRPoint point;
				roofPoly->getExteriorRing()->getPoint(pointNr,&point);
				if(point.Intersects(ridgeLine)){
					point.setZ(buildingCell.ridgeHeight); //set ridgeHeight
				}else point.setZ(buildingCell.eavesHeight); //set eavesHeight
				roofPoly->getExteriorRing()->setPoint(pointNr,&point);
			}
		}
	}else{
		printf("addHeightToGabledRoof: roof is not a MultiPolygon with 2 rings\n");
	}
}

bool calculateGabledRoof(BuildingCell &cell){
	OGRLinearRing* roofRing = (OGRLinearRing*) ((OGRPolygon*)cell.roof)->getExteriorRing();
	//roofRing->dumpReadable(NULL);
		
	if(roofRing->getNumPoints() > 5){
		OGRGeometry* test = roofRing->Simplify(0.01);
		if(roofRing->getNumPoints() > 5){
			printf("calculateGabledRoof: ring has >5 points\n");
		}
	}
	
	bool foundStartPoint = false;
	int startPoint = 0;
	for(startPoint;startPoint<roofRing->getNumPoints()-1;startPoint++){
		OGRPoint point1;
		OGRPoint point2;
		roofRing->getPoint(0+startPoint,&point1);
		roofRing->getPoint(1+startPoint,&point2);

		if(linePerpendicularToDirection(cell.roofDirection,&point1,&point2,10)){
			foundStartPoint = true;
			break;
		}
	}
	if(!foundStartPoint){
		printf("ERROR: CalculateGabledRoof StartPoint not found");
		return false;
		//throw "ERROR: CalculateGabledRoof StartPoint not found";
	}
	if(startPoint > roofRing->getNumPoints()-3){
		printf("ERROR: CalculateGabledRoof StartPoint too large");
		return false;
		//throw "ERROR: CalculateGabledRoof StartPoint too large";
	}
	
	double split1x,split1y,split2x,split2y = 0;

	split1x = (roofRing->getX(0+startPoint) + roofRing->getX(1+startPoint)) /2;
	split1y = (roofRing->getY(0+startPoint) + roofRing->getY(1+startPoint)) /2;
	split2x = (roofRing->getX(2+startPoint) + roofRing->getX(3+startPoint)) /2;
	split2y = (roofRing->getY(2+startPoint) + roofRing->getY(3+startPoint)) /2;

	OGRLinearRing ring1;
	OGRLinearRing ring2;

	ring1.addPoint(roofRing->getX(1+startPoint),roofRing->getY(1+startPoint));
	ring1.addPoint(roofRing->getX(2+startPoint),roofRing->getY(2+startPoint));
	ring1.addPoint(split2x,split2y);
	ring1.addPoint(split1x,split1y);
	ring1.addPoint(roofRing->getX(1+startPoint),roofRing->getY(1+startPoint));
	//ring1->dumpReadable(NULL);

	ring2.addPoint(roofRing->getX(0+startPoint),roofRing->getY(0+startPoint));
	ring2.addPoint(split1x,split1y);
	ring2.addPoint(split2x,split2y);
	ring2.addPoint(roofRing->getX(3+startPoint),roofRing->getY(3+startPoint));
	ring2.addPoint(roofRing->getX(0+startPoint),roofRing->getY(0+startPoint));
	//ring2->dumpReadable(NULL);

	OGRMultiPolygon* splittedRoof = (OGRMultiPolygon*)OGRGeometryFactory::createGeometry(wkbMultiPolygon);
	OGRPolygon newCell;
	newCell.addRing(&ring1);
	splittedRoof->addGeometry(&newCell);
	newCell.empty();
	newCell.addRing(&ring2);
	splittedRoof->addGeometry(&newCell);

	//splittedRoof->dumpReadable(NULL);
	
	cell.roof = splittedRoof;
	return true;
}

OGRMultiPolygon* calculateWalls(OGRMultiPolygon* roof){
	//Calculate walls
	OGRMultiPolygon* walls = (OGRMultiPolygon*)OGRGeometryFactory::createGeometry(wkbMultiPolygon);
	for(int geomNr=0;geomNr<roof->getNumGeometries();geomNr++){
		OGRLinearRing* roofRing = ((OGRPolygon*)roof->getGeometryRef(geomNr))->getExteriorRing();

		for(int pointNr=0;pointNr < roofRing->getNumPoints()-1;pointNr++){
			OGRPoint point1;
			OGRPoint point2;
			OGRLinearRing ring;
			OGRPolygon wall;
			roofRing->getPoint(pointNr,&point1);
			roofRing->getPoint(pointNr+1,&point2);
			ring.addPoint(&point1);
			ring.addPoint(&point2);
			ring.addPoint(point2.getX(),point2.getY(),0);
			ring.addPoint(point1.getX(),point1.getY(),0);
			ring.addPoint(&point1);
			wall.addRing(&ring);
			walls->addGeometry(&wall);
		}
	}
	return walls;
}

//OGRMultiPolygon* calculateWalls(OGRMultiPolygon* roof){
//	//Calculate walls
//	OGRMultiPolygon* walls = new OGRMultiPolygon;
//	//if(roof->getNumGeometries() == 1){
//	for(int geomNr=0;geomNr<roof->getNumGeometries();geomNr++){
//		OGRLinearRing* roofRing = ((OGRPolygon*)roof->getGeometryRef(geomNr))->getExteriorRing();
//
//		for(int pointNr=0;pointNr < roofRing->getNumPoints()-1;pointNr++){
//			OGRPoint point1;
//			OGRPoint point2;
//			OGRLinearRing ring;
//			OGRPolygon wall;
//			roofRing->getPoint(pointNr,&point1);
//			roofRing->getPoint(pointNr+1,&point2);
//			ring.addPoint(&point1);
//			ring.addPoint(&point2);
//			ring.addPoint(point2.getX(),point2.getY(),0);
//			ring.addPoint(point1.getX(),point1.getY(),0);
//			ring.addPoint(&point1);
//			wall.addRing(&ring);
//			walls->addGeometry(&wall);
//		}
//	}
//	return walls;
//}

OGRPoint* findPoint(OGRPoint* inPoint, OGRMultiPolygon* roof){
	for(int geom=0;geom<roof->getNumGeometries();geom++){
		OGRLinearRing* roofRing = ((OGRPolygon*)roof->getGeometryRef(geom))->getExteriorRing();
		for(int pointnr=0;pointnr < roofRing->getNumPoints();pointnr++){
			OGRPoint* point = (OGRPoint*)OGRGeometryFactory::createGeometry(wkbPoint);
			roofRing->getPoint(pointnr,point);
			if(point->Intersects(inPoint))
				return point;
		}
	}
	return NULL;
}

OGRMultiPolygon* calculateGabledWalls(OGRMultiPolygon* roof){
	//Calculate walls
	OGRMultiPolygon* walls = (OGRMultiPolygon*)OGRGeometryFactory::createGeometry(wkbMultiPolygon);
	OGRPolygon* flattenedRoof = (OGRPolygon*)roof->UnionCascaded();
	OGRLinearRing* roofRing = flattenedRoof->getExteriorRing();

	for(int pointNr=0;pointNr < roofRing->getNumPoints()-1;pointNr++){
		OGRPoint* point1 = (OGRPoint*)OGRGeometryFactory::createGeometry(wkbPoint);
		OGRPoint* point2 = (OGRPoint*)OGRGeometryFactory::createGeometry(wkbPoint);
		OGRLinearRing ring;
		OGRPolygon wall;
		roofRing->getPoint(pointNr,point1);
		roofRing->getPoint(pointNr+1,point2);
	
		if(point1->Intersects(roof->getGeometryRef(0)) && point1->Intersects(roof->getGeometryRef(1))){
			printf("calculateGabledWalls point1 is on ridge\n");
		}
		if(point2->Intersects(roof->getGeometryRef(0)) && point2->Intersects(roof->getGeometryRef(1))){ //ridge point
			OGRPoint* point3 = (OGRPoint*)OGRGeometryFactory::createGeometry(wkbPoint);
			if(pointNr+2>roofRing->getNumPoints()-1) //if ridge point is endpoint take point 1 (NOT 0!)
				roofRing->getPoint(1,point3);
			else
				roofRing->getPoint(pointNr+2,point3);

			point1 = findPoint(point1,roof);
			point2 = findPoint(point2,roof);
			point3 = findPoint(point3,roof);
			ring.addPoint(point1);
			ring.addPoint(point2);
			ring.addPoint(point3);
			ring.addPoint(point3->getX(),point3->getY(),0);
			ring.addPoint(point1->getX(),point1->getY(),0);
			ring.addPoint(point1);
			wall.addRing(&ring);
			pointNr++;
		}else{
			point1 = findPoint(point1,roof);
			point2 = findPoint(point2,roof);
			ring.addPoint(point1);
			ring.addPoint(point2);
			ring.addPoint(point2->getX(),point2->getY(),0);
			ring.addPoint(point1->getX(),point1->getY(),0);
			ring.addPoint(point1);
			wall.addRing(&ring);
		}
		walls->addGeometry(&wall);
	}
	return walls;
}

void correctRoofOrientation(OGRMultiPolygon* roof){
	if(wkbFlatten(roof->getGeometryType()) == wkbPolygon){
		((OGRPolygon*)roof)->getExteriorRing()->reverseWindingOrder();
	}else{
		for(int roofPart=0;roofPart<roof->getNumGeometries();roofPart++){
			((OGRPolygon*)roof->getGeometryRef(roofPart))->getExteriorRing()->reverseWindingOrder();
		}
	}
}

void createBuildingGeometryLOD1(vector<BuildingCell> &building){
	for(vector<BuildingCell>::iterator part = building.begin();part!=building.end();++part){
		if(!((OGRPolygon*)part->roof)->getExteriorRing()->isClockwise()){
			//printf("ERROR: unknown roof is oriented CCW\n");
			correctRoofOrientation(part->roof);
		}
		part->floor = (OGRPolygon*)part->roof->clone();
		addHeightToCell(*part);
		part->walls = calculateWalls(part->roof);
		correctRoofOrientation(part->roof);
		if(((OGRPolygon*)part->roof->getGeometryRef(0))->getExteriorRing()->isClockwise()){
			printf("ERROR: Roof is oriented CW after correction\n");
		}
	}
}

void createBuildingGeometryLOD2(vector<BuildingCell> &building){
	for(vector<BuildingCell>::iterator part = building.begin();part!=building.end();++part){
		if(part->roofType==Flat){
			if(!((OGRPolygon*)part->roof)->getExteriorRing()->isClockwise()){
				//printf("ERROR: flat roof is oriented CCW\n");
				correctRoofOrientation(part->roof);
			}
			part->floor = (OGRPolygon*)part->roof->clone();
			addHeightToCell(*part);
			part->walls = calculateWalls(part->roof);
			correctRoofOrientation(part->roof);
			if(((OGRPolygon*)part->roof->getGeometryRef(0))->getExteriorRing()->isClockwise()){
				printf("ERROR: flat roof is oriented CW after correction\n");
			}
		}else if(part->roofType==Shed){
			if(!((OGRPolygon*)part->roof)->getExteriorRing()->isClockwise()){
				//printf("ERROR: shed roof is oriented CCW\n");
				correctRoofOrientation(part->roof);
			}
			part->floor = (OGRPolygon*)part->roof->clone();
			if(addHeightToShedRoof(*part)){
				part->walls = calculateWalls(part->roof);
				correctRoofOrientation(part->roof);
				if(((OGRPolygon*)part->roof->getGeometryRef(0))->getExteriorRing()->isClockwise()){
					printf("ERROR: shed roof is oriented CW after correction\n");
				}
			}else{
				part->roofType=UnknownType;
				part--;
			}
		}else if(part->roofType==Gabled){
			if(!((OGRPolygon*)part->roof)->getExteriorRing()->isClockwise()){
				//printf("ERROR: gabled roof is oriented CCW\n");
				correctRoofOrientation(part->roof);
			}
			part->floor = (OGRPolygon*)part->roof->clone();
			if(calculateGabledRoof(*part)){
				addHeightToGabledRoof(*part);
				if(!((OGRPolygon*)part->roof->getGeometryRef(0))->getExteriorRing()->isClockwise()){
					correctRoofOrientation((OGRMultiPolygon*)part->roof->getGeometryRef(0));
				}
				if(!((OGRPolygon*)part->roof->getGeometryRef(1))->getExteriorRing()->isClockwise()){
					correctRoofOrientation((OGRMultiPolygon*)part->roof->getGeometryRef(1));
				}
				if(!((OGRPolygon*)part->roof->getGeometryRef(0))->getExteriorRing()->isClockwise()
					&& !((OGRPolygon*)part->roof->getGeometryRef(1))->getExteriorRing()->isClockwise()){
					printf("ERROR: gabled roof is oriented CCW after correction\n");
				}
				//part->walls = calculateWalls(part->roof);
				part->walls = calculateGabledWalls(part->roof);
				correctRoofOrientation(part->roof);
				if(((OGRPolygon*)part->roof->getGeometryRef(0))->getExteriorRing()->isClockwise()
					|| ((OGRPolygon*)part->roof->getGeometryRef(1))->getExteriorRing()->isClockwise()){
					printf("ERROR: gabled roof is oriented CW after correction\n");
				}
			}else{
				part->roofType=UnknownType;
				part--;
			}
		}else if(part->roofType==NoSegments){
			//Do nothing
		}else{
			if(!((OGRPolygon*)part->roof)->getExteriorRing()->isClockwise()){
				//printf("ERROR: unknown roof is oriented CCW\n");
				correctRoofOrientation(part->roof);
			}
			part->floor = (OGRPolygon*)part->roof->clone();
			addHeightToCell(*part);
			part->walls = calculateWalls(part->roof);
			correctRoofOrientation(part->roof);
			if(((OGRPolygon*)part->roof->getGeometryRef(0))->getExteriorRing()->isClockwise()){
				printf("ERROR: unknown roof is oriented CW after correction\n");
			}
		}
	}
}