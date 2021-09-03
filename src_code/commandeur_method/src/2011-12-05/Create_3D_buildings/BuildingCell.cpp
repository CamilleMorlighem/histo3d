#include "BuildingCell.h"

using namespace std;

double roundDouble(double in){
	if(in-floor(in) >= 0.5){
		return ceil(in);
	}else return floor(in);
}

//Find RoofType
double angleWithHorizontalVector(vector<int> vector1){
	double lenA = sqrt((double)vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

	return abs(acos(vector1[2] / (lenA)) * (180/PI));
}

//double angleWithHorizontalVector(vector<int> vector1){
//	vector<int> vector2;
//	vector2.push_back(0);
//	vector2.push_back(0);
//	vector2.push_back(1);
//
//	double lenA = sqrt((double)vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);
//	double lenB = sqrt((double)vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);
//
//	double dotAB = vector1[0]*vector2[0] + vector1[1]*vector2[1] + vector1[2]*vector2[2];
//
//	return abs(acos(dotAB / (lenA*lenB)) * (180/PI));
//}

double angleXYBetweenVectors(vector<int> vector1, vector<int> vector2){
	double lenA = sqrt((double)vector1[0]*vector1[0] + vector1[1]*vector1[1]);
	double lenB = sqrt((double)vector2[0]*vector2[0] + vector2[1]*vector2[1]);

	double dotAB = vector1[0]*vector2[0] + vector1[1]*vector2[1];
	
	if( abs(abs(dotAB)-abs(lenA*lenB)) < 0.001){
		return 180;
	}

	//double temp = abs(acos(dotAB / (lenA*lenB)) * (180/PI));
	return abs(acos(dotAB / (lenA*lenB)) * (180/PI));
}

vector<int> getRoofDirection(map<int,Segment> segments, vector<int> segNrs){
	if(segments.at(segNrs.at(0)).normal.at(1) > 0 ) //Check if normalY is positive, then vector points up
		return segments.at(segNrs.at(0)).normal;
	else
		return segments.at(segNrs.at(1)).normal;
}

bool isFlatRoof(BuildingCell cell, map<int,Segment> segments){
	if(cell.segmentNumbers.size()==1
		&& angleWithHorizontalVector(segments.at(cell.segmentNumbers.at(0)).normal) < 10 ){
		return true;
	}
	return false;
}

bool isShedRoof(BuildingCell cell, map<int,Segment> segments){
	if(cell.segmentNumbers.size()==1
		&& angleWithHorizontalVector(segments.at(cell.segmentNumbers.at(0)).normal) > 10 ){
		return true;
	}
	return false;
}

bool isGabledRoof(BuildingCell cell, map<int,Segment> segments){
	if(cell.segmentNumbers.size()==2){
		int segmentNr1 = cell.segmentNumbers.at(0);
		int segmentNr2 = cell.segmentNumbers.at(1);
		double segmentAngle1 = angleWithHorizontalVector(segments.at(segmentNr1).normal);
		double segmentAngle2 = angleWithHorizontalVector(segments.at(segmentNr2).normal);
		double segmentsAngle = angleXYBetweenVectors(segments.at(segmentNr1).normal,segments.at(segmentNr2).normal);
		if( segmentAngle1 > 5
			&& segmentAngle2 > 5
			&& segmentsAngle > 175
			&& segmentsAngle < 185){
			return true;
		}
	}
	return false;
}

bool isRectangularCell(OGRMultiPolygon* cellPolygon){
	//before the number of points is checked, known is the poly has 5 points
	OGRLineString line;
	OGRPoint point;
	vector<double> lengths;

	for(int lineNr=0;lineNr<((OGRPolygon*)cellPolygon)->getExteriorRing()->getNumPoints()-1;lineNr++){
		((OGRPolygon*)cellPolygon)->getExteriorRing()->getPoint(lineNr,&point);
		line.addPoint(&point);
		((OGRPolygon*)cellPolygon)->getExteriorRing()->getPoint(lineNr+1,&point);
		line.addPoint(&point);
		lengths.push_back(line.get_Length());
		line.empty();
	}

	//angles are not checked, angle 0+2 should be ~180 deg.
	if((abs(lengths.at(0)-lengths.at(2)) / lengths.at(2)) < 0.1
		&& (abs(lengths.at(1)-lengths.at(3)) / lengths.at(3)) < 0.1){
		return true;
	}
	
	
	return false;
}

void processBuildingCell(BuildingCell &buildingCell, map<int,Segment> segments, OGRMultiPolygon* cellPolygon){
	buildingCell.roof = cellPolygon;
	if(buildingCell.segmentNumbers.size()==0){ //If cell has no segments, write attribute RoofType "NoSegments"
		buildingCell.roofType = NoSegments;
		buildingCell.roofAngle = 0;
	}else if(isFlatRoof(buildingCell, segments)){
		buildingCell.roofType = Flat;
		buildingCell.roofAngle = 0;
		buildingCell.roofHeight = segments.at(buildingCell.segmentNumbers.at(0)).averageHeight;
		//printf("Roof part %d of polygon %d is flat\n",buildingPart,iBuilding);
	}else if(isShedRoof(buildingCell, segments) 
		&& segments.at(buildingCell.segmentNumbers.at(0)).cells.size()==1
		&& ((OGRPolygon*)cellPolygon)->getExteriorRing()->getNumPoints()==5
		&& isRectangularCell(cellPolygon)){
		buildingCell.roofType = Shed;
		buildingCell.roofDirection = segments.at(buildingCell.segmentNumbers.at(0)).normal;
		buildingCell.roofAngle = roundDouble(angleWithHorizontalVector(segments.at(buildingCell.segmentNumbers.at(0)).normal));
		buildingCell.roofHeight = segments.at(buildingCell.segmentNumbers.at(0)).averageHeight;
		buildingCell.eavesHeight = calculateEavesHeight(segments,buildingCell.segmentNumbers);
		buildingCell.ridgeHeight = calculateRidgeHeight(segments,buildingCell.segmentNumbers);
		//printf("Roof part %d of polygon %d is shed\n",buildingPart,iBuilding);
	}else if(isGabledRoof(buildingCell,segments)
		&& segments.at(buildingCell.segmentNumbers.at(0)).cells.size()==1 //segments do only cover this cell
		&& segments.at(buildingCell.segmentNumbers.at(1)).cells.size()==1
		&& ((OGRPolygon*)cellPolygon)->getExteriorRing()->getNumPoints()==5
		&& isRectangularCell(cellPolygon)){
		buildingCell.roofType = Gabled;
		buildingCell.roofDirection = getRoofDirection(segments, buildingCell.segmentNumbers);
		buildingCell.roofAngle = roundDouble(angleWithHorizontalVector(segments.at(buildingCell.segmentNumbers.at(0)).normal));
		buildingCell.roofHeight = calculateWeightedAverageCellHeight(segments, buildingCell.segmentNumbers, buildingCell.roof);
		buildingCell.eavesHeight = calculateEavesHeight(segments,buildingCell.segmentNumbers);
		buildingCell.ridgeHeight = calculateRidgeHeight(segments,buildingCell.segmentNumbers);
		//printf("Roof part %d of polygon %d is gabled\n",buildingPart,iBuilding);
	}else{
		buildingCell.roofType = UnknownType;
		buildingCell.roofAngle = 0;
		buildingCell.roofHeight = calculateHeightMultipleSegments(segments, buildingCell.segmentNumbers);
		buildingCell.roofHeight = calculateWeightedAverageCellHeight(segments, buildingCell.segmentNumbers, buildingCell.roof);
	}
}