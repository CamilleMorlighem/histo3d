#include "Segments.h"

using namespace std;

//Calculate Heights
double calculateEavesHeight(map<int,Segment> segments,vector<int> segmentNumbers){
	double minimumHeight = 9999;
	for(vector<int>::iterator segNr = segmentNumbers.begin();segNr!=segmentNumbers.end();++segNr){
		for(int pointNr=0;pointNr<segments.at(*segNr).points->getNumGeometries();pointNr++){
			if(minimumHeight > ((OGRPoint*)segments.at(*segNr).points->getGeometryRef(pointNr))->getZ()){
				minimumHeight = ((OGRPoint*)segments.at(*segNr).points->getGeometryRef(pointNr))->getZ();
			}
		}
	}
	return ceil(minimumHeight*100)/100;
}

double calculateRidgeHeight(map<int,Segment> segments,vector<int> segmentNumbers){
	double maximumHeight = 0;
	for(vector<int>::iterator segNr = segmentNumbers.begin();segNr!=segmentNumbers.end();++segNr){
		for(int pointNr=0;pointNr<segments.at(*segNr).points->getNumGeometries();pointNr++){
			if(maximumHeight < ((OGRPoint*)segments.at(*segNr).points->getGeometryRef(pointNr))->getZ()){
				maximumHeight = ((OGRPoint*)segments.at(*segNr).points->getGeometryRef(pointNr))->getZ();
			}
		}
	}
	return ceil(maximumHeight*100)/100;
}

double calculateAverageHeight(Segment segment){
	double totalHeight = 0;
	for(int pointNr=0;pointNr<segment.points->getNumGeometries();pointNr++){
		totalHeight += ((OGRPoint*)segment.points->getGeometryRef(pointNr))->getZ();
	}
	return ceil((totalHeight/segment.points->getNumGeometries())*100)/100;
}

void calculateSegmentHeights(map<int,Segment> &segments){
	for(map<int,Segment>::iterator segment = segments.begin();segment!=segments.end();++segment){
		segment->second.averageHeight = calculateAverageHeight(segment->second);
	}
}

double calculateWeightedAverageCellHeight(map<int,Segment> &segments, vector<int> segmentNumbers, OGRMultiPolygon* roof){
	//OGRMultiPoint* points = new OGRMultiPoint;
	double totalHeight = 0;
	int totalPoints = 0;

	for(vector<int>::iterator segNr = segmentNumbers.begin();segNr!=segmentNumbers.end();++segNr){
		for(int pointNr=0;pointNr<segments.at(*segNr).points->getNumGeometries();pointNr++){
			if(segments.at(*segNr).points->getGeometryRef(pointNr)->Within(roof)){
				totalHeight += ((OGRPoint*)segments.at(*segNr).points->getGeometryRef(pointNr))->getZ();
				totalPoints++;
			}
		}
	}
	return ceil((totalHeight/totalPoints)*100)/100;
}

double calculateHeightMultipleSegments(map<int,Segment> &segments, vector<int> segmentNumbers){
	double averageHeight = 0;
	int pointCount = 0;

	for(vector<int>::iterator segNr = segmentNumbers.begin();segNr!=segmentNumbers.end();++segNr){
		averageHeight += segments.at(*segNr).averageHeight*segments.at(*segNr).points->getNumGeometries();
		pointCount += segments.at(*segNr).points->getNumGeometries();
	}
	return ceil((averageHeight/pointCount)*100)/100;
}

map<int,Segment> findSegmentsInBuilding(OGRLayer* pointsLayer, OGRGeometry* buildingGeometry){
	map<int,Segment> segments;
	if(buildingGeometry != NULL){
		pointsLayer->SetSpatialFilter(buildingGeometry->UnionCascaded());

		OGRFeature* feat = OGRFeature::CreateFeature(pointsLayer->GetLayerDefn());
		while( (feat = pointsLayer->GetNextFeature()) != NULL){
			Segment tempSeg;
			tempSeg.id = feat->GetFieldAsInteger("Segment_id");
			tempSeg.normal.push_back(feat->GetFieldAsInteger("Normal_X"));
			tempSeg.normal.push_back(feat->GetFieldAsInteger("Normal_Y"));
			tempSeg.normal.push_back(feat->GetFieldAsInteger("Normal_Z"));
			tempSeg.points = (OGRMultiPoint*)feat->GetGeometryRef();
			tempSeg.convexHull = (OGRPolygon*)((OGRMultiPoint*)feat->GetGeometryRef())->ConvexHull();
			segments[tempSeg.id] = tempSeg;
		}
		OGRFeature::DestroyFeature(feat);
	}
	return segments;
}