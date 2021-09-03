#include <vld.h>

//Memory usage
#include <Windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib") //added

#include <algorithm>
#include <ctime>
#include <sys/stat.h>
#include "geos_c.h"
#include "Cell_Decomposition.h"

using namespace std;

bool write_adjacent_buildings = false;
bool write_gen_buf_lines = false;
bool write_split_bb = false;
bool write_gen_polygon = true;
bool write_decomp_poly = true;

const double PI = 4*atan(1.0);

typedef enum {
	NE = 1,
	SE = -1,
	SW = -2,
	NW = 2
} Direction;

struct Line{
	int id;
	double A,B;
	double length;
	vector<OGRPoint> point;
	Direction direction;
	void clear() {
		id=0;
		A=B=length=0;
		point.clear();
	}
};

//class Line{
//public:
//	int id;
//	double A,B;
//	double length;
//	vector<OGRPoint> point;
//	Direction direction;
//	void clear() {
//		id=0;
//		A=B=length=0;
//		point.clear();
//	}
//};

struct Buffer{
	int id;
	double A,B,Bmin,Bmax;
	double length;
	bool marked;
	vector<Line> parallelLines;
	vector<Line> nonParallelLines;
	void clear() {
		id=0;
		A=B=Bmin=Bmax=0;
		length=0;
		marked=false;
		parallelLines.clear();
		nonParallelLines.clear();
	}
};

//class Buffer{
//public:
//	int id;
//	double A,B,Bmin,Bmax;
//	double length;
//	bool marked;
//	vector<Line> parallelLines;
//	vector<Line> nonParallelLines;
//	void clear() {
//		id=0;
//		A=B=Bmin=Bmax=0;
//		length=0;
//		marked=false;
//		parallelLines.clear();
//		nonParallelLines.clear();
//	}
//};

int bufferLengthSorter(Buffer &i,Buffer &j) {return i.length > j.length;}

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

OGRMultiPolygon* cutPolygon(OGRMultiPolygon* Polygon,OGRLineString* CutLine){
	if(CutLine == NULL){
		CPLError( CE_Failure, CPLE_ObjectNull, "Cut Line is null" );
		return NULL;
	}

	//1. Check intersects
	if(!Polygon->Intersects(CutLine)){
		CPLError( CE_Warning, CPLE_AppDefined, "Cut Line does not intersect polygon" );
		return NULL;
	}
	//2. Boundary
	OGRGeometry *pPolyBoundary = Polygon->Boundary();
	if(NULL == pPolyBoundary){
 		CPLError( CE_Failure, CPLE_ObjectNull, "Error get Boundary" );
 		return NULL;
	}
	//3. Union
	OGRGeometry *pNewPoly = pPolyBoundary->Union(CutLine);
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
	OGRMultiPolygon* pOutputGeom = (OGRMultiPolygon*) OGRGeometryFactory::createGeometry(wkbMultiPolygon);
 		 
	for ( int i = 0; i < GetGEOSNumberOfGeometries(polygons); i++ ){
 		OGRGeometry* pGeom = OGRGeometryFactory::createFromGEOS((GEOSGeom)GEOSGetGeometryN(polygons, i));
 		pOutputGeom->addGeometry(pGeom);
	}
	GEOSGeom_destroy(polygons);

	return pOutputGeom;
}

double OGRGeometry_Get_Area(OGRGeometry* inputGeom){
	OGRwkbGeometryType type = wkbFlatten(inputGeom->getGeometryType());
	switch(type){
		case wkbPolygon: {
			OGRPolygon* inputPoly = (OGRPolygon*) inputGeom;
			return inputPoly->get_Area();
		}
		case wkbMultiPolygon: {
			OGRMultiPolygon* inputMultiPoly = (OGRMultiPolygon*) inputGeom;
			return inputMultiPoly->get_Area();
		}
		case wkbGeometryCollection: {
			OGRGeometryCollection* inputGeometryCollection = (OGRGeometryCollection*) inputGeom;
			return inputGeometryCollection->get_Area();
		}
		default:
			return 0;
	}
}

double minBBuffer(double A,Buffer buffer){
	double tempB=0,minB=999999999999;

	for(vector<Line>::iterator line = buffer.parallelLines.begin(); line!=buffer.parallelLines.end();++line){
		for(vector<OGRPoint>::iterator point = line->point.begin(); point != line->point.end(); ++point){
			tempB = point->getY() - A * point->getX();
			if(tempB<minB) minB=tempB;
		}
	}
	return minB;
}

double maxBBuffer(double A,Buffer buffer){
	double tempB=0,maxB=-999999999999;

	for(vector<Line>::iterator line = buffer.parallelLines.begin(); line!=buffer.parallelLines.end();++line){
		for(vector<OGRPoint>::iterator point = line->point.begin(); point != line->point.end(); ++point){
			tempB = point->getY() - A * point->getX();
			if(tempB>maxB) maxB=tempB;
		}
	}

	return maxB;
}

bool joinableBuffers(Buffer buffer1, Buffer buffer2, double epsilon, double alpha){
	double angle = abs(atan(buffer1.A) - atan(buffer2.A));
	alpha = alpha / (180/PI); //From degrees to radians

	if(angle < alpha){
		//double weightangle = cos(atan(buffer2.A) - atan(buffer1.A));
	
		//double weight1 = buffer1.length						/ (buffer1.length + abs(buffer2.length*weightangle));
		//double weight2 = abs(buffer2.length*weightangle)	/ (buffer1.length + abs(buffer2.length*weightangle));
		double weight1 = buffer1.length	/ (buffer1.length + buffer2.length);
		double weight2 = buffer2.length	/ (buffer1.length + buffer2.length);

		//Update Ax+B
		double tempA = tan(atan(buffer1.A) * weight1 + atan(buffer2.A) * weight2);
	
		double minB1=minBBuffer(tempA, buffer1);
		double maxB1=maxBBuffer(tempA, buffer1);
		double minB=minBBuffer(tempA, buffer2);
		double maxB=maxBBuffer(tempA, buffer2);

		if(minB<minB1) minB1=minB;
		if(maxB>maxB1) maxB1=maxB;
	
		if(abs((maxB1-minB1) * cos(atan(tempA))) < epsilon){
			return true;
		}
	}
	return false;
}

bool joinableBuffersAdjacentBuilding(Buffer buffer1, Buffer buffer2, double epsilon, double alpha){
	double angle = abs(atan(buffer1.A) - atan(buffer2.A));
	alpha = alpha / (180/PI); //From degrees to radians

	if(angle < alpha || abs(angle-PI) < alpha){
		double weight1 = buffer1.length	/ (buffer1.length + buffer2.length);
		double weight2 = buffer2.length	/ (buffer1.length + buffer2.length);

		//Update Ax+B
		double tempA = tan(atan(buffer1.A) * weight1 + atan(buffer2.A) * weight2);
		
		double minB1=minBBuffer(tempA, buffer1);
		double maxB1=maxBBuffer(tempA, buffer1);
		double minB=minBBuffer(tempA, buffer2);
		double maxB=maxBBuffer(tempA, buffer2);

		if(minB<minB1) minB1=minB;
		if(maxB>maxB1) maxB1=maxB;
	
		if(abs((maxB1-minB1) * cos(atan(tempA))) < epsilon){
			//Check if adjacent building line is within epsilon distance
			OGRLineString buffer1Line;
			OGRLineString buffer2Line;
			for(vector<Line>::iterator line=buffer1.parallelLines.begin();line<buffer1.parallelLines.end();line++){
				buffer1Line.addPoint(&line->point.at(0));
				buffer1Line.addPoint(&line->point.at(1));
			}
			for(vector<Line>::iterator line=buffer2.parallelLines.begin();line<buffer2.parallelLines.end();line++){
				buffer2Line.addPoint(&line->point.at(0));
				buffer2Line.addPoint(&line->point.at(1));
			}
		
			//if(buffer1Line->Buffer(epsilon,1)->Intersects(buffer2Line)){
			if(buffer1Line.Distance(&buffer2Line) < epsilon){
				return true;
			}
		}
	}
	return false;
}

int joinBuffers(Buffer &buffer1,Buffer buffer2){
	double weight1,weight2;

	////Join parallel lines from both buffers into new buffer
	//for(vector<Line>::iterator it = buffer1.parallelLines.begin(); it!=buffer1.parallelLines.end();++it){
	//	bufferout.parallelLines.push_back(*it);
	//}
	//for(vector<Line>::iterator it = buffer2.parallelLines.begin(); it!=buffer2.parallelLines.end();++it){
	//	bufferout.parallelLines.push_back(*it);
	//}

	////Join all lines from both buffers into new buffer
	//for(vector<Line>::iterator it = buffer1.allLines.begin(); it!=buffer1.allLines.end();++it){
	//	bufferout.allLines.push_back(*it);
	//}
	//for(vector<Line>::iterator it = buffer2.allLines.begin(); it!=buffer2.allLines.end();++it){
	//	bufferout.allLines.push_back(*it);
	//}

	for(vector<Line>::iterator it = buffer2.parallelLines.begin(); it!=buffer2.parallelLines.end();++it){
		buffer1.parallelLines.push_back(*it);
	}

	for(vector<Line>::iterator it = buffer2.nonParallelLines.begin(); it!=buffer2.nonParallelLines.end();++it){
		buffer1.nonParallelLines.push_back(*it);
	}

	//Calculate weights for both lines
	//weightangle = cos(atan(buffer2.A) - atan(buffer1.A));
	
	//weight1 = buffer1.length					/ (buffer1.length + abs(buffer2.length*weightangle));
	//weight2 = abs(buffer2.length*weightangle)	/ (buffer1.length + abs(buffer2.length*weightangle));

	weight1 = buffer1.length / (buffer1.length + buffer2.length);
	weight2 = buffer2.length / (buffer1.length + buffer2.length);
	
	//Set id for new line to id of line1
	//if(buffer1.id>buffer2.id) buffer1.id=buffer2.id;

	//Set marked to false
	//bufferout.marked = false;
	buffer1.marked = false;

	//Update Ax+B
	//bufferout.A = buffer1.A * weight1 + buffer2.A * weight2;
	buffer1.A = tan(atan(buffer1.A) * weight1 + atan(buffer2.A) * weight2);

	//bufferout.Bmin=minBBuffer(bufferout.A,bufferout);
	//bufferout.Bmax=maxBBuffer(bufferout.A,bufferout);
	buffer1.Bmin=minBBuffer(buffer1.A,buffer1);
	buffer1.Bmax=maxBBuffer(buffer1.A,buffer1);

	//weightangle = cos(atan(buffer2.A) - atan(buffer1.A));
	//Calculate length of new line
	//bufferout.length = buffer1.length + buffer2.length*weightangle;
	buffer1.length = buffer1.length + buffer2.length;

	return 1;
}

bool includableBuffers(Buffer buffer1, Buffer buffer2, double epsilon){
	//Include line if it keeps the buffer smaller then the maximum generalization distance
	//Lines that are not within two parallel lines are also included
	double minB=minBBuffer(buffer1.A, buffer2);
	double maxB=maxBBuffer(buffer1.A, buffer2);

	if(buffer1.Bmin<minB) minB=buffer1.Bmin;
	if(buffer1.Bmax>maxB) maxB=buffer1.Bmax;
	
	if( abs((maxB-minB) * cos(atan(buffer1.A))) < epsilon ){
		return true;
	}

	return false;
}

bool includableBuffers2(Buffer buffer1, Buffer buffer2){
	//Include line when it fits within the buffer as is
	//Lines that are not within two parallel lines are NOT included
	if(buffer1.Bmin <= minBBuffer(buffer1.A, buffer2) && buffer1.Bmax >= maxBBuffer(buffer1.A, buffer2)){
		return true;
	}
	return false;
}

bool includableLineWithBuffer(Buffer buffer1, Line line1, double epsilon){
	double tempB=0,minB=999999999,maxB=-999999999;

	tempB = line1.point.at(0).getY() - buffer1.A * line1.point.at(0).getX();
	if(tempB<minB) minB=tempB;
	if(tempB>maxB) maxB=tempB;
	tempB = line1.point.at(1).getY() - buffer1.A * line1.point.at(1).getX();
	if(tempB<minB) minB=tempB;
	if(tempB>maxB) maxB=tempB;

	if(buffer1.Bmin<minB && buffer1.Bmax>maxB){
		return true;
	}

	//if(buffer1.Bmin<minB) minB=buffer1.Bmin;
	//if(buffer1.Bmax>maxB) maxB=buffer1.Bmax;
	//
	//if( abs((maxB-minB) * cos(atan(buffer1.A))) < epsilon ){
	//	return true;
	//}
	return false;
}

int includeBuffers(Buffer &buffer1,Buffer buffer2){
	//Join all lines from both buffers into new buffer
	//for(vector<Line>::iterator it = buffer1.allLines.begin(); it!=buffer1.allLines.end();++it){
	//	bufferout.allLines.push_back(*it);
	//}
	//for(vector<Line>::iterator it = buffer2.allLines.begin(); it!=buffer2.allLines.end();++it){
	//	bufferout.allLines.push_back(*it);
	//}
	for(vector<Line>::iterator it = buffer2.parallelLines.begin(); it!=buffer2.parallelLines.end();++it){
		buffer1.nonParallelLines.push_back(*it);
	}
	for(vector<Line>::iterator it = buffer2.nonParallelLines.begin(); it!=buffer2.nonParallelLines.end();++it){
		buffer1.nonParallelLines.push_back(*it);
	}
	
	//double weightangle = cos(atan(buffer2.A) - atan(buffer1.A));

	//Set id for new line to id of line1
	//if(buffer1.id>buffer2.id) buffer1.id=buffer2.id;

	//Set marked to false
	//bufferout.marked = false;
	buffer1.marked = false;

	//Update Ax+B
	//bufferout.A = buffer1.A;

	//bufferout.Bmin=minBBuffer(bufferout.A,bufferout);
	//bufferout.Bmax=maxBBuffer(bufferout.A,bufferout);
	//NOT needed, buffer B is calculated only from parallel lines
	//buffer1.Bmin=minBBuffer(buffer1.A,buffer1);
	//buffer1.Bmax=maxBBuffer(buffer1.A,buffer1);

	//Calculate length of new line
	//bufferout.length = buffer1.length + buffer2.length*weightangle;
	//NOT needed, length stays the same
	//buffer1.length = buffer1.length + buffer2.length*weightangle;

	return 1;
}

bool writeBufferLine(double A, double B, double xIn1, double xIn2, OGRPolygon *BB, OGRLayer *poLayerOut, int polyid, int lineid){
	double x[2],y[2];
	OGRLineString bufferline;
	int i=0;

	x[0] = xIn1; 
	y[0] = A* xIn1 + B;
	x[1] = xIn2; 
	y[1] = A* xIn2 + B;

	bufferline.addPoint(x[0],y[0]);
	bufferline.addPoint(x[1],y[1]);
	
	OGRFeature *poFeatureOut;
		
	poFeatureOut = OGRFeature::CreateFeature(poLayerOut->GetLayerDefn());
	poFeatureOut->SetField("ID", polyid);
	poFeatureOut->SetField("Line", lineid);
	poFeatureOut->SetGeometry(bufferline.Intersection(BB));
	poLayerOut->CreateFeature(poFeatureOut);

	OGRFeature::DestroyFeature(poFeatureOut);

	return true;
}

void calculateNewBWeighted(vector<Buffer> &generalizedBuffers){
	//Calculate weighted B values for generalized buffers
	for(vector<Buffer>::iterator buffer = generalizedBuffers.begin();buffer!=generalizedBuffers.end();++buffer){
		if(buffer->parallelLines.size()>1){
			buffer->B=0;
			double totallength=0;
			double weight=0;
			double weighttest=0;
			for(vector<Line>::iterator line = buffer->parallelLines.begin();line!=buffer->parallelLines.end();++line){
				totallength += line->length * cos(atan(line->A)-atan(buffer->A));
			}
			for(vector<Line>::iterator line = buffer->parallelLines.begin();line!=buffer->parallelLines.end();++line){
				weight = (line->length * cos(atan(line->A)-atan(buffer->A))) / totallength;
				double tempB = (line->point.at(0).getY() - buffer->A * line->point.at(0).getX()
								+line->point.at(1).getY() - buffer->A * line->point.at(1).getX())/2;
				buffer->B += tempB*weight;
				weighttest += weight;
			}
			if(fabs(weighttest-1) > numeric_limits<double>::epsilon()){
				printf("Weights are incorrect, combined weight is: %.20lf\n",weighttest);
			}
		}
	}
}

void calculateNewBPeter(vector<Buffer> &generalizedBuffers){
	//Calculate B values for generalized buffers using the method of Peter
	for(vector<Buffer>::iterator buffer = generalizedBuffers.begin();buffer!=generalizedBuffers.end();++buffer){
		if(buffer->parallelLines.size()>1){
			vector<Line>::iterator maxLine = buffer->parallelLines.begin();
			for(vector<Line>::iterator line = buffer->parallelLines.begin();line!=buffer->parallelLines.end();++line){
				if(line->length>maxLine->length){
					maxLine=line;
				}
			}
			buffer->B = (maxLine->point.at(0).getY() - buffer->A * maxLine->point.at(0).getX()
							+maxLine->point.at(1).getY() - buffer->A * maxLine->point.at(1).getX())/2;
		}
	}
}

int calculateBufferLine(vector<Buffer> &buffers,OGRRawPoint *points,int numPoints){
	Buffer tempBuffer;
	Line tempLine;
	for(int i=0 ; i<numPoints-1; i++){
		tempLine.id = i;
		tempLine.length = sqrt(pow(points[i+1].x-points[i].x,2)
							+ pow(points[i+1].y-points[i].y,2));

		OGRPoint aPoint;
		aPoint.setX(points[i].x);
		aPoint.setY(points[i].y);
		//Add the first point to the point list
		tempLine.point.push_back(aPoint);
		
		aPoint.setX(points[i+1].x);
		aPoint.setY(points[i+1].y);
		//Add the second point to the point list
		tempLine.point.push_back(aPoint);

		//y=Ax+B gives A=(y1-y2)/-(x1-x2) and B=y1-Ax1
		if(fabs(points[i+1].x-points[i].x) < numeric_limits<double>::epsilon() ){
			if(points[i+1].y < points[i].y) tempLine.A = 9999;	//Line is oriented downwards
			else tempLine.A = -9999;							//Line is oriented upwards
		}else{
			tempLine.A = (points[i+1].y-points[i].y) / (points[i+1].x-points[i].x);
		}
		tempLine.B = points[i].y - tempLine.A*points[i].x;

		//Direction of the line by quadrant
		if(  (points[i+1].y-points[i].y)>0	&&  (points[i+1].x-points[i].x)>0 ) tempLine.direction=NE;
		if( -(points[i+1].y-points[i].y)>0	&&  (points[i+1].x-points[i].x)>0 ) tempLine.direction=SE;
		if( -(points[i+1].y-points[i].y)>0	&& -(points[i+1].x-points[i].x)>0 ) tempLine.direction=SW;
		if(  (points[i+1].y-points[i].y)>0	&& -(points[i+1].x-points[i].x)>0 ) tempLine.direction=NW;

		tempBuffer.A = tempLine.A;
		tempBuffer.B = tempLine.B;
		tempBuffer.Bmin = tempLine.B;
		tempBuffer.Bmax = tempLine.B;
		tempBuffer.id = tempLine.id;
		tempBuffer.length = tempLine.length;
		tempBuffer.marked = false;

		tempBuffer.parallelLines.push_back(tempLine);
		//tempBuffer.allLines.push_back(tempLine);
		buffers.push_back(tempBuffer);

		tempBuffer.clear();
		tempLine.clear();
	}
	return 1;
}

bool linesAreParallelnonCollinear(Buffer buffer, Line line1, Line line2){
	if(abs(buffer.Bmax-buffer.Bmin)>numeric_limits<double>::epsilon()){
		OGRLineString lineString1;
		lineString1.addPoint(&line1.point.at(0));
		lineString1.addPoint(&line1.point.at(1));
		GEOSGeometry* geosline1 = lineString1.exportToGEOS();
		GEOSGeometry* geosbuffer = GEOSBufferWithStyle(geosline1,abs(buffer.Bmax-buffer.Bmin),0,GEOSBUF_CAP_FLAT,GEOSBUF_JOIN_MITRE,0);
		OGRGeometry* linebuffer = OGRGeometryFactory::createFromGEOS(geosbuffer);
		GEOSGeom_destroy(geosline1);
		GEOSGeom_destroy(geosbuffer);
		OGRLineString lineString2;
		lineString2.addPoint(&line2.point.at(0));
		lineString2.addPoint(&line2.point.at(1));
		if(lineString2.Within(linebuffer)){
			return true;
		}
	}
	return false;
}

//bool linesAreParallelnonCollinearWidth(Line line1, Line line2, double bufferWidth){
//	OGRLineString* lineString1 = new OGRLineString;
//	lineString1->addPoint(&line1.point.at(0));
//	lineString1->addPoint(&line1.point.at(1));
//	GEOSGeometry* geosline1 = lineString1->exportToGEOS();
//	GEOSGeometry* geosbuffer = GEOSBufferWithStyle(geosline1,bufferWidth,0,GEOSBUF_CAP_FLAT,GEOSBUF_JOIN_MITRE,0);
//	OGRGeometry* linebuffer = OGRGeometryFactory::createFromGEOS(geosbuffer);
//
//	OGRLineString* lineString2 = new OGRLineString;
//	lineString2->addPoint(&line2.point.at(0));
//	lineString2->addPoint(&line2.point.at(1));
//	if(lineString2->Within(linebuffer)){
//		return true;
//	}
//	return false;
//}

bool linesHaveSameDirection(Line line1, Line line2){
	if(line1.direction == line2.direction) return true;
	
	if(abs(line1.A) < 1){
		if(line1.direction == NE && line2.direction == SE) return true;
		if(line1.direction == SE && line2.direction == NE) return true;
		if(line1.direction == SW && line2.direction == NW) return true;
		if(line1.direction == NW && line2.direction == SW) return true;
	}

	if(abs(line1.A) >= 1){
		if(line1.direction == NE && line2.direction == NW) return true;
		if(line1.direction == SE && line2.direction == SW) return true;
		if(line1.direction == SW && line2.direction == SE) return true;
		if(line1.direction == NW && line2.direction == NE) return true;
	}

	return false;
}

int populateBuffers(OGRPolygon *Polygon, vector<Buffer> &buffers){
	OGRRawPoint* points = new OGRRawPoint[Polygon->getExteriorRing()->getNumPoints()];
	Polygon->getExteriorRing()->getPoints(points);
	calculateBufferLine(buffers, points, Polygon->getExteriorRing()->getNumPoints());

	return 1;
}

int populateBuffersWithInterior(OGRPolygon *Polygon, vector<Buffer> &buffers){
	int totalPoints=Polygon->getExteriorRing()->getNumPoints();
	//printf("interior rings: %d\n",Polygon->getNumInteriorRings());
	for(int i=0;i<Polygon->getNumInteriorRings();i++){
		//printf("interior ring %d has %d points\n",i,Polygon->getInteriorRing(i)->getNumPoints());
		totalPoints += Polygon->getInteriorRing(i)->getNumPoints();
	}

	OGRRawPoint* points = new OGRRawPoint[Polygon->getExteriorRing()->getNumPoints()];
	Polygon->getExteriorRing()->getPoints(points);
	calculateBufferLine(buffers, points,Polygon->getExteriorRing()->getNumPoints());

	for(int i=0;i<Polygon->getNumInteriorRings();i++){
		points = new OGRRawPoint[Polygon->getInteriorRing(i)->getNumPoints()];
		Polygon->getInteriorRing(i)->getPoints(points);
		calculateBufferLine(buffers,points,Polygon->getInteriorRing(i)->getNumPoints());
	}
	
	return 1;
}

int populateBuffersWithInterior(OGRGeometry *Geometry, vector<Buffer> &buffers){
	if(wkbFlatten(Geometry->getGeometryType()) == wkbPolygon){
		OGRPolygon* Polygon = (OGRPolygon*) Geometry;

		OGRRawPoint* points = new OGRRawPoint[Polygon->getExteriorRing()->getNumPoints()];
		Polygon->getExteriorRing()->getPoints(points);
		calculateBufferLine(buffers, points,Polygon->getExteriorRing()->getNumPoints());

		for(int i=0;i<Polygon->getNumInteriorRings();i++){
			points = new OGRRawPoint[Polygon->getInteriorRing(i)->getNumPoints()];
			Polygon->getInteriorRing(i)->getPoints(points);
			calculateBufferLine(buffers,points,Polygon->getInteriorRing(i)->getNumPoints());
		}
		delete [] points;
	}else if(wkbFlatten(Geometry->getGeometryType()) == wkbMultiPolygon){
		OGRMultiPolygon* MultiPolygon = (OGRMultiPolygon*) Geometry;
		for(int i=0;i<MultiPolygon->getNumGeometries();i++){
			if(wkbFlatten(MultiPolygon->getGeometryRef(i)->getGeometryType()) == wkbPolygon){
				OGRPolygon* Polygon = (OGRPolygon*) MultiPolygon->getGeometryRef(i);

				OGRRawPoint* points = new OGRRawPoint[Polygon->getExteriorRing()->getNumPoints()];
				Polygon->getExteriorRing()->getPoints(points);
				calculateBufferLine(buffers, points,Polygon->getExteriorRing()->getNumPoints());

				for(int j=0;j<Polygon->getNumInteriorRings();j++){
					points = new OGRRawPoint[Polygon->getInteriorRing(j)->getNumPoints()];
					Polygon->getInteriorRing(j)->getPoints(points);
					calculateBufferLine(buffers,points,Polygon->getInteriorRing(j)->getNumPoints());
				}
			}
		}
	}else if(wkbFlatten(Geometry->getGeometryType()) == wkbGeometryCollection){
		OGRGeometryCollection* geomCollection = (OGRGeometryCollection*) Geometry;
		for(int geom=0;geom<geomCollection->getNumGeometries();geom++){
			populateBuffersWithInterior(geomCollection->getGeometryRef(geom),buffers);
		}

	}else throw "PopulateBuffersWithInterior, not a Polygon or MultiPolygon";

	return 1;
}

void joinCollinearBuffers(vector<Buffer> &generalizedBuffers,double epsilon,double alpha){
	//TODO: INCOMPLETE! Only joining buffer on epsilon now, same as joinBuffers
	int index=0;
	//generalizedBuffers.at(index)
	//generalizedBuffers.at(index+1)
	while(index<(int)generalizedBuffers.size()-1){
		if(joinableBuffers(generalizedBuffers.at(index),generalizedBuffers.at(index+1),epsilon,alpha)){
			joinBuffers(generalizedBuffers.at(index),generalizedBuffers.at(index+1));
			generalizedBuffers.erase(generalizedBuffers.begin()+index+1);
		}else index++;
	}
	//test for last line with first line
	if(joinableBuffers(generalizedBuffers.at(0),generalizedBuffers.at(index),epsilon,alpha)){
		joinBuffers(generalizedBuffers.at(0),generalizedBuffers.at(index));
		generalizedBuffers.erase(generalizedBuffers.begin()+index);
	}
}

int generalizeBuffers(vector<Buffer> &generalizedBuffers,double epsilon,double alpha){
	//1: test all joinable, 2: test all includable, 3: test all lines in smaller buffers if includable
	sort(generalizedBuffers.begin(),generalizedBuffers.end(),bufferLengthSorter);

	int index = 0;
	//Test if buffers can be joined based on alpha and epsilon
	while(index<(int)generalizedBuffers.size()-1){
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			if(joinableBuffers(generalizedBuffers.at(index),*buffer,epsilon,alpha)){
				joinBuffers(generalizedBuffers.at(index),*buffer);
				buffer = generalizedBuffers.erase(buffer);
			}else buffer++;
		}
		index++;
	}

	index=0;
	while(index<(int)generalizedBuffers.size()-1){
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			if(includableBuffers(generalizedBuffers.at(index),*buffer,epsilon)){
				includeBuffers(generalizedBuffers.at(index),*buffer);
				buffer = generalizedBuffers.erase(buffer);
			}
			else buffer++;
		}
		index++;
	}
		
	index=0;
	//Test if lines from other buffers can be included based on epsilon
	while(index<(int)generalizedBuffers.size()-1){
		sort(generalizedBuffers.begin(),generalizedBuffers.end(),bufferLengthSorter);

		//index=0;
		while(generalizedBuffers.at(0).marked && generalizedBuffers.at(index).marked && index<(int)generalizedBuffers.size()-1){
			index++;
		}

		generalizedBuffers.at(index).marked=true;

		//Loop over buffers
		for(vector<Buffer>::iterator k=generalizedBuffers.begin();k!=generalizedBuffers.end();){
			bool updateBufferit = true;
			if(k->id!=generalizedBuffers.at(index).id){
				//Loop over allLines to find includable lines
				for(vector<Line>::iterator l=k->nonParallelLines.begin();l!=k->nonParallelLines.end();){
				//for(vector<Line>::iterator l=k->allLines.begin();l!=k->allLines.end();++l){
					bool updateLineit = true;
					if(includableLineWithBuffer(generalizedBuffers.at(index),*l,epsilon)){
						if(k->nonParallelLines.size()==1){
							k = generalizedBuffers.erase(k);
							updateBufferit = false;
							break;
						}
						
						//Loop over parallelLines to remove line if in there
						for(vector<Line>::iterator p=k->parallelLines.begin();p!=k->parallelLines.end();++p){
							if(p->id == l->id){
								k->parallelLines.erase(p);
								break;
							}
						}

						l = k->nonParallelLines.erase(l);
						updateLineit = false;

					}
					if(l!=k->nonParallelLines.end() && updateLineit) l++;
				}
				if(updateBufferit){
					k->nonParallelLines.shrink_to_fit();
					k->parallelLines.shrink_to_fit();
				}
			}
			if(k!=generalizedBuffers.end() && updateBufferit) k++;
		}
		generalizedBuffers.shrink_to_fit();
	}
	//calculateNewBWeighted(generalizedBuffers);
	calculateNewBPeter(generalizedBuffers);
	return 1;
}

int generalizeBuffers2(vector<Buffer> &generalizedBuffers,double epsilon,double alpha){
	//1: test all joinable, 2: test all includable   //, 3: test all lines in smaller buffers if includable
	sort(generalizedBuffers.begin(),generalizedBuffers.end(),bufferLengthSorter);

	int index = 0;
	//Test if buffers can be joined based on alpha and epsilon
	while(index<(int)generalizedBuffers.size()-1){
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			if(joinableBuffers(generalizedBuffers.at(index),*buffer,epsilon,alpha)){
				joinBuffers(generalizedBuffers.at(index),*buffer);
				buffer = generalizedBuffers.erase(buffer);
			}else buffer++;
		}
		index++;
	}

	index=0;
	//Test if buffers can be included, if they fall completely within the larger buffer
	while(index<(int)generalizedBuffers.size()-1){
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			if(includableBuffers2(generalizedBuffers.at(index),*buffer)){
				includeBuffers(generalizedBuffers.at(index),*buffer);
				buffer = generalizedBuffers.erase(buffer);
			}
			else buffer++;
		}
		index++;
	}
	
	//index=0;
	//while(index<(int)generalizedBuffers.size()-1){
	//	sort(generalizedBuffers.begin(),generalizedBuffers.end(),bufferLengthSorter);

	//	for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
	//		bool updateBufferit = true;
	//		for(vector<Line>::iterator l=buffer->allLines.begin();l!=buffer->allLines.end();){
	//			if(includableLineWithBuffer(generalizedBuffers.at(index),*l,epsilon)){
	//				//If buffer is empty after line removal delete buffer
	//				if(buffer->allLines.size()==1){
	//					buffer = generalizedBuffers.erase(buffer);
	//					updateBufferit = false;
	//					break;
	//				}

	//				//Remove line from parallelLines vector
	//				for(vector<Line>::iterator p=buffer->parallelLines.begin();p!=buffer->parallelLines.end();++p){
	//					if(p->id == l->id){
	//						buffer->parallelLines.erase(p);
	//						break;
	//					}
	//				}
	//				//Remove line from allLines vector
	//				l = buffer->allLines.erase(l);
	//			}
	//			else l++;
	//		}
	//		if(updateBufferit){
	//			buffer++;
	//		}
	//	}
	//	index++;
	//}

	//calculateNewBWeighted(generalizedBuffers);
	calculateNewBPeter(generalizedBuffers);
	return 1;
}

int generalizeBuffers3(vector<Buffer> &generalizedBuffers,double epsilon,double alpha){
	//1: test all joinable with same directions, test all includable 
	//2: test all joinable with both directions, test all includable 
	sort(generalizedBuffers.begin(),generalizedBuffers.end(),bufferLengthSorter);

	int index = 0;
	//Test if buffers can be joined based on alpha and epsilon
	while(index<(int)generalizedBuffers.size()-1){
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			if(linesHaveSameDirection(generalizedBuffers.at(index).parallelLines.at(0),buffer->parallelLines.at(0))
				&& joinableBuffers(generalizedBuffers.at(index),*buffer,epsilon,alpha)){
					joinBuffers(generalizedBuffers.at(index),*buffer);
					buffer = generalizedBuffers.erase(buffer);
			}else buffer++;
		}
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			if(includableBuffers2(generalizedBuffers.at(index),*buffer)){
				includeBuffers(generalizedBuffers.at(index),*buffer);
				buffer = generalizedBuffers.erase(buffer);
			}
			else buffer++;
		}
		index++;
	}

	index=0;
	while(index<(int)generalizedBuffers.size()-1){
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			if(joinableBuffers(generalizedBuffers.at(index),*buffer,epsilon,alpha)){
				joinBuffers(generalizedBuffers.at(index),*buffer);
				buffer = generalizedBuffers.erase(buffer);
			}else buffer++;
		}
		//Test if buffers can be included, if they fall completely within the buffer
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			if(includableBuffers2(generalizedBuffers.at(index),*buffer)){
				includeBuffers(generalizedBuffers.at(index),*buffer);
				buffer = generalizedBuffers.erase(buffer);
			}
			else buffer++;
		}
		index++;
	}

	//calculateNewBWeighted(generalizedBuffers);
	calculateNewBPeter(generalizedBuffers);
	return 1;
}

int generalizeBuffersWithAdjacentBuildings(vector<Buffer> &generalizedBuffers,vector<Buffer> &adjacentBuffers,double epsilon,double alpha){
	//1: test all joinable with same directions, test all includable with adjacent buildings
	//2: test all joinable with both directions, test all includable with adjacent buildings
	sort(generalizedBuffers.begin(),generalizedBuffers.end(),bufferLengthSorter);

	int index = 0;
	//Test if buffers can be joined based on alpha and epsilon
	while(index<(int)generalizedBuffers.size()){
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			if(linesHaveSameDirection(generalizedBuffers.at(index).parallelLines.at(0),buffer->parallelLines.at(0))
				&& joinableBuffers(generalizedBuffers.at(index),*buffer,epsilon,alpha)){
					joinBuffers(generalizedBuffers.at(index),*buffer);
					buffer = generalizedBuffers.erase(buffer);
			}else buffer++;
		}
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			if(includableBuffers2(generalizedBuffers.at(index),*buffer)){
				includeBuffers(generalizedBuffers.at(index),*buffer);
				buffer = generalizedBuffers.erase(buffer);
			}
			else buffer++;
		}
		index++;
	}

	index=0;
	while(index<(int)generalizedBuffers.size()){
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			if(joinableBuffers(generalizedBuffers.at(index),*buffer,epsilon,alpha)){
				joinBuffers(generalizedBuffers.at(index),*buffer);
				buffer = generalizedBuffers.erase(buffer);
			}else buffer++;
		}
		for(vector<Buffer>::iterator adjacentBuildingBuffer=adjacentBuffers.begin();adjacentBuildingBuffer!=adjacentBuffers.end();){
			if(joinableBuffersAdjacentBuilding(generalizedBuffers.at(index),*adjacentBuildingBuffer,epsilon,alpha)){
					joinBuffers(generalizedBuffers.at(index),*adjacentBuildingBuffer);
					adjacentBuildingBuffer = adjacentBuffers.erase(adjacentBuildingBuffer);
			}else adjacentBuildingBuffer++;
		}
		//Test if buffers can be included, if they fall completely within the buffer
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			if(includableBuffers2(generalizedBuffers.at(index),*buffer)){
				includeBuffers(generalizedBuffers.at(index),*buffer);
				buffer = generalizedBuffers.erase(buffer);
			}
			else buffer++;
		}
		for(vector<Buffer>::iterator adjacentBuildingBuffer=adjacentBuffers.begin();adjacentBuildingBuffer!=adjacentBuffers.end();){
			if(includableBuffers2(generalizedBuffers.at(index),*adjacentBuildingBuffer)){
				includeBuffers(generalizedBuffers.at(index),*adjacentBuildingBuffer);
				adjacentBuildingBuffer = adjacentBuffers.erase(adjacentBuildingBuffer);
			}
			else adjacentBuildingBuffer++;
		}
		index++;
	}

	//calculateNewBWeighted(generalizedBuffers);
	calculateNewBPeter(generalizedBuffers);
	return 1;
}

int rightAnglesBuffers(vector<Buffer> &generalizedBuffers,double rightAlpha){
	double rightAngle = 90 / (180/PI);
	rightAlpha = rightAlpha / (180/PI);
	
	sort(generalizedBuffers.begin(),generalizedBuffers.end(),bufferLengthSorter);

	int index=0;
	while(index<(int)generalizedBuffers.size()-1){
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();++buffer){
			double angle = abs(atan(generalizedBuffers.at(index).A) - atan(buffer->A));
			if(abs(rightAngle - angle) < rightAlpha){
				if(abs(atan(generalizedBuffers.at(index).A)) < numeric_limits<double>::epsilon() ){ //for vertical lines test for horizontal + 90deg
					if(buffer->A < 0) buffer->A = -9999;
					else buffer->A = 9999;
				}else{
					if(buffer->A < 0) buffer->A = tan( atan(generalizedBuffers.at(index).A) - rightAngle);
					else buffer->A = tan( atan(generalizedBuffers.at(index).A) + rightAngle);
				}

				////Updating B
				//double totallength=0;
				//double weight=0;
				//double weighttest=0;
				//for(vector<Line>::iterator line = buffer->parallelLines.begin();line!=buffer->parallelLines.end();++line){
				//	totallength += line->length * cos(atan(line->A)-atan(buffer->A));
				//}
				//buffer->B=0;
				//for(vector<Line>::iterator line = buffer->parallelLines.begin();line!=buffer->parallelLines.end();++line){
				//	weight = (line->length * cos(atan(line->A)-atan(buffer->A))) / totallength;
				//	double tempB = ( line->point.at(0).getY() - buffer->A * line->point.at(0).getX()
				//					+line->point.at(1).getY() - buffer->A * line->point.at(1).getX()) /2;
				//	buffer->B += tempB*weight;
				//	weighttest += weight;
				//}
				//if(fabs(weighttest-1) > numeric_limits<double>::epsilon()){
				//	printf("Weights are incorrect, combined weight is: %.20lf\n",weighttest);
				//}

				//Updating B Peter
				vector<Line>::iterator maxLine = buffer->parallelLines.begin();
				for(vector<Line>::iterator line = buffer->parallelLines.begin();line!=buffer->parallelLines.end();++line){
					if(line->length>maxLine->length){
						maxLine=line;
					}
				}
				buffer->B = (maxLine->point.at(0).getY() - buffer->A * maxLine->point.at(0).getX()
							+maxLine->point.at(1).getY() - buffer->A * maxLine->point.at(1).getX())/2;
			}
		}
		index++;
	}
	return 1;
}

int parallelizeBuffers(vector<Buffer> &generalizedBuffers,double parallelAlpha){
	parallelAlpha = parallelAlpha / (180/PI);

	sort(generalizedBuffers.begin(),generalizedBuffers.end(),bufferLengthSorter);

	int index=0;
	while(index<(int)generalizedBuffers.size()-1){
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();++buffer){
			double angle = abs(atan(generalizedBuffers.at(index).A) - atan(buffer->A));
			if(angle < parallelAlpha){
				buffer->A = generalizedBuffers.at(index).A;

				////Updating B Kada
				//double totallength=0;
				//double weight=0;
				//double weighttest=0;
				//for(vector<Line>::iterator line = buffer->parallelLines.begin();line!=buffer->parallelLines.end();++line){
				//	totallength += line->length * cos(atan(line->A)-atan(buffer->A));
				//}

				//buffer->B=0;

				//for(vector<Line>::iterator line = buffer->parallelLines.begin();line!=buffer->parallelLines.end();++line){
				//	weight = (line->length * cos(atan(line->A)-atan(buffer->A))) / totallength;
				//	double tempB = ( line->point.at(0).getY() - buffer->A * line->point.at(0).getX()
				//					+line->point.at(1).getY() - buffer->A * line->point.at(1).getX()) /2;
				//	buffer->B += tempB*weight;
				//	weighttest += weight;
				//}
				//if(fabs(weighttest-1) > numeric_limits<double>::epsilon()){
				//	printf("Weights are incorrect, combined weight is: %.20lf\n",weighttest);
				//}

				//Updating B Peter
				vector<Line>::iterator maxLine = buffer->parallelLines.begin();
				for(vector<Line>::iterator line = buffer->parallelLines.begin();line!=buffer->parallelLines.end();++line){
					if(line->length>maxLine->length){
						maxLine=line;
					}
				}
				buffer->B = (maxLine->point.at(0).getY() - buffer->A * maxLine->point.at(0).getX()
							+maxLine->point.at(1).getY() - buffer->A * maxLine->point.at(1).getX())/2;
			}
		}
		index++;
	}
	return 1;
}

int removeCollapsedLines(vector<Buffer> &generalizedBuffers, double epsilon){
	for(vector<Buffer>::iterator buffer = generalizedBuffers.begin();buffer<generalizedBuffers.end();){
		double length=0;
		int index=0;
		length = buffer->parallelLines.at(index).length;
		while(index<(int)buffer->parallelLines.size()-1){
			for(vector<Line>::iterator line = buffer->parallelLines.begin()+index+1;line<buffer->parallelLines.end();++line){
				if(linesAreParallelnonCollinear(*buffer,buffer->parallelLines.at(index),*line)){
					printf("Buffer %d: lines %d and %d are parallel and non-collinear\n",buffer->id,buffer->parallelLines.at(index).id,line->id);	
				}else{
					length += line->length * cos(atan(line->A)-atan(buffer->A));
				}
			}
			index++;
		}
		//Remove buffers with single lines when smaller then epsilon
		//if(buffer->parallelLines.size()>1 && length<epsilon){
		if(length<epsilon){
			printf("Removing small buffer %d\n",buffer->id);
			buffer = generalizedBuffers.erase(buffer);
		}else ++buffer;
	}
	return 1;
}

double findParallelBase(vector<Buffer> generalizedBuffers){
	double parallelAlpha = 10 / (180/PI);
	double rightAngle = 90 / (180/PI);
	double rightAlpha = 10 / (180/PI);
	int maxCount=0;
	int id=0;
	double A=0;

	for(vector<Buffer>::iterator bufferBase=generalizedBuffers.begin();bufferBase!=generalizedBuffers.end();++bufferBase){
		int count=0;
		for(vector<Buffer>::iterator buffer=generalizedBuffers.begin();buffer!=generalizedBuffers.end();++buffer){
			if(bufferBase->id!=buffer->id){
				//double angle = abs(atan(generalizedBuffers.at(0).A) - atan(buffer->A));
				double angle = abs(atan(bufferBase->A) - atan(buffer->A));
				if(angle < parallelAlpha || abs(rightAngle - angle) < rightAlpha){
					count++;
				}
			}
		}
		if(count>maxCount){
			if(count==generalizedBuffers.size()){
				return bufferBase->A;
			}else{
				maxCount=count;
				A=bufferBase->A;
			}
		}	
	}
	return A;
}

OGRMultiPolygon* cutBoundingBoxSlantLast(OGRPolygon* DoubleBBox, vector<Buffer> &generalizedBuffers, double epsilon){
	double BaseBufferA = findParallelBase(generalizedBuffers);
	//Mark parallel lines and lines with right angles
	//SHPWriter* testWriter = new SHPWriter("test.shp",wkbMultiPolygon);
	//SHPWriter* lineWriter = new SHPWriter("lines.shp",wkbLineString);
	double parallelAlpha = 10 / (180/PI);
	double rightAngle = 90 / (180/PI);
	double rightAlpha = 10 / (180/PI);
	
	OGRMultiPolygon* splitDoubleBBox = (OGRMultiPolygon*) OGRGeometryFactory::createGeometry(wkbMultiPolygon);
	splitDoubleBBox->addGeometry(DoubleBBox);
	for(vector<Buffer>::iterator buffer=generalizedBuffers.begin();buffer!=generalizedBuffers.end();++buffer){
		//double angle = abs(atan(generalizedBuffers.at(0).A) - atan(buffer->A));
		double angle = abs(atan(BaseBufferA) - atan(buffer->A));
		if(angle < parallelAlpha || abs(rightAngle - angle) < rightAlpha){
			buffer->marked = true;

			double x[2],y[2];
			OGRLineString splitline;

			x[0] = DoubleBBox->getExteriorRing()->getX(0);
			y[0] = buffer->A * x[0] + buffer->B;
			x[1] = DoubleBBox->getExteriorRing()->getX(2);
			y[1] = buffer->A * x[1] + buffer->B;

			splitline.addPoint(x[0],y[0]);
			splitline.addPoint(x[1],y[1]);

			splitDoubleBBox = cutPolygon(splitDoubleBBox,&splitline);
			//splitDoubleBBox->dumpReadable(NULL);
		}
	}
	
	OGRMultiPolygon intersectCell;
	//Split double Bounding Box with remainder lines using overlay of the origional lines with the cells
	for(vector<Buffer>::iterator buffer = generalizedBuffers.begin();buffer!=generalizedBuffers.end();++buffer){
		if(!buffer->marked){
			OGRLineString inputLine;
			//Create line to find cells containing this buffer
			//test if line is vertical
			if(buffer->A < -9998 || buffer->A > 9998){
				double min_y=999999999;
				double max_y=0;
				for(vector<Line>::iterator line = buffer->parallelLines.begin();line!=buffer->parallelLines.end();++line){
					if(line->point.at(0).getY() < min_y) min_y = line->point.at(0).getY();
					if(line->point.at(1).getY() < min_y) min_y = line->point.at(1).getY();
					if(line->point.at(0).getY() > max_y) max_y = line->point.at(0).getY();
					if(line->point.at(1).getY() > max_y) max_y = line->point.at(1).getY();
				
				}
				inputLine.addPoint((max_y - buffer->B)/buffer->A ,min_y);
				inputLine.addPoint((max_y - buffer->B)/buffer->A ,max_y);
			}else{
				double min_x=999999999;
				double max_x=0;
				for(vector<Line>::iterator line = buffer->parallelLines.begin();line!=buffer->parallelLines.end();++line){
					if(line->point.at(0).getX() < min_x) min_x = line->point.at(0).getX();
					if(line->point.at(1).getX() < min_x) min_x = line->point.at(1).getX();
					if(line->point.at(0).getX() > max_x) max_x = line->point.at(0).getX();
					if(line->point.at(1).getX() > max_x) max_x = line->point.at(1).getX();
				}
				inputLine.addPoint(min_x,buffer->A * min_x + buffer->B);
				inputLine.addPoint(max_x,buffer->A * max_x + buffer->B);
			}

			//Test for intersecting cell
			for(int cell=0;cell<splitDoubleBBox->getNumGeometries();cell++){
				if(inputLine.Intersects(splitDoubleBBox->getGeometryRef(cell))){
					OGRGeometry* intersectLines = inputLine.Intersection(splitDoubleBBox->getGeometryRef(cell));
					if(wkbFlatten(intersectLines->getGeometryType()) == wkbLineString){
						if(((OGRLineString*)intersectLines)->get_Length() > 0.01){
							intersectCell.addGeometry(splitDoubleBBox->getGeometryRef(cell));
						}
					}else if(wkbFlatten(intersectLines->getGeometryType()) == wkbMultiLineString){
						OGRMultiLineString* multiLines = (OGRMultiLineString*)intersectLines;
						for(int k=0;k<multiLines->getNumGeometries();k++){
							if(((OGRLineString*)multiLines->getGeometryRef(k))->get_Length() > 0.01){
								intersectCell.addGeometry(splitDoubleBBox->getGeometryRef(cell));
							}
						}
					}

				}
			}

			if(!intersectCell.IsEmpty() && intersectCell.getNumGeometries()>0){
				//Find line to split the double BB
				double x[2],y[2];
				OGRLineString splitlineDBB;

				OGREnvelope bufferEnvelope;
				DoubleBBox->Buffer(0.01,0)->getEnvelope(&bufferEnvelope);
				x[0] = bufferEnvelope.MinX;
				//x[0] = DoubleBBox->getExteriorRing()->getX(0);
				y[0] = buffer->A * x[0] + buffer->B;
				x[1] = bufferEnvelope.MaxX;
				//x[1] = DoubleBBox->getExteriorRing()->getX(2);
				y[1] = buffer->A * x[1] + buffer->B;

				splitlineDBB.addPoint(x[0],y[0]);
				splitlineDBB.addPoint(x[1],y[1]);
				//splitlineDBB.dumpReadable(NULL);
				
				OGRGeometry* tempoGeom = ((OGRPolygon*)intersectCell.UnionCascaded())->Buffer(0.001,0)->Boundary()->Intersection(&splitlineDBB);

				//tempoGeom->dumpReadable(NULL);
				if(wkbFlatten(tempoGeom->getGeometryType()) == wkbMultiPoint){
					OGRMultiPoint* points = (OGRMultiPoint*) tempoGeom;
					OGRLineString splitline;
					splitline.addPoint((OGRPoint*)points->getGeometryRef(0));
					splitline.addPoint((OGRPoint*)points->getGeometryRef(1));

					//printf("Splitline \n");
					//splitline.dumpReadable(NULL);
					//OGRFeature* linefeat = OGRFeature::CreateFeature(lineWriter->Layer->GetLayerDefn());
					//linefeat->SetGeometry(&splitline);
					//lineWriter->writeFeature(linefeat);
					//lineWriter->Layer->SyncToDisk();

					splitDoubleBBox = cutPolygon(splitDoubleBBox,&splitline);

					//OGRFeature* feat = OGRFeature::CreateFeature(testWriter->Layer->GetLayerDefn());
					//feat->SetGeometry(splitDoubleBBox);
					//testWriter->writeFeature(feat);
					//testWriter->Layer->SyncToDisk();

				}else printf("cutBoundingBoxSlantLast: tempoGeom is not MultiPoint but %s \n",tempoGeom->getGeometryName());
			}else printf("cutBoundingBoxSlantLast: intersectCell is empty \n");

			intersectCell.empty();
		}
	}

	return splitDoubleBBox;
}

OGRMultiPolygon* cutBoundingBox(OGRPolygon* DoubleBBox, vector<Buffer> &generalizedBuffers){
	OGRMultiPolygon* splitDoubleBBox = new OGRMultiPolygon();
	splitDoubleBBox->addGeometry(DoubleBBox);
	//Cutting double Bounding Box up using generalized lines
	for(vector<Buffer>::iterator buffer = generalizedBuffers.begin();buffer!=generalizedBuffers.end();++buffer){
		double x[2],y[2];
		OGRLineString bufferline;
		int i=0;

		x[0] = DoubleBBox->getExteriorRing()->getX(0); 
		y[0] = buffer->A * x[0] + buffer->B;
		x[1] = DoubleBBox->getExteriorRing()->getX(2); 
		y[1] = buffer->A * x[1] + buffer->B;

		bufferline.addPoint(x[0],y[0]);
		bufferline.addPoint(x[1],y[1]);

		splitDoubleBBox = cutPolygon(splitDoubleBBox,&bufferline);
	}
	return splitDoubleBBox;
}

bool overlapTest(OGRGeometry* intersectGeom,OGRPolygon* decompositionCell,OGRGeometry* input){
	OGRwkbGeometryType type = wkbFlatten(intersectGeom->getGeometryType());
	switch(type){
		case wkbPolygon: {
			OGRPolygon* intersectCell = (OGRPolygon*) intersectGeom;
			if(intersectCell->get_Area() / decompositionCell->get_Area() > 0.5){
				return true;
			}
			return false;
		}
		case wkbMultiPolygon: {
			OGRMultiPolygon* intersectCell = (OGRMultiPolygon*) intersectGeom;
			if(intersectCell->get_Area() / decompositionCell->get_Area() > 0.5){
				return true;
			}
			return false;
		}
		case wkbGeometryCollection: {
			OGRGeometryCollection* geomColl = (OGRGeometryCollection*) intersectGeom;
			for(int k=0;k<geomColl->getNumGeometries();k++){
				//Test geometries in GeometryCollection with overlap test, if true for any geometry, return true
				if(overlapTest(geomColl->getGeometryRef(k),decompositionCell,input)){
					return true;
				}
			}
			return false;
		}
		default:
			return false;
	}
}

int overlapTestSplitBBandInput(OGRMultiPolygon* splitDoubleBBox, OGRGeometry* input){
	//Removing unused cells using overlap test of splitted double Bounding Box with input polygon
	for(int i=0;i<splitDoubleBBox->getNumGeometries();i++){
		if(splitDoubleBBox->getGeometryRef(i) != NULL && wkbFlatten(splitDoubleBBox->getGeometryRef(i)->getGeometryType()) == wkbPolygon){
			OGRPolygon* decompCell = (OGRPolygon*) splitDoubleBBox->getGeometryRef(i);
			OGRGeometry* intersectGeom;

			//Test for input geometry
			if(wkbFlatten(input->getGeometryType()) == wkbMultiPolygon){
				intersectGeom = decompCell->Intersection(input->UnionCascaded());
			}else{
				intersectGeom = decompCell->Intersection(input);
			}

			if(intersectGeom != NULL){
				if(!overlapTest(intersectGeom,decompCell,input)){
					splitDoubleBBox->removeGeometry(i);
					i--;
				}
			}else throw "ERROR: overlapTestSplitBBandInput intersectGeom is NULL";
		}else throw "ERROR: overlapTestSplitBBandInput Subcell is not a polygon";	//printf("ERROR: Subcell is not a polygon"); //ERROR
	}

	return 1;
}

void helpMessage(){
	printf("\nOptions:		Description (default)\n");
	printf(" -i arg			Input shape file\n");
	printf(" -odir arg		Output directory\n");
	printf(" -distance		Generalization distance\n");
	printf(" -angle			Generalization angle for joining as parallel\n");
	printf(" -parallel		Force generalized lines to be parallel\n");
	printf(" -rightAngles		Force generalized lines to have right angles\n");
	printf(" -overwrite		Overwrite existing files\n");
	exit(1);
}

bool removeIfFileExists(string filename){
	FILE *istream;
	if ((fopen_s(&istream, filename.c_str(),"r")) == NULL){
		fclose(istream);
		if(remove(filename.c_str())){
			printf("Error deleting %s",filename.c_str());
			exit(1);
		}
	}
	return true;
}

bool removeExistingFiles(string outputDirectory){
	printf("Overwriting existing files\n");
	
	if(write_adjacent_buildings){
		removeIfFileExists(outputDirectory + "out_AdjacentBuildings.shp");
		removeIfFileExists(outputDirectory + "out_AdjacentBuildings.shx");
		removeIfFileExists(outputDirectory + "out_AdjacentBuildings.dbf");
	}
	if(write_gen_buf_lines){
		removeIfFileExists(outputDirectory + "out_generalizedbufferlines.shp");
		removeIfFileExists(outputDirectory + "out_generalizedbufferlines.shx");
		removeIfFileExists(outputDirectory + "out_generalizedbufferlines.dbf");
	}
	if(write_split_bb){
		removeIfFileExists(outputDirectory + "out_splitbbox.shp");
		removeIfFileExists(outputDirectory + "out_splitbbox.shx");
		removeIfFileExists(outputDirectory + "out_splitbbox.dbf");
	}
	if(write_gen_polygon){
		removeIfFileExists(outputDirectory + "out_Generalized_Polygon.shp");
		removeIfFileExists(outputDirectory + "out_Generalized_Polygon.shx");
		removeIfFileExists(outputDirectory + "out_Generalized_Polygon.dbf");
	}
	if(write_decomp_poly){
		removeIfFileExists(outputDirectory + "out_Decomposed_Polygon.shp");
		removeIfFileExists(outputDirectory + "out_Decomposed_Polygon.shx");
		removeIfFileExists(outputDirectory + "out_Decomposed_Polygon.dbf");
	}

	return false;
}

int main(int argc, char* argv[]){
	string filename;
	//string outputFilename = "out_Generalized_Polygon.shp";
	string outputDirectory = "";
	bool overwrite_files = false;
	int iPolygon = 0;
	int failedPolygons = 0;
	double epsilon = 1;				//generalization distance in meters
	double alpha = 10;				//generalization angle in degrees
	bool makeParallel = false;		//force generalized lines to be parallel
	bool makeRightAngled = false;	//force generalized lines to have right angles

	if(argc == 1) helpMessage();
	
	for(int a = 1; a < argc; a++){
		if(strcmp(argv[a],"-i") == 0){
			a++;
			filename = string(argv[a]);
		}
		else if(strcmp(argv[a],"-odir") == 0){
			a++;
			outputDirectory = string(argv[a]);
		}
		else if(strcmp(argv[a],"-distance") == 0 || strcmp(argv[a],"-d") == 0){
			a++;
			epsilon = atof(argv[a]);
		}
		else if(strcmp(argv[a],"-angle") == 0 || strcmp(argv[a],"-a") == 0){
			a++;
			alpha = atof(argv[a]);
		}
		else if(strcmp(argv[a],"-parallel") == 0){
			makeParallel = true;
		}
		else if(strcmp(argv[a],"-right") == 0 || strcmp(argv[a],"-rightAngles") == 0){
			makeRightAngled = true;
		}
		else if(strcmp(argv[a],"-overwrite") == 0){
			overwrite_files = true;
		}
		else helpMessage();
	}

	printf("Using:\nDistance %.2lf meter \nAngle %.2lf degree \n",epsilon,alpha);
	
	PROCESS_MEMORY_COUNTERS pmc;

	clock_t timer_begin_polygon,timer_total;
	timer_total = clock();
	RegisterOGRShape();

	SHPWriter *bufferLineWriter = new SHPWriter;
	SHPWriter *splitBBoxWriter = new SHPWriter;
	SHPWriter *decompPolyWriter = new SHPWriter;
	SHPWriter *generalizedPolygonWriter = new SHPWriter;
	SHPWriter *adjacentBuildingsWriter = new SHPWriter;

	//Create directory if not existing
	if(!outputDirectory.empty()){
		//Check if directory ends with a '\'
		if ( outputDirectory[ strlen( outputDirectory.c_str() ) - 1 ] != '\\' ){
			outputDirectory += '\\';
		}

		struct stat results;
		if(stat(outputDirectory.c_str(),&results) != 0){
			//create dir
			printf("Output directory does not exist, creating %s\n", outputDirectory.c_str());
			_mkdir(outputDirectory.c_str());
		}
	}

	if(overwrite_files) removeExistingFiles(outputDirectory);

	if(write_adjacent_buildings){
		adjacentBuildingsWriter->open(outputDirectory + "out_AdjacentBuildings.shp",wkbPolygon);
		adjacentBuildingsWriter->addField("ID",OFTInteger);
	}

	if(write_gen_buf_lines){
		bufferLineWriter->open(outputDirectory + "out_generalizedbufferlines.shp",wkbLineString);
		bufferLineWriter->addField("ID",OFTInteger);
		bufferLineWriter->addField("Line",OFTInteger);
	}
	if(write_split_bb){
		splitBBoxWriter->open(outputDirectory + "out_splitbbox.shp",wkbMultiPolygon);
		splitBBoxWriter->addField("ID",OFTInteger);
	}
	if(write_gen_polygon){
		generalizedPolygonWriter->open(outputDirectory + "out_Generalized_Polygon.shp",wkbPolygon);
		generalizedPolygonWriter->addField("FeaID",OFTInteger);
	}
	if(write_decomp_poly){
		decompPolyWriter->open(outputDirectory + "out_Decomposed_Polygon.shp",wkbMultiPolygon);
		decompPolyWriter->addField("FeaID",OFTInteger);
		decompPolyWriter->addField("Points_Org",OFTInteger);
		decompPolyWriter->addField("Points_New",OFTInteger);
		decompPolyWriter->addField("Area_Org",OFTReal);
		decompPolyWriter->addField("Area_New",OFTReal);
		decompPolyWriter->addField("Sym_dif_m",OFTReal);
		decompPolyWriter->addField("Sym_dif_p",OFTReal);
		decompPolyWriter->addField("Area_dif_m",OFTReal);
		decompPolyWriter->addField("Area_dif_p",OFTReal);
		decompPolyWriter->addField("Cont_Org",OFTReal);
		decompPolyWriter->addField("Cont_New",OFTReal);
		decompPolyWriter->addField("Cont_true",OFTReal);
		decompPolyWriter->addField("Haus_dist",OFTReal);
		decompPolyWriter->addField("Correct",OFTInteger);
		decompPolyWriter->addField("Proc_time",OFTReal);
	}

	SHPReader *reader = new SHPReader(filename);

	////////////////////START iteration over all input polygons
	long startPoly	= 0;
	long nb; 
	long endPoly	= reader->Layer->GetFeatureCount();
	for(iPolygon=startPoly;iPolygon<endPoly;iPolygon++){
	//for(iPolygon=0;iPolygon<reader->Layer->GetFeatureCount();iPolygon++){
		try{
			timer_begin_polygon = clock();
			nb = reader->Layer->GetFeature(iPolygon)->GetFieldAsInteger("FeaID");

			printf("\n============ Polygon nr: %d ============\n",nb);
			reader->readPolygon(iPolygon);
			reader->calculateBoundingBox();
			reader->calculateBoundingBoxDouble();
			reader->findAdjacentAdjacentBuildings();


			if(write_adjacent_buildings){
				for(int building=0;building<reader->adjacentBuildings->getNumGeometries();building++){
					OGRFeature *adjacentBuildingFeature = OGRFeature::CreateFeature(adjacentBuildingsWriter->Layer->GetLayerDefn());
					adjacentBuildingFeature->SetGeometry(reader->adjacentBuildings->getGeometryRef(building));
					adjacentBuildingFeature->SetField("ID", iPolygon);
					adjacentBuildingsWriter->writeFeature(adjacentBuildingFeature);
				}
			}

			//Create vector of buffers
			vector<Buffer> generalizedBuffers;
			vector<Buffer> adjacentBuffers;

			//Populate buffers
			//populateBuffers(reader->Polygon,generalizedBuffers);
			//populateBuffersWithInterior(reader->Polygon,generalizedBuffers);
			populateBuffersWithInterior(reader->Geometry,generalizedBuffers);
			populateBuffersWithInterior(reader->adjacentBuildings,adjacentBuffers);
			
			//TODO: join perfect collinear lines into single line
			//joinCollinearBuffers(generalizedBuffers,epsilon,alpha);
			//if(!adjacentBuffers.empty()){
			//	joinCollinearBuffers(adjacentBuffers,epsilon,alpha);
			//}

			//Generalize buffers given epsilon and aplha
			//generalizeBuffers2(generalizedBuffers,epsilon,alpha);
			//generalizeBuffers3(generalizedBuffers,epsilon,alpha);
			generalizeBuffersWithAdjacentBuildings(generalizedBuffers,adjacentBuffers,epsilon,alpha);
			
			removeCollapsedLines(generalizedBuffers,epsilon);

			if(generalizedBuffers.size()<=2) throw "Not enough Buffers after generalization";

			if(makeParallel){
				parallelizeBuffers(generalizedBuffers,10);
			}
			if(makeRightAngled){
				rightAnglesBuffers(generalizedBuffers,10);
			}
			//join parallel buffers that are within epsilon
			//int index=0;
			//while(index<(int)generalizedBuffers.size()-1){
			//	for(vector<Buffer>::iterator buffer=generalizedBuffers.begin()+index+1;buffer!=generalizedBuffers.end();){
			//		if( abs(generalizedBuffers.at(index).B - buffer->B) * cos(atan(generalizedBuffers.at(index).A)) < epsilon){
			//		//if(joinableBuffers(generalizedBuffers.at(index),*buffer,epsilon,alpha)){
			//			joinBuffers(generalizedBuffers.at(index),*buffer);
			//			buffer = generalizedBuffers.erase(buffer);
			//		}else buffer++;
			//	}
			//	index++;
			//}
			//CalculateNewB(generalizedBuffers);

			if(write_gen_buf_lines){
				//Write bufferlines
				//Write generalized line by buffer Bmin
				for(vector<Buffer>::iterator buffer = generalizedBuffers.begin();buffer!=generalizedBuffers.end();++buffer){
					writeBufferLine(buffer->A, buffer->Bmin, 
						reader->BoundingBoxDouble->getExteriorRing()->getX(0), reader->BoundingBoxDouble->getExteriorRing()->getX(2), 
							reader->BoundingBoxDouble, bufferLineWriter->Layer, iPolygon, buffer->id);
				}
				//Write generalized line by buffer Bmax
				for(vector<Buffer>::iterator buffer = generalizedBuffers.begin();buffer!=generalizedBuffers.end();++buffer){
					writeBufferLine(buffer->A, buffer->Bmax, 
						reader->BoundingBoxDouble->getExteriorRing()->getX(0), reader->BoundingBoxDouble->getExteriorRing()->getX(2), 
							reader->BoundingBoxDouble, bufferLineWriter->Layer, iPolygon, buffer->id);
				}
				bufferLineWriter->Layer->SyncToDisk();
			}

			//Cut the double Bounding Box using the bufferlines
			//OGRMultiPolygon* splitDoubleBBox = cutBoundingBox(reader->BoundingBoxDouble,generalizedBuffers);
			//TODO: Find base direction of building > longest line could be wrong direction
			OGRMultiPolygon* splitDoubleBBox = cutBoundingBoxSlantLast(reader->BoundingBoxDouble,generalizedBuffers,epsilon);
			

			//Write splitted double Bounding Box
			if(write_split_bb){
				OGRFeature *splitBBoxFeature = OGRFeature::CreateFeature(splitBBoxWriter->Layer->GetLayerDefn());
				splitBBoxFeature->SetGeometry(splitDoubleBBox);
				splitBBoxFeature->SetField("ID", iPolygon);
				splitBBoxWriter->writeFeature(splitBBoxFeature);
				splitBBoxWriter->Layer->SyncToDisk();
			}

			//Removing unused cells using overlap test of splitted double Bounding Box with input polygon
			overlapTestSplitBBandInput(splitDoubleBBox,reader->Geometry);

			//Test if the overlap test removed all parts
			if(splitDoubleBBox->IsEmpty()) throw "Decomposed polygon IsEmpty";

			//Test if output geometry is valid, MultiPolygon is marked self-intersection so test all parts
			for(int part=0;part<splitDoubleBBox->getNumGeometries();part++){
				if(!GEOSisValid(splitDoubleBBox->getGeometryRef(part)->exportToGEOS())){
					throw GEOSisValidReason(splitDoubleBBox->getGeometryRef(part)->exportToGEOS());
				}
			}

			//test UnionCascaded
			if(!GEOSisValid(splitDoubleBBox->UnionCascaded()->exportToGEOS())){
				printf("UnionCascadedInvalid: %s\n",GEOSisValidReason(splitDoubleBBox->UnionCascaded()->exportToGEOS()));
			}
			
			//Write decomposed polygon created by the overlap test
			if(write_decomp_poly){
				//Create generalized polygon for statistics
				OGRPolygon* genPolygon = (OGRPolygon*) splitDoubleBBox->UnionCascaded()->Simplify(0.000001);
				if( genPolygon==NULL || genPolygon->IsEmpty()) throw "Generalized polygon IsEmpty";

				OGRFeature *decompPolyFeature = OGRFeature::CreateFeature(decompPolyWriter->Layer->GetLayerDefn());
				decompPolyFeature->SetGeometry(splitDoubleBBox);
				decompPolyFeature->SetField("FeaID", nb);

				decompPolyFeature->SetField("Points_Org",reader->Polygon->getExteriorRing()->getNumPoints());
				decompPolyFeature->SetField("Points_New",genPolygon->getExteriorRing()->getNumPoints());
				decompPolyFeature->SetField("Area_Org", reader->Polygon->get_Area());
				decompPolyFeature->SetField("Area_New", splitDoubleBBox->get_Area());

				double Intrusion = OGRGeometry_Get_Area(reader->Geometry->Difference(splitDoubleBBox->UnionCascaded()));
				double Extrusion = OGRGeometry_Get_Area(splitDoubleBBox->UnionCascaded()->Difference(reader->Geometry));
				double symmetricDifference = OGRGeometry_Get_Area(splitDoubleBBox->UnionCascaded()->SymDifference(reader->Geometry));

				decompPolyFeature->SetField("Sym_dif_m", symmetricDifference);
				decompPolyFeature->SetField("Sym_dif_p", symmetricDifference/OGRGeometry_Get_Area(reader->Geometry));
				decompPolyFeature->SetField("Area_dif_m", (Intrusion-Extrusion));
				decompPolyFeature->SetField("Area_dif_p", (Intrusion-Extrusion)/OGRGeometry_Get_Area(reader->Geometry));

				//Calculate contour trueness (Not completely working, still has small overlap, thus better results, due to the buffer use)
				OGRMultiLineString* OrgBoundary = (OGRMultiLineString*) reader->Geometry->Boundary();
				OGRMultiLineString* GenBoundary = (OGRMultiLineString*) splitDoubleBBox->UnionCascaded()->Boundary();
				OGRMultiLineString* IntBoundary = (OGRMultiLineString*) GenBoundary->Intersection(OrgBoundary->Buffer(0.01,1));
				
				decompPolyFeature->SetField("Cont_Org",OrgBoundary->get_Length());
				decompPolyFeature->SetField("Cont_New",GenBoundary->get_Length());
				decompPolyFeature->SetField("Cont_true",IntBoundary->get_Length()/OrgBoundary->get_Length());

				//Calculate Hausdorff Distance
				double hausdorff_distance = 0;
				GEOSHausdorffDistance(reader->Geometry->exportToGEOS(),genPolygon->exportToGEOS(),&hausdorff_distance);
				decompPolyFeature->SetField("Haus_dist",hausdorff_distance);


				//Calculate quasi-euclidean distance of values
				if((Intrusion-Extrusion)/OGRGeometry_Get_Area(reader->Geometry) < 1 //< 1 Percent Area difference
						&& hausdorff_distance < epsilon) {
					decompPolyFeature->SetField("Correct",1);
				}else decompPolyFeature->SetField("Correct",0);

				decompPolyFeature->SetField("Proc_time",(double)(clock()-timer_begin_polygon)/CLOCKS_PER_SEC);
				decompPolyWriter->writeFeature(decompPolyFeature);
				//decompPolyWriter->Layer->SyncToDisk();
			}
			
			//Writing generalized polygon, created with UnionCascaded of decomposed polygon
			if(write_gen_polygon){
				//Convert from multipolygon to polygon and write
				OGRPolygon* temppoly = (OGRPolygon*) splitDoubleBBox->UnionCascaded()->Simplify(0.000001);
				if( temppoly==NULL || temppoly->IsEmpty()) throw "Generalized polygon IsEmpty";
				OGRFeature *generalizedPolygonFeature = OGRFeature::CreateFeature(generalizedPolygonWriter->Layer->GetLayerDefn());
				generalizedPolygonFeature->SetGeometry(temppoly);
				generalizedPolygonFeature->SetField("FeaID", nb);
				//OGRGeometry* symdif = splitDoubleBBox->UnionCascaded()->SymmetricDifference(reader->Geometry);

				////if is multipoly or poly
				//OGRMultiPolygon* symdifMultiPoly = (OGRMultiPolygon*) symdif;
				//double HD = 0;
				//for(int geom = 0;geom<symdifMultiPoly->getNumGeometries();++geom){
				//	if(OGRGeometry_Get_Area(symdifMultiPoly->getGeometryRef(geom)) > 0.1){
				//		
				//		OGRGeometry* symdifINTorg = reader->Geometry->Boundary()->Buffer(0.000001,1)->Intersection(symdifMultiPoly->getGeometryRef(geom));
				//		OGRGeometry* symdifINTgen = splitDoubleBBox->UnionCascaded()->Boundary()->Buffer(0.000001,1)->Intersection(symdifMultiPoly->getGeometryRef(geom));
				//
				//		printf("symdif is %s\n",symdifMultiPoly->getGeometryRef(geom)->getGeometryName());
				//		symdifMultiPoly->getGeometryRef(geom)->dumpReadable(NULL);
				//		printf("symdifINTorg is %s\n",symdifINTorg->getGeometryName());
				//		symdifINTorg->dumpReadable(NULL);
				//		printf("symdifINTgen is %s\n",symdifINTgen->getGeometryName());
				//		symdifINTgen->dumpReadable(NULL);
				//		if(wkbFlatten(symdifINTorg->getGeometryType())!=wkbPoint){

				//			double HDtemp=0;
				//			GEOSHausdorffDistance(symdifINTorg->exportToGEOS(),symdifMultiPoly->getGeometryRef(geom)->exportToGEOS(),&HDtemp);
				//			if(HDtemp>HD) HD=HDtemp;
				//			printf("HD is %lf \n",HD);
				//			GEOSHausdorffDistance(symdifINTgen->exportToGEOS(),symdifMultiPoly->getGeometryRef(geom)->exportToGEOS(),&HDtemp);
				//			if(HDtemp>HD) HD=HDtemp;
				//		}
				//		printf("HD is %lf \n",HD);
				//	}
				//}
				//printf("New Hausdorff distance is %lf \n",HD);

				generalizedPolygonWriter->writeFeature(generalizedPolygonFeature);
			}
			printf("Polygon processing time: %g seconds\n", (double)(clock()-timer_begin_polygon)/CLOCKS_PER_SEC);
			GetProcessMemoryInfo(GetCurrentProcess(), &pmc,sizeof(pmc));
			printf("Memory footprint: %5lu MB \n", pmc.WorkingSetSize/1048576);

		}catch(char* err){
			printf("Building %d processing failure: %s\n",iPolygon,err);
			failedPolygons++;
		}
		if(iPolygon%50==0){
			if(write_gen_buf_lines) bufferLineWriter->Layer->SyncToDisk();
			if(write_split_bb) splitBBoxWriter->Layer->SyncToDisk();
			if(write_decomp_poly) decompPolyWriter->Layer->SyncToDisk();
			if(write_gen_polygon) generalizedPolygonWriter->Layer->SyncToDisk();
		}
	}
	////////////////////END iteration over all input polygons
	
	if(write_gen_buf_lines) bufferLineWriter->close();
	if(write_split_bb) splitBBoxWriter->close();
	if(write_decomp_poly) decompPolyWriter->close();
	if(write_gen_polygon) generalizedPolygonWriter->close();
	if(write_adjacent_buildings) adjacentBuildingsWriter->close();
	
	delete bufferLineWriter;
	delete splitBBoxWriter;
	delete decompPolyWriter;
	delete generalizedPolygonWriter;
	delete adjacentBuildingsWriter;

	reader->close();
	delete reader;
	printf("\n\nTotal number of polygons: %d \n", iPolygon-1);
	printf("Total number of failed polygons: %d \n", failedPolygons);
	printf("Total processing time: %g seconds\n", (double)(clock()-timer_total)/CLOCKS_PER_SEC);
	printf("Average processing time per building: %g seconds\n", ((double)(clock()-timer_total)/CLOCKS_PER_SEC)/(iPolygon-1));
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc,sizeof(pmc));
	printf("Maximum memory footprint: %5lu MB \n", pmc.PeakWorkingSetSize/1048576);
}