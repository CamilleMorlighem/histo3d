//#include <iostream>
//#include <sstream>
//#include <string>
//#include <vector>
//#include <algorithm>
////#include "ogr\ogrsf_frmts.h"
//#include "ogr\ogr_core.h"
//#include "Cell_Decomposition.h"
//
//using namespace std;
//
//const double PI = 4*atan(1.0);
//
//typedef struct Line{
//	int id;
//	double A,B,C,Cmin,Cmax;
//	double x[2],y[2],bearing,length;
//	bool marked;
//	OGRLineString line;
//	OGRPoint *intersection;
//	vector<OGRPoint> point;
//} Line;
//
//bool pairSorter(pair<int,double> &i,pair<int,double> &j){
//	return i.second > j.second;
//};
//
//int vecSorter(vector<int> &i,vector<int> &j){
//	return i.at(0) < j.at(0);
//};
//
//int CalculateIntersections(OGRPoint one,OGRPoint two,OGRLinearRing bb,OGRPoint *outPoints){
//	OGRPoint three,four;
//	int linenr=0;
//	double ua,ub;
//
//	for(int i=0;i<4;i++){
//		bb.getPoint(i,&three);
//		bb.getPoint(i+1,&four);
//		ua = ((four.getX()-three.getX())*(one.getY()-three.getY()) - (four.getY()-three.getY())*(one.getX()-three.getX())) 
//				/ ((four.getY()-three.getY())*(two.getX()-one.getX()) - (four.getX()-three.getX())*(two.getY()-one.getY()));
//		ub = ((two.getX()-one.getX())   *(one.getY()-three.getY()) - (two.getY()-one.getY())*(one.getX()-three.getX()))
//				/ ((four.getY()-three.getY())*(two.getX()-one.getX()) - (four.getX()-three.getX())*(two.getY()-one.getY()));
//	
//		if((ua>0 && ua<1)||(ub>0 && ub<1)){
//			outPoints[linenr].setX(one.getX() + ua * (two.getX() - one.getX()));
//			outPoints[linenr].setY(one.getY() + ua * (two.getY() - one.getY()));
//			linenr++;
//		}
//	}
//	return 1;
//}
//
//double minimumC(Line &line1, Line &line2){
//	double tempC=0,minC=999999999;
//
//	//minC=line1.Cmin;
//	for(vector<OGRPoint>::iterator it = line2.point.begin(); it != line2.point.end(); ++it){
//		tempC = -line1.A * it->getX() - line1.B * it->getY();
//		if(tempC<minC) minC=tempC;
//	}
//		return minC;
//}
//
//double maximumC(Line &line1, Line &line2){
//	double tempC=0,maxC=-999999999;
//	
//	//maxC=line1.Cmax;
//	for(vector<OGRPoint>::iterator it = line2.point.begin(); it != line2.point.end(); ++it){
//		tempC = -line1.A * it->getX() - line1.B * it->getY();
//		if(tempC>maxC) maxC=tempC;
//	}
//		return maxC;
//}
//
//bool joinable(Line &line1, Line &line2, double epsilon, double alpha){
//	double angle;
//	double tempC=0,minC,maxC;
//	alpha = alpha / (180/PI);
//	angle = acos( (line1.A*line2.A + line1.B*line2.B) / ( sqrt(line1.A*line1.A + line1.B*line1.B) * sqrt(line2.A*line2.A + line2.B*line2.B) ) );
//
//	minC=line1.Cmin;
//	maxC=line1.Cmax;
//	for(vector<OGRPoint>::iterator it = line2.point.begin(); it != line2.point.end(); ++it){
//		tempC = -(line1.A * it->getX()) - (line1.B * it->getY());
//		if(tempC<minC) minC=tempC;
//		if(tempC>maxC) maxC=tempC;
//	}
//
//	if( (angle < alpha) && (abs(((maxC-minC) / -line1.B) * cos(atan(line1.A/-line1.B))) < epsilon) ){
//		return true;
//	}
//	return false;
//}
//
//int joinLinesAB(Line &line1, Line &line2){
//	double weight1,weight2,weightangle;
//
//	weightangle = cos(atan( (line2.A/-line2.B) - ( line1.A/-line1.B) ));
//	weight1 = line1.length*weightangle		/ (line1.length*weightangle + abs(line2.length*weightangle));
//	weight2 = abs(line2.length*weightangle) / (line1.length*weightangle + abs(line2.length*weightangle));
//	line1.A = line1.A * weight1 + line2.A * weight2;
//	line1.B = line1.B * weight1 + line2.B * weight2;
//	//line1.length = line1.length + line2.length*weightangle;
//
//	return 1;
//}
//
//int joinLinesC(Line &line1, Line &line2){
//	double weight1,weight2,weightangle;
//
//	//For the weight the total length is taken because this is updated for the AB updates
//	weightangle = cos(atan( (line2.A/-line2.B) - ( line1.A/-line1.B) ));
//	weight1 = (line1.length-abs(line2.length*weightangle))	/ line1.length;
//	weight2 = abs(line2.length*weightangle)					/ line1.length;
//	line1.Cmin = line1.Cmin * weight1 + ((maximumC(line1,line2) + minimumC(line1,line2))/2) * weight2;
//
//	return 1;
//}
//
//bool includable(Line &line1, Line &line2, double epsilon){
//	double tempC=0,minC,maxC;
//	minC=line1.Cmin;
//	maxC=line1.Cmax;
//	for(vector<OGRPoint>::iterator it = line2.point.begin(); it != line2.point.end(); ++it){
//		tempC = -(line1.A * it->getX()) - (line1.B * it->getY());
//		if(tempC<minC) minC=tempC;
//		if(tempC>maxC) maxC=tempC;
//	}
//
//	if( abs(((maxC-minC) / -line1.B) * cos(atan(line1.A/-line1.B))) < epsilon ){
//		return true;
//	}
//	return false;
//}
//
//bool includableCreateBuffer(Line &line1, Line &line2, double epsilon){
//	double tempC=0,minC=0,maxC=0;
//	double x[4],y[4];
//	OGRLineString *bufferline = new OGRLineString[4];
//	int i=0;
//	minC=line1.Cmin;
//	maxC=line1.Cmax;
//	x[0] = line1.point[0].getX()-10; 
//	y[0] = (line1.A* (line1.point[0].getX()-10) + minC) / (-line1.B);
//	x[1] = line1.point[0].getX()+10; 
//	y[1] = (line1.A* (line1.point[0].getX()+10) + minC) / (-line1.B);
//	x[2] = line1.point[0].getX()-10; 
//	y[2] = (line1.A* (line1.point[0].getX()-10) + maxC) / (-line1.B);
//	x[3] = line1.point[0].getX()+10; 
//	y[3] = (line1.A* (line1.point[0].getX()+10) + maxC) / (-line1.B);
//
//	for(vector<OGRPoint>::iterator it = line2.point.begin(); it != line2.point.end(); ++it){
//		tempC = -(line1.A * it->getX()) - (line1.B * it->getY());
//		if(tempC<minC){
//			minC=tempC;
//			x[0] = it->getX()-10; 
//			y[0] = (line1.A* (it->getX()-10) + minC) / (-line1.B);
//			x[1] = it->getX()+10; 
//			y[1] = (line1.A* (it->getX()+10) + minC) / (-line1.B);
//		}
//		if(tempC>maxC){
//			maxC=tempC;
//			x[2] = it->getX()-10; 
//			y[2] = (line1.A* (it->getX()-10) + maxC) / (-line1.B);
//			x[3] = it->getX()+10; 
//			y[3] = (line1.A* (it->getX()+10) + maxC) / (-line1.B);
//		}
//	}
//
//	bufferline[0].addPoint(x[0],y[0]);
//	bufferline[0].addPoint(x[1],y[1]);
//	bufferline[1].addPoint(x[2],y[2]);
//	bufferline[1].addPoint(x[3],y[3]);
//
//	OGRSFDriver *poDriverOut;
//	OGRDataSource *poDSOut;
//	OGRLayer *poLayerOut;
//	OGRFeature *poFeatureOut;
//	OGRFieldDefn oField("ID", OFTInteger);
//	oField.SetWidth(5);
//	try{
//		poDriverOut = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
//		poDSOut = poDriverOut->CreateDataSource("out_bufferlines.shp", NULL);
//		poLayerOut = poDSOut->CreateLayer( "line", NULL, wkbLineString, NULL );
//		poLayerOut->CreateField( &oField );
//		poFeatureOut = OGRFeature::CreateFeature(poLayerOut->GetLayerDefn());
//		poFeatureOut->SetField("ID", 1);
//		poFeatureOut->SetGeometry(&bufferline[0]);
//		poLayerOut->CreateFeature(poFeatureOut);
//		poFeatureOut->SetGeometry(&bufferline[1]);
//		poLayerOut->CreateFeature(poFeatureOut);
//	}catch(char *err){
//		printf("Exception raised: %s\n",err);
//	}
//
//	OGRFeature::DestroyFeature(poFeatureOut);
//	OGRDataSource::DestroyDataSource(poDSOut);
//	
//	cout << "abs(((maxC-minC) / -line1.B) * cos(atan(line1.A/-line1.B))) " << abs(((maxC-minC) / -line1.B) * cos(atan(line1.A/-line1.B)))<< endl;
//
//	if( abs(((maxC-minC) / -line1.B) * cos(atan(line1.A/-line1.B))) < epsilon ){
//		return true;
//	}
//	return false;
//}
//
//bool CreateBufferLine(double A, double B, double C, double xIn1, double xIn2, OGRLayer *poLayerOut, int id){
//	double x[2],y[2];
//	OGRLineString bufferline;
//	int i=0;
//
//	x[0] = xIn1; 
//	y[0] = (A* xIn1 + C) / (-B);
//	x[1] = xIn2; 
//	y[1] = (A* xIn2 + C) / (-B);
//
//	bufferline.addPoint(x[0],y[0]);
//	bufferline.addPoint(x[1],y[1]);
//	
//	OGRFeature *poFeatureOut;
//		
//	try{
//		poFeatureOut = OGRFeature::CreateFeature(poLayerOut->GetLayerDefn());
//		poFeatureOut->SetField("ID", id);
//		poFeatureOut->SetGeometry(&bufferline);
//		poLayerOut->CreateFeature(poFeatureOut);
//	}catch(char *err){
//		printf("Exception raised: %s\n",err);
//	}
//
//	OGRFeature::DestroyFeature(poFeatureOut);
//
//	return true;
//}
//
//bool CreateBufferLine(double A, double B, double C, double xIn1, double xIn2, OGRLayer *poLayerOut, int id, double weight){
//	double x[2],y[2];
//	OGRLineString bufferline;
//	int i=0;
//
//	x[0] = xIn1; 
//	y[0] = (A* xIn1 + C) / (-B);
//	x[1] = xIn2; 
//	y[1] = (A* xIn2 + C) / (-B);
//
//	bufferline.addPoint(x[0],y[0]);
//	bufferline.addPoint(x[1],y[1]);
//	
//	OGRFeature *poFeatureOut;
//		
//	try{
//		poFeatureOut = OGRFeature::CreateFeature(poLayerOut->GetLayerDefn());
//		poFeatureOut->SetField("ID", id);
//		poFeatureOut->SetField("Weight", weight);
//		poFeatureOut->SetGeometry(&bufferline);
//		poLayerOut->CreateFeature(poFeatureOut);
//	}catch(char *err){
//		printf("Exception raised: %s\n",err);
//	}
//
//	OGRFeature::DestroyFeature(poFeatureOut);
//
//	return true;
//}
//
//int CalculateIntersectionABC(Line a, Line b, OGRPoint *point){
//	double det = a.A*b.B - b.A*a.B;
//	double x,y;
//	if(det == 0){
//		cout << "Lines are parallel" << endl;
//		return 0;
//	}else{
//		x = -(b.B*a.Cmin - a.B*b.Cmin)/det;
//		y = -(a.A*b.Cmin - b.A*a.Cmin)/det;
//	}
//	point->setX(x);
//	point->setY(y);
//	//printf("det: %lf\n ",det);
//	//printf("x,y = %10lf,%10lf\n",x,y);
//
//	return 1;
//}
//
//int main(int argc, char* argv[]){
//	string filename;
//	int iPolygon = 0;
//	double bb[4],*x,*y;
//	Line *lines;
//
//	for(int a = 1; a < argc; a++){
//		if(strcmp(argv[a],"-i") == 0){
//			a++;
//			filename = string(argv[a]);
//		}
//		else if(strcmp(argv[a],"-p") == 0){
//			a++;
//			iPolygon = atoi(argv[a]);
//		}
//	}
//
//	OGRRegisterAll();
//
//	//Reading input file as polygons
//	OGRDataSource *poDS;
//	OGRLayer  *poLayer;
//	OGRFeature *poFeature;
//	OGRGeometry *poGeometry;
//	OGRPolygon *poPolygon;
//	OGREnvelope *psEnvelope = new OGREnvelope;
//	double org_area;
//
//	poDS = OGRSFDriverRegistrar::Open(filename.c_str(), FALSE);
//    if(poDS == NULL){
//        printf("Open failed.\n");
//        exit(1);
//    }
//	  
//    poLayer = poDS->GetLayer(0);
//	printf("Feature count is: %d\n", poLayer->GetFeatureCount());
//	
//	poFeature = poLayer->GetFeature(iPolygon);
//	poGeometry = poFeature->GetGeometryRef();
//
//	if( poGeometry != NULL 
//            && wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon ){
//        poPolygon = (OGRPolygon *) poGeometry;
//		poPolygon->getEnvelope(psEnvelope);
//		printf( "Bounding box:    %lf,%lf\n                 %lf,%lf\n",
//			psEnvelope->MinX,psEnvelope->MinY,psEnvelope->MaxX,psEnvelope->MaxY);
//		org_area = poPolygon->get_Area();
//    }else{
//        printf( "No polygon geometry\n" );
//    }       
//
//	bb[0] = psEnvelope->MinX-((psEnvelope->MaxX-psEnvelope->MinX)/2);
//	bb[1] = psEnvelope->MinY-((psEnvelope->MaxY-psEnvelope->MinY)/2);
//	bb[2] = psEnvelope->MaxX+((psEnvelope->MaxX-psEnvelope->MinX)/2);
//	bb[3] = psEnvelope->MaxY+((psEnvelope->MaxY-psEnvelope->MinY)/2);
//	x = new double[5];
//	y = new double[5];
//	x[0]=bb[0];
//	x[1]=bb[0];
//	x[2]=bb[2];
//	x[3]=bb[2];
//	x[4]=bb[0];
//	y[0]=bb[1];
//	y[1]=bb[3];
//	y[2]=bb[3];
//	y[3]=bb[1];
//	y[4]=bb[1];
//
//	printf( "\nBounding box 2x: %lf,%lf\n                 %lf,%lf\n",
//		bb[0],bb[1],bb[2],bb[3]);
//
//	OGRPolygon boundingbox;
//	OGRLinearRing ring;
//	ring.setPoints(5,x,y,NULL);
//	boundingbox.addRing(&ring);
//
//
//	//writing the 2x bounding box as polygon to shp
//	OGRSFDriver *poDriverBBOut;
//	OGRDataSource *poDSBBOut;
//	OGRLayer *poLayerBBOut;
//	OGRFeature *poFeatureBBOut;
//	OGRFieldDefn oBBField("ID", OFTInteger);
//	oBBField.SetWidth(5);
//
//	bool write_bb = false;
//	bool write_ip = false;
//	bool write_sp = false;
//	
//	if(write_bb){
//		try{
//			poDriverBBOut = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
//			poDSBBOut = poDriverBBOut->CreateDataSource("out_boundingbox.shp", NULL);
//			poLayerBBOut = poDSBBOut->CreateLayer( "polygon", NULL, wkbPolygon, NULL );
//			poLayerBBOut->CreateField( &oBBField );
//			poFeatureBBOut = OGRFeature::CreateFeature(poLayerBBOut->GetLayerDefn());
//			poFeatureBBOut->SetField("ID", 1);
//			poFeatureBBOut->SetGeometry(&boundingbox);
//			poLayerBBOut->CreateFeature(poFeatureBBOut);
//		}catch(char *err){
//			printf("Exception raised: %s\n",err);
//		}
//
//		OGRFeature::DestroyFeature(poFeatureBBOut);
//		OGRDataSource::DestroyDataSource(poDSBBOut);
//	}
//
//	OGRSFDriver *poDriverPointOut;
//	OGRDataSource *poDSPointOut;
//	OGRLayer *poLayerPointOut;
//	OGRFeature *poFeaturePointOut;
//	OGRFieldDefn oPointField("ID", OFTInteger);
//	oPointField.SetWidth(5);
//	
//	if(write_ip){
//		try{
//			poDriverPointOut = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
//			poDSPointOut = poDriverPointOut->CreateDataSource("out_intersectionpoints.shp", NULL);
//			poLayerPointOut = poDSPointOut->CreateLayer("point", NULL, wkbPoint, NULL);
//			poLayerPointOut->CreateField(&oPointField);
//		}catch(char *err){
//			printf("Exception raised: %s\n",err);
//		}
//	}
//
//	//Calculate lines
//	int iLine=0;
//	lines = new Line[poPolygon->getExteriorRing()->getNumPoints()-1];
//	OGRRawPoint *points = new OGRRawPoint[poPolygon->getExteriorRing()->getNumPoints()];
//	poPolygon->getExteriorRing()->getPoints(points,NULL);
//	for(int i=0 ; i<poPolygon->getExteriorRing()->getNumPoints()-1; i++){
//		lines[i].marked = false;
//		lines[i].id = iLine;
//		iLine++;
//		printf("ID: %3d, ", lines[i].id);
//		
//		lines[i].length = sqrt(pow(points[i+1].x-points[i].x,2)
//							+ pow(points[i+1].y-points[i].y,2));
//		printf("Length is %10lf, ", lines[i].length);
//
//		lines[i].bearing = atan ( (points[i+1].y-points[i].y) 
//								/ (points[i+1].x-points[i].x)) * (180/PI);
//		printf("Bearing is %10lf\n", lines[i].bearing);
//
//		
//		lines[i].line.addPoint(points[i].x,points[i].y);
//		lines[i].line.addPoint(points[i+1].x,points[i+1].y);
//		
//		OGRPoint aPoint;
//		aPoint.setX(points[i].x);
//		aPoint.setY(points[i].y);
//		lines[i].point.push_back(aPoint);
//		aPoint.setX(points[i+1].x);
//		aPoint.setY(points[i+1].y);
//		lines[i].point.push_back(aPoint);
//
//		lines[i].A = points[i].y-points[i+1].y;
//		lines[i].B = points[i+1].x-points[i].x;
//		lines[i].Cmin = points[i].x*points[i+1].y - points[i+1].x*points[i].y;
//		lines[i].Cmax = points[i].x*points[i+1].y - points[i+1].x*points[i].y;
//
//		printf("A is %10lf, B is %10lf, C is %10lf \n", lines[i].A,lines[i].B,lines[i].Cmin);
//		printf("Slope is %10lf, Offset is %10lf \n\n", lines[i].A/-lines[i].B,lines[i].Cmin/-lines[i].B);
//		//printf("test offset: %lf \n",-(lines[i].Cmin/lines[i].B)/bb[1]);
//		//printf("ABC test %lf \n",lines[i].A*points[i].x + lines[i].B*points[i].y + lines[i].Cmin);
//		//printf("x1,y1: %lf,%lf x2,y2: %lf,%lf \n\n", points[i].x,points[i+1].x,points[i].y,points[i+1].y);
//	}
//	
//	//Calculate intersection points
//	for(int i=0;i<sizeof(lines);i++){
//		OGRPoint *intersections = new OGRPoint[2];
//		OGRPoint point0,point1;
//		lines[i].line.getPoint(0,&point0);
//		lines[i].line.getPoint(1,&point1);
//
//		CalculateIntersections(point0,point1,boundingbox.getExteriorRing(),intersections);
//		lines[i].intersection = intersections;
//
//		//Write intersection points to shp
//		if(write_ip){
//			for(int j=0;j<2;j++){
//				try{
//					poFeaturePointOut = OGRFeature::CreateFeature(poLayerPointOut->GetLayerDefn());
//					poFeaturePointOut->SetField("ID", 1);
//					poFeaturePointOut->SetGeometry(&intersections[j]);
//					poLayerPointOut->CreateFeature(poFeaturePointOut);
//				}catch(char *err){
//					printf("Exception raised: %s\n",err);
//				}
//			}
//		}
//	}
//	if(write_ip){
//		OGRFeature::DestroyFeature(poFeaturePointOut);
//		OGRDataSource::DestroyDataSource(poDSPointOut);
//	}
//	
//	OGRSFDriver *poDriverSPOut;
//	OGRDataSource *poDSSPOut;
//	OGRLayer *poLayerSPOut;
//	OGRFeature *poFeatureSPOut;
//	OGRFieldDefn oSPField("ID", OFTInteger);
//	oSPField.SetWidth(5);
//
//	if(write_sp){
//		poDriverSPOut = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
//		poDSSPOut = poDriverSPOut->CreateDataSource("out_simplepolygons.shp", NULL);
//		poLayerSPOut = poDSSPOut->CreateLayer( "polygon", NULL, wkbPolygon, NULL );
//		poLayerSPOut->CreateField( &oSPField );
//	}
//	
//	//Simplification with ogr simplify
//	for(int i=0;i<6;i++){
//		poFeature = poLayer->GetFeature(i);
//		poGeometry = poFeature->GetGeometryRef();
//	
//		if( poGeometry != NULL 
//				&& wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon ){
//			poPolygon = (OGRPolygon *) poGeometry;
//		}else{
//			printf( "No polygon geometry\n" );
//		}
//
//		OGRGeometry *simple = poPolygon->Simplify(0.20);
//		
//		//Writing the simplified polygons to shp
//		if(write_sp){
//			try{
//				poFeatureSPOut = OGRFeature::CreateFeature(poLayerSPOut->GetLayerDefn());
//				poFeatureSPOut->SetField("ID", 1);
//				poFeatureSPOut->SetGeometry(simple);
//				poLayerSPOut->CreateFeature(poFeatureSPOut);
//			}catch(char *err){
//				printf("Exception raised: %s\n",err);
//			}
//		}
//	}
//	if(write_sp){
//		OGRFeature::DestroyFeature(poFeatureSPOut);
//		OGRDataSource::DestroyDataSource(poDSSPOut);
//	}
//
//	//Create half-spaces
//	pair<int,double> *sortedlength = new pair<int,double>[iLine];
//	vector<vector<int>> bucketVec(iLine);
//	int j,numBuckets=0;
//	double epsilon = 1; //buffer in meters
//	double alpha = 5;   //angle difference in degrees
//
//	//Create pairs of line id and length
//	for(int i=0;i<iLine;i++){
//		sortedlength[i]= make_pair(lines[i].id,lines[i].length);
//	}
//
//	//Sort lines on length
//	sort(&sortedlength[0],&sortedlength[(iLine)],pairSorter);
//
//	//Create buckets from lines based on epsilon and alpha
//	for(int i=0;i<iLine;i++){
//		if(!lines[sortedlength[i].first].marked){
//			lines[sortedlength[i].first].marked = true;
//			bucketVec.at(numBuckets).push_back(sortedlength[i].first);
//			if(sortedlength[i].first==iLine-1){
//				j=0;
//			}else{ 
//				j = sortedlength[i].first+1;
//			}
//
//			//while( ( abs(lines[sortedlength[i].first].bearing - lines[j].bearing) < alpha || lines[j].length < epsilon/2 ) 
//			//		&& j!=sortedlength[i].first && !lines[j].marked){
//			while( (joinable(lines[sortedlength[i].first],lines[j],epsilon,alpha) 
//					|| includable(lines[sortedlength[i].first],lines[j],epsilon))
//						&& j!=sortedlength[i].first && !lines[j].marked){
//				printf("pass+: id: %d matches id: %d\n",sortedlength[i].first,lines[j].id);
//				lines[j].marked=true;
//				bucketVec.at(numBuckets).push_back(lines[j].id);
//				//if( joinable(lines[sortedlength[i].first],lines[j],epsilon,alpha) ) printf("Joinable+\n");
//				//if( includable(lines[sortedlength[i].first],lines[j],epsilon) ) printf("Includable+\n");
//				
//				if(j<iLine-1){
//					j++;
//				}else{
//					j=0;
//				}
//			}
//
//			if(sortedlength[i].first==0){
//				j=iLine-1;
//			}else{ 
//				j = sortedlength[i].first-1;
//			}
//
//			//while( (abs(lines[sortedlength[i].first].bearing - lines[j].bearing) < alpha || (lines[j].length < epsilon/2)) 
//			//		&& j!=sortedlength[i].first && !lines[j].marked){
//			while( (joinable(lines[sortedlength[i].first],lines[j],epsilon,alpha) 
//					|| includable(lines[sortedlength[i].first],lines[j],epsilon))
//						&& j!=sortedlength[i].first && !lines[j].marked){
//				printf("pass-: id: %d matches id: %d\n",sortedlength[i].first,lines[j].id);
//				lines[j].marked=true;
//				bucketVec.at(numBuckets).push_back(lines[j].id);
//				//if( joinable(lines[sortedlength[i].first],lines[j],epsilon,alpha) ) printf("Joinable-\n");
//				//if( includable(lines[sortedlength[i].first],lines[j],epsilon) ) printf("Includable-\n");
//
//				if(j>0){
//					j--;
//				}else{
//					j=iLine-1;
//				}
//			}
//			numBuckets++;
//		}
//	}
//
//	//resize the vector to remove empty buckets
//	bucketVec.resize(numBuckets);
//
//	//sort inside buckets by linenumber and outside by first linenumber
//	//for(vector< vector<int> >::iterator ite = bucketVec.begin(); ite != bucketVec.end(); ++ite){
//	//	sort(ite->begin(),ite->end());
//	//}
//	sort(bucketVec.begin(),bucketVec.end(),vecSorter);
//
//	Line *generalizedLines = new Line[bucketVec.size()];
//	
//	//
//	char buffer[10];
//	string outfilename = "out_bufferlines";
//	_itoa_s(iPolygon,buffer,10);
//	outfilename.append(buffer);
//	outfilename.append(".shp");
//	OGRSFDriver *poDriverBufferOut;
//	OGRDataSource *poDSBufferOut;
//	OGRLayer *poLayerBufferOut;
//	poDriverBufferOut = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
//	poDSBufferOut = poDriverBufferOut->CreateDataSource(outfilename.c_str(), NULL);
//	poLayerBufferOut = poDSBufferOut->CreateLayer( "line", NULL, wkbLineString, NULL );
//	OGRFieldDefn oBufferField("ID", OFTInteger);
//	oBufferField.SetWidth(5);
//	poLayerBufferOut->CreateField( &oBufferField );
//	OGRFieldDefn oBufferField2("Weight", OFTReal);
//	oBufferField2.SetWidth(5);
//	oBufferField2.SetPrecision(5);
//	poLayerBufferOut->CreateField( &oBufferField2 );
//	//
//	
//	//Calculate the generalized lines
//	vector< vector<int> >::iterator it = bucketVec.begin();
//	vector< vector<int> >::iterator it_end = bucketVec.end();
//	int i=0;
//	double totalLength=0,weight1=0,weight2=0,weightangle=0;
//	double midC=0;
//	for(it, i; it != it_end; ++it,i++){
//		vector<int>::iterator it2 = it->begin();
//		vector<int>::iterator it2_end = it->end();
//		OGRPoint linePoint;
//		linePoint.setX(lines[*it2].line.getX(0));
//		linePoint.setY(lines[*it2].line.getY(0));
//		generalizedLines[i].point.push_back(linePoint);
//		generalizedLines[i].id = i;
//		generalizedLines[i].A = lines[*it2].A;
//		generalizedLines[i].B = lines[*it2].B;
//		generalizedLines[i].Cmin = lines[*it2].Cmin;
//		generalizedLines[i].Cmax = lines[*it2].Cmax;
//		generalizedLines[i].length = lines[*it2].length;
//
//		//Update A and B
//		for(it2; it2 != it2_end; ++it2){
//			linePoint.setX(lines[*it2].line.getX(1));
//			linePoint.setY(lines[*it2].line.getY(1));
//			generalizedLines[i].point.push_back(linePoint);
//			weightangle = cos(atan( (lines[*it2].A/-lines[*it2].B) - ( generalizedLines[i].A/-generalizedLines[i].B) ));
//
//			if(joinable(generalizedLines[i],lines[*it2],epsilon,alpha)){
//				joinLinesAB(generalizedLines[i],lines[*it2]);
//			}
//
//			if(it2 != it->begin()){
//				generalizedLines[i].length=generalizedLines[i].length + lines[*it2].length*weightangle;
//			}
//		}
//
//		//Update C
//		it2 = it->begin();
//		//generalizedLines[i].length = lines[*it2].length;
//		generalizedLines[i].Cmin = ( maximumC(generalizedLines[i],lines[*it2]) + minimumC(generalizedLines[i],lines[*it2])) /2;
//		for(it2; it2 != it2_end; ++it2){
//			joinLinesC(generalizedLines[i],lines[*it2]);
//
//			//CreateBufferLine(generalizedLines[i].A,generalizedLines[i].B,generalizedLines[i].Cmin,bb[0],bb[2],poLayerBufferOut,*it2,weight2);
//			midC= (minimumC(generalizedLines[i],lines[*it2])+maximumC(generalizedLines[i],lines[*it2])) /2;
//			CreateBufferLine(generalizedLines[i].A,generalizedLines[i].B,midC,bb[0],bb[2],poLayerBufferOut,*it2,weight2);
//			//CreateBufferLine(generalizedLines[i].A,generalizedLines[i].B,minimumC(generalizedLines[i],lines[*it2]),bb[0],bb[2],poLayerBufferOut,*it2);
//			//CreateBufferLine(generalizedLines[i].A,generalizedLines[i].B,maxC(generalizedLines[i],lines[*it2]),bb[0],bb[2],poLayerBufferOut,*it2);
//		}
//		printf("Slope is %10lf, Offset is %10lf \n", generalizedLines[i].A/-(generalizedLines[i].B),generalizedLines[i].Cmin/-(generalizedLines[i].B));
//	}
//	
//	OGRDataSource::DestroyDataSource(poDSBufferOut);
//
//	outfilename="out_intersectionpoints";
//	_itoa_s(iPolygon,buffer,10);
//	outfilename.append(buffer);
//	outfilename.append(".shp");
//	poDriverPointOut = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
//	poDSPointOut = poDriverPointOut->CreateDataSource(outfilename.c_str(), NULL);
//	poLayerPointOut = poDSPointOut->CreateLayer("point", NULL, wkbPoint, NULL);
//	poLayerPointOut->CreateField(&oPointField);
//
//	//Calculate intersection of generalized lines
//	OGRPoint tempoint;
//	for(int i=0; i<numBuckets; i++){
//		if(i==numBuckets-1){
//			CalculateIntersectionABC(generalizedLines[0], generalizedLines[numBuckets-1],&tempoint);
//		}else{
//			CalculateIntersectionABC(generalizedLines[i], generalizedLines[i+1], &tempoint);
//		}
//		if(tempoint.Within(&boundingbox)){
//			try{
//				poFeaturePointOut = OGRFeature::CreateFeature(poLayerPointOut->GetLayerDefn());
//				poFeaturePointOut->SetField("ID", 1);
//				poFeaturePointOut->SetGeometry(&tempoint);
//				poLayerPointOut->CreateFeature(poFeaturePointOut);
//				OGRFeature::DestroyFeature(poFeaturePointOut);
//			}catch(char *err){
//				printf("Exception raised: %s\n",err);
//			}
//		}
//	}
//
//	OGRDataSource::DestroyDataSource(poDSPointOut);
//
//
//	outfilename="out_generalizedpoly";
//	_itoa_s(iPolygon,buffer,10);
//	outfilename.append(buffer);
//	outfilename.append(".shp");
//	OGRSFDriver *poDriverPolyOut;
//	OGRDataSource *poDSPolyOut;
//	OGRLayer *poLayerPolyOut;
//	OGRFeature *poFeaturePolyOut;
//	OGRFieldDefn oPolyField("ID", OFTInteger);
//	oPolyField.SetWidth(5);
//	poDriverPolyOut = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
//	poDSPolyOut = poDriverPolyOut->CreateDataSource(outfilename.c_str(), NULL);
//	poLayerPolyOut = poDSPolyOut->CreateLayer("poly", NULL, wkbPolygon, NULL);
//	poLayerPolyOut->CreateField(&oPolyField);
//
//	//double *xtemp,*ytemp;
//	//xtemp = new OGRPoint[numBuckets];
//	//ytemp = new double[numBuckets];
//	//OGRPoint tempoint;
//	OGRLinearRing linering;
//	OGRPolygon genpoly;
//	for(int i=0; i<numBuckets; i++){
//		if(i==numBuckets-1){
//			CalculateIntersectionABC(generalizedLines[0], generalizedLines[numBuckets-1],&tempoint);
//		}else{
//			CalculateIntersectionABC(generalizedLines[i], generalizedLines[i+1], &tempoint);
//		}
//		if(tempoint.Within(&boundingbox)){
//			linering.addPoint(&tempoint);
//		}
//	}
//	genpoly.addRing(&linering);
//	try{
//		poFeaturePolyOut = OGRFeature::CreateFeature(poLayerPolyOut->GetLayerDefn());
//		poFeaturePolyOut->SetField("ID", 1);
//		poFeaturePolyOut->SetGeometry(&genpoly);
//		poLayerPolyOut->CreateFeature(poFeaturePolyOut);
//		OGRFeature::DestroyFeature(poFeaturePolyOut);
//		OGRDataSource::DestroyDataSource(poDSPolyOut);
//	}catch(char *err){
//		printf("Exception raised: %s\n",err);
//	}
//	
//	printf("Old area; %lf\nNew area: %lf\nDiff area: %lf\n",org_area,genpoly.get_Area(),genpoly.get_Area()-org_area);
//
//	//OGRFeature *FeatureOut2;
//	//SHPWriter *writer = new SHPWriter("blabla.shp",wkbPolygon);
//	//writer->addField("ID",OFTInteger);
//
//	//FeatureOut2 = OGRFeature::CreateFeature(writer->Layer->GetLayerDefn());
//	//FeatureOut2->SetGeometry(&boundingbox);
//	//FeatureOut2->SetField("ID", 1);
//
//	//writer->writeFeature(FeatureOut2);
//	//writer->close();
//
//	return 0;
//}