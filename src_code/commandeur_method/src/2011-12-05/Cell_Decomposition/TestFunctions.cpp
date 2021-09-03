//#include "Cell_Decomposition.h"
//
//using namespace std;
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
////Calculate intersection points
//void WriteOrigionalIntersectionPoints(vector<Line> lines,OGRPolygon *boundingbox){
//	OGRSFDriver *poDriverPointOut;
//	OGRDataSource *poDSPointOut;
//	OGRLayer *poLayerPointOut;
//	OGRFeature *poFeaturePointOut;
//	OGRFieldDefn oPointField("ID", OFTInteger);
//	oPointField.SetWidth(5);
//
//	poDriverPointOut = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
//	poDSPointOut = poDriverPointOut->CreateDataSource("out_origionalintersectionpoints.shp", NULL);
//	poLayerPointOut = poDSPointOut->CreateLayer("point", NULL, wkbPoint, NULL);
//	poLayerPointOut->CreateField(&oPointField);
//
//	for(int i=0;i<(int)lines.size();i++){
//		OGRPoint *intersections = new OGRPoint[2];
//		OGRPoint point0,point1;
//		lines[i].line->getPoint(0,&point0);
//		lines[i].line->getPoint(1,&point1);
//
//		CalculateIntersections(point0,point1,boundingbox->getExteriorRing(),intersections);
//		lines[i].intersection = intersections;
//
//		//Write intersection points to shp
//		for(int j=0;j<2;j++){
//			try{
//				poFeaturePointOut = OGRFeature::CreateFeature(poLayerPointOut->GetLayerDefn());
//				poFeaturePointOut->SetField("ID", 1);
//				poFeaturePointOut->SetGeometry(&intersections[j]);
//				poLayerPointOut->CreateFeature(poFeaturePointOut);
//			}catch(char *err){
//				printf("Exception raised: %s\n",err);
//			}
//		}
//	}
//	OGRFeature::DestroyFeature(poFeaturePointOut);
//	OGRDataSource::DestroyDataSource(poDSPointOut);
//}
//
//void WriteOGRSimplifiedPolygons(OGRLayer *poLayer){
////Simplification with ogr simplify
//	OGRSFDriver *poDriverSPOut;
//	OGRDataSource *poDSSPOut;
//	OGRLayer *poLayerSPOut;
//	OGRFeature *poFeatureSPOut;
//	OGRFieldDefn oSPField("ID", OFTInteger);
//	oSPField.SetWidth(5);
//
//	OGRFeature *poFeature;
//	OGRGeometry *poGeometry;
//	OGRPolygon *poPolygon;
//
//	poDriverSPOut = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
//	poDSSPOut = poDriverSPOut->CreateDataSource("out_simplepolygons.shp", NULL);
//	poLayerSPOut = poDSSPOut->CreateLayer( "polygon", NULL, wkbPolygon, NULL );
//	poLayerSPOut->CreateField( &oSPField );
//	
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
//		try{
//			poFeatureSPOut = OGRFeature::CreateFeature(poLayerSPOut->GetLayerDefn());
//			poFeatureSPOut->SetField("ID", 1);
//			poFeatureSPOut->SetGeometry(simple);
//			poLayerSPOut->CreateFeature(poFeatureSPOut);
//		}catch(char *err){
//			printf("Exception raised: %s\n",err);
//		}
//	}
//	OGRFeature::DestroyFeature(poFeatureSPOut);
//	OGRDataSource::DestroyDataSource(poDSSPOut);
//
//}