#include "SHPReader.h"

using namespace std;

SHPReader::SHPReader(string filename){
	this->filename = filename;
	this->DataSource = OGRSFDriverRegistrar::Open(filename.c_str(), FALSE);
    if(this->DataSource == NULL){
		printf("Opening file %s failed...\n", filename.c_str());
        exit(1);
    }
	Layer = this->DataSource->GetLayer(0);
}

void SHPReader::close(){
	OGRDataSource::DestroyDataSource(this->DataSource);
}

OGRMultiPolygon* correctInvalidRing(OGRPolygon* inputPoly){
	OGRMultiPolygon* correctGeom = (OGRMultiPolygon*)OGRGeometryFactory::createGeometry(wkbMultiPolygon);
	OGRPolygon partGeom;

	partGeom.addRing(inputPoly->getExteriorRing());
	correctGeom->addGeometry(&partGeom);
	for(int ringNr=0;ringNr<inputPoly->getNumInteriorRings();ringNr++){
		partGeom.empty();
		partGeom.addRing(inputPoly->getInteriorRing(ringNr));
		correctGeom->addGeometry(&partGeom);
	}
	return correctGeom;
}

OGRMultiPolygon* correctInvalidMultiPolygon(OGRGeometry* inputGeom){
	if(wkbFlatten(inputGeom->getGeometryType()) == wkbPolygon){
		if(((OGRPolygon*)inputGeom)->getNumInteriorRings() > 0){
			return correctInvalidRing((OGRPolygon*)inputGeom);
		}
	}else if(wkbFlatten(inputGeom->getGeometryType()) == wkbMultiPolygon){
		for(int poly=0;poly<((OGRMultiPolygon*)inputGeom)->getNumGeometries();poly++){
			if(((OGRPolygon*)((OGRMultiPolygon*)inputGeom)->getGeometryRef(poly))->getNumInteriorRings() > 0){
				OGRMultiPolygon* tempGeom = correctInvalidRing((OGRPolygon*)((OGRMultiPolygon*)inputGeom)->getGeometryRef(poly));
				((OGRMultiPolygon*)inputGeom)->removeGeometry(poly);
				for(int tempPoly=0;tempPoly<tempGeom->getNumGeometries();tempPoly++){
					((OGRMultiPolygon*)inputGeom)->addGeometry(tempGeom->getGeometryRef(tempPoly));
				}
			}
		}
		return (OGRMultiPolygon*)inputGeom;
	}
	return NULL;
}

int SHPReader::readPolygon(){
	OGRFeature* feature = OGRFeature::CreateFeature(Layer->GetLayerDefn());
	if((feature = Layer->GetNextFeature()) != NULL){
		this->Geometry = feature->GetGeometryRef()->clone();
		if(!this->Geometry->IsValid()){
			this->Geometry = correctInvalidMultiPolygon(this->Geometry);
		}else{
			this->Geometry = OGRGeometryFactory::forceToMultiPolygon(this->Geometry);
		}
		OGRFeature::DestroyFeature(feature);
		return 1;
	}
	OGRFeature::DestroyFeature(feature);
	return 0;
}