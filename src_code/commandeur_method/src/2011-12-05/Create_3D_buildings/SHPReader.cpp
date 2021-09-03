#include "SHPReader.h"

using namespace std;

SHPReader::SHPReader(string filename){
	this->filename = filename;
	this->DataSource = OGRSFDriverRegistrar::Open(filename.c_str(), FALSE);
    if(DataSource == NULL){
        printf("Opening file %s failed...\n", filename.c_str());
        exit(1);
    }

	Layer = this->DataSource->GetLayer(0);
}

SHPReader::~SHPReader(){};

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

vector<long> SHPReader::getAdjacentBuildings(OGRMultiPolygon* inputGeom){
	vector<long> buildingData;

	string pszSQLCommand = "SELECT * FROM \"";
	pszSQLCommand.append(Layer->GetName()).append("\"");
	OGRLayer* result = DataSource->ExecuteSQL(pszSQLCommand.c_str(),inputGeom->Buffer(0.5,1),"OGRSQL");

	OGRFeature* feature = OGRFeature::CreateFeature(Layer->GetLayerDefn());
	while(feature = result->GetNextFeature()){
		OGRGeometry* geom = feature->GetGeometryRef()->clone();
		if(!geom->IsValid()){
			geom = correctInvalidMultiPolygon(geom);
		}else{
			geom = OGRGeometryFactory::forceToMultiPolygon(geom);
		}

		if(inputGeom->Distance(geom) < 0.01){
			buildingData.push_back(feature->GetFID());
		}
	}
	this->DataSource->ReleaseResultSet(result);
	OGRFeature::DestroyFeature(feature);
	return buildingData;
}

OGRMultiPolygon* SHPReader::readBuilding(long iBuilding){
	OGRFeature* feature = Layer->GetFeature(iBuilding);
	OGRGeometry* geom = feature->GetGeometryRef()->clone();
	if(!geom->IsValid()){
		geom = correctInvalidMultiPolygon(geom);
	}else{
		geom = OGRGeometryFactory::forceToMultiPolygon(geom);
	}

	OGRFeature::DestroyFeature(feature);
	return (OGRMultiPolygon*)geom;
}

map<long,OGRGeometry*> SHPReader::readBuildingBlock(long iBuilding){
	map<long,OGRGeometry*> blockData;
	OGRGeometryCollection* blockGeom = new OGRGeometryCollection;

	string pszSQLCommand = "SELECT * FROM \"";
	pszSQLCommand.append(Layer->GetName()).append("\"");

	OGRFeature* feat = Layer->GetFeature(iBuilding);
	blockGeom->addGeometry(feat->GetGeometryRef());
	blockData[feat->GetFID()] = feat->GetGeometryRef();

	int numGeometries = 0;
	while(numGeometries!=blockGeom->getNumGeometries()){
		numGeometries = blockGeom->getNumGeometries();
		OGRLayer* result = DataSource->ExecuteSQL(pszSQLCommand.c_str(),blockGeom->Buffer(0.5),"OGRSQL");

		while(OGRFeature* feature = result->GetNextFeature()){
			if(blockData.find(feature->GetFID()) == blockData.end() && blockGeom->Distance(feature->GetGeometryRef()) < 0.01 ){
				blockGeom->addGeometry(feature->GetGeometryRef());
				blockData[feature->GetFID()] = feature->GetGeometryRef();
			}
		}
		this->DataSource->ReleaseResultSet(result);
	}
	printf("Total adjacent adjacent features is %d \n",blockGeom->getNumGeometries());
	return blockData;
}