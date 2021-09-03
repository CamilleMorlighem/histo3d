#include "Cell_Decomposition.h"

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

void SHPReader::close(){
	OGRDataSource::DestroyDataSource(this->DataSource);
}

SHPReader::~SHPReader(){};

int SHPReader::readPolygon(long iPolygon){
	this->Geometry = Layer->GetFeature(iPolygon)->GetGeometryRef();

	if(this->Geometry != NULL && wkbFlatten(Geometry->getGeometryType()) == wkbPolygon ){
		this->Polygon = (OGRPolygon*) Geometry;
		if(!this->Polygon->IsValid()){
			this->Polygon->closeRings();
			if(this->Polygon->IsValid()){
				printf("Polygon %d is valid after closing rings\n",iPolygon);
			}else{
				throw "Input Polygon is invalid\n";
			}
		}
	}else{
		printf("No Polygon geometry but %s for %d\n",this->Geometry->getGeometryName(),iPolygon);
		return 0;
    }
	return 1;
}

void SHPReader::calculateBoundingBox(){
	OGREnvelope envelope;
	OGRLinearRing ring;
	this->BoundingBox = (OGRPolygon*) OGRGeometryFactory::createGeometry(wkbPolygon);
	this->Geometry->getEnvelope(&envelope);

	ring.addPoint(envelope.MinX,envelope.MinY);
	ring.addPoint(envelope.MinX,envelope.MaxY);
	ring.addPoint(envelope.MaxX,envelope.MaxY);
	ring.addPoint(envelope.MaxX,envelope.MinY);
	ring.addPoint(envelope.MinX,envelope.MinY);

	this->BoundingBox->addRing(&ring);
}

void SHPReader::calculateBoundingBoxDouble(){
	OGREnvelope envelope;
	OGRLinearRing ring;
	this->BoundingBoxDouble = (OGRPolygon*) OGRGeometryFactory::createGeometry(wkbPolygon);
	this->Geometry->getEnvelope(&envelope);

	ring.addPoint(envelope.MinX-((envelope.MaxX-envelope.MinX)/2),envelope.MinY-((envelope.MaxY-envelope.MinY)/2));
	ring.addPoint(envelope.MinX-((envelope.MaxX-envelope.MinX)/2),envelope.MaxY+((envelope.MaxY-envelope.MinY)/2));
	ring.addPoint(envelope.MaxX+((envelope.MaxX-envelope.MinX)/2),envelope.MaxY+((envelope.MaxY-envelope.MinY)/2));
	ring.addPoint(envelope.MaxX+((envelope.MaxX-envelope.MinX)/2),envelope.MinY-((envelope.MaxY-envelope.MinY)/2));
	ring.addPoint(envelope.MinX-((envelope.MaxX-envelope.MinX)/2),envelope.MinY-((envelope.MaxY-envelope.MinY)/2));

	this->BoundingBoxDouble->addRing(&ring);
}

void SHPReader::findAdjacentBuildings(){
	if(this->Geometry != NULL){
		//printf("Org features is %d \n",this->Layer->GetFeatureCount());
		this->adjacentBuildings = new OGRGeometryCollection;
		string pszSQLCommand = "SELECT * FROM \"";
		pszSQLCommand.append(Layer->GetName()).append("\"");
		OGRLayer* result = DataSource->ExecuteSQL(pszSQLCommand.c_str(),this->Geometry,"OGRSQL");
		//DataSource->ExecuteSQL("CREATE SPATIAL INDEX ON \"bld - kopie\"",NULL,"");
		//this->BoundingBoxDouble->dumpReadable(NULL);
		//printf("Found features is %d \n",result->GetFeatureCount());

		while(OGRFeature* feature = result->GetNextFeature()){
			if(this->Geometry->Touches(feature->GetGeometryRef())){
				this->adjacentBuildings->addGeometry(feature->GetGeometryRef());
			}
		}
		//this->adjacentBuildings->dumpReadable(NULL);
		//printf("Touching features is %d \n",this->adjacentBuildings->getNumGeometries());
		this->DataSource->ReleaseResultSet(result);
	}else printf("Adjacent buildings: Geometry is NULL");
}

void SHPReader::findAdjacentAdjacentBuildings(){
	string pszSQLCommand = "SELECT * FROM \"";
	pszSQLCommand.append(Layer->GetName()).append("\"");

	if(this->Geometry != NULL){
		//printf("Org features is %d \n",this->Layer->GetFeatureCount());
		this->adjacentBuildings = (OGRGeometryCollection*) OGRGeometryFactory::createGeometry(wkbGeometryCollection);

		OGRLayer* result = DataSource->ExecuteSQL(pszSQLCommand.c_str(),this->Geometry,"OGRSQL");
		//printf("Found features is %d \n",result->GetFeatureCount());

		while(OGRFeature* feature = result->GetNextFeature()){
			bool validFeature=false;
			if(feature->GetGeometryRef()->IsValid()){
				validFeature=true;
			}else{
				((OGRPolygon*)feature->GetGeometryRef())->closeRings();
				if(feature->GetGeometryRef()->IsValid()) validFeature=true;
			}
			if(validFeature && this->Geometry->Touches(feature->GetGeometryRef()) && !feature->GetGeometryRef()->Equals(this->Geometry)){
				this->adjacentBuildings->addGeometry(feature->GetGeometryRef());
			}
		}
		//this->adjacentBuildings->dumpReadable(NULL);
		//printf("Touching features is %d \n",this->adjacentBuildings->getNumGeometries());
		this->DataSource->ReleaseResultSet(result);
	}else printf("Adjacent buildings: Geometry is NULL");

	int numGeometries = 0;
	while(numGeometries!=this->adjacentBuildings->getNumGeometries()){
		numGeometries = this->adjacentBuildings->getNumGeometries();
		OGRLayer* result = DataSource->ExecuteSQL(pszSQLCommand.c_str(),this->adjacentBuildings,"OGRSQL");
		//printf("Found features is %d \n",result->GetFeatureCount());

		while(OGRFeature* feature = result->GetNextFeature()){
			bool validFeature=false;
			if(feature->GetGeometryRef()->IsValid()){
				validFeature=true;
			}else{
				((OGRPolygon*)feature->GetGeometryRef())->closeRings();
				if(feature->GetGeometryRef()->IsValid()) validFeature=true;
			}
			if(validFeature && !feature->GetGeometryRef()->Equals(this->Geometry) && !this->adjacentBuildings->Contains(feature->GetGeometryRef()) && this->adjacentBuildings->Touches(feature->GetGeometryRef())){
				this->adjacentBuildings->addGeometry(feature->GetGeometryRef());
			}
		}
		//this->adjacentBuildings->dumpReadable(NULL);
		//printf("Touching features is %d \n",this->adjacentBuildings->getNumGeometries());
		this->DataSource->ReleaseResultSet(result);
	}
	printf("Total adjacent adjacent features is %d \n",this->adjacentBuildings->getNumGeometries());
}