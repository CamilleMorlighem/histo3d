#include "SHPWriter.h"

using namespace std;

SHPWriter::SHPWriter(){}

SHPWriter::SHPWriter(string filename,OGRwkbGeometryType type){
	SHPWriter::open(filename,type);
}

SHPWriter::~SHPWriter(){};

int SHPWriter::open(string filename,OGRwkbGeometryType type){
	this->filename = filename;
	this->geometryType = type;
	this->layername = filename;
	this->Driver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
	
	this->DataSource = Driver->CreateDataSource(this->filename.c_str(), NULL);
	if(DataSource == NULL){
        printf("Creation of output file %s failed\n",this->filename.c_str());
		exit(1);
    }

	this->Layer = DataSource->CreateLayer(this->layername.c_str(), NULL, this->geometryType, NULL);
	if(Layer == NULL){
        printf("Layer creation failed for output file %s\n",this->filename.c_str());
        return 0;
    }
	return 1;
}

int SHPWriter::close(){
	OGRDataSource::DestroyDataSource(this->DataSource);

	return 1;
}

int SHPWriter::writeFeature(OGRFeature *feature){
	if(this->Layer->CreateFeature(feature) != OGRERR_NONE){
        printf("Failed to create feature in shapefile\n");
        return 0;
    }
	OGRFeature::DestroyFeature(feature);
	return 1;
}

int SHPWriter::addField(string fieldName,OGRFieldType fieldType){
	OGRFieldDefn Field(fieldName.c_str(), fieldType);
	//FieldOut->SetWidth(5);
	//FieldOut->SetPrecision(5);
	if( Layer->CreateField( &Field ) != OGRERR_NONE ){
        printf("Creating field %s failed\n",fieldName.c_str());
        return 0;
    }

	return 1;
}

void SHPWriter::setFilename(string filename){
	this->filename = filename;
}

void SHPWriter::setLayername(string layername){
	this->layername = layername;
}

void SHPWriter::setGeometrytype(OGRwkbGeometryType geometryType){
	this->geometryType = geometryType;
}

string SHPWriter::getFilename(){
	return this->filename;
}

string SHPWriter::getLayername(){
	return this->layername;
}

OGRwkbGeometryType SHPWriter::getGeometrytype(){
	return this->geometryType;
}