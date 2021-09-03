#include <iostream>
#include <vector>
#include "ogrsf_frmts.h"

class SHPReader{
	public:
		std::string filename;
		OGRLayer* Layer;
		OGRGeometry* Geometry;
		OGRPolygon* Polygon;
		OGRPolygon* BoundingBox;
		OGRPolygon* BoundingBoxDouble;
		OGRGeometryCollection* adjacentBuildings;

		SHPReader();
		SHPReader(std::string);
		virtual ~SHPReader();
		void SHPReader::close();
		int SHPReader::readPolygon(long);
		void SHPReader::calculateBoundingBox();
		void SHPReader::calculateBoundingBoxDouble();
		void SHPReader::findAdjacentBuildings();
		void SHPReader::findAdjacentAdjacentBuildings();
		
	private:
		OGRSFDriver* Driver;
		OGRDataSource* DataSource;
};

class SHPWriter{
	public:
		OGRLayer* Layer;

		SHPWriter();
		SHPWriter(std::string,OGRwkbGeometryType);
		virtual ~SHPWriter();

		int open(std::string,OGRwkbGeometryType);
		int openUpdate();
		int close();
		int writeFeature(OGRFeature*);
		int addField(std::string,OGRFieldType);
		void setFilename(std::string);
		void setLayername(std::string);
		void setGeometrytype(OGRwkbGeometryType);
		std::string getFilename();
		std::string getLayername();
		OGRwkbGeometryType getGeometrytype();

	private:
		OGRSFDriver* Driver;
		OGRDataSource* DataSource;
		std::string filename;
		std::string layername;
		OGRwkbGeometryType geometryType;
};

////TestFunctions
//void WriteOGRSimplifiedPolygons(OGRLayer*);
//void WriteOrigionalIntersectionPoints(std::vector<Line>,OGRPolygon*);