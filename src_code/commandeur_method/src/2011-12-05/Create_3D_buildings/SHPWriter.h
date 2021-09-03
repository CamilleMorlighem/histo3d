#include "ogrsf_frmts.h"

class SHPWriter{
	public:
		OGRLayer* Layer;

		SHPWriter();
		SHPWriter(std::string,OGRwkbGeometryType);
		~SHPWriter();

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