#include "ogrsf_frmts.h"
#include "ogr_api.h"

class SHPReader{
	public:
		std::string filename;
		OGRLayer* Layer;
		OGRGeometry* Geometry;

		SHPReader();
		SHPReader(std::string);
	
		void SHPReader::close();
		int SHPReader::readPolygon();
		
	private:
		OGRSFDriver* Driver;
		OGRDataSource* DataSource;
};