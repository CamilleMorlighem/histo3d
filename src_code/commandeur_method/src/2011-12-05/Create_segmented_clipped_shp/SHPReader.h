#include "ogrsf_frmts.h"
#include "ogr_api.h"

class SHPReader{
	public:
		long polygonID;
		OGRLayer* Layer;
		OGRGeometry* Geometry;

		SHPReader();
		SHPReader(std::string);
		virtual ~SHPReader();
	
		void SHPReader::close();
		int SHPReader::readNextGeometry();
		
	private:
		std::string filename;
		OGRSFDriver* Driver;
		OGRDataSource* DataSource;
};