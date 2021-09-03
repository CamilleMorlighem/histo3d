#include "ogrsf_frmts.h"
#include <vector>
#include <map>

class SHPReader{
	public:
		std::string filename;
		OGRLayer* Layer;

		SHPReader();
		SHPReader(std::string);
		~SHPReader();
	
		void SHPReader::close();
		std::vector<long> SHPReader::getAdjacentBuildings(OGRMultiPolygon*);
		//std::map<long,OGRGeometry*> 
		OGRMultiPolygon* SHPReader::readBuilding(long);
		std::map<long,OGRGeometry*> SHPReader::readBuildingBlock(long);
		
	private:
		OGRSFDriver* Driver;
		OGRDataSource* DataSource;
};