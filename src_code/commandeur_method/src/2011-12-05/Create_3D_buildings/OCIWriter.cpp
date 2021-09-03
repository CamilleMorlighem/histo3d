#include "OCIWriter.h"
#include <vector>

using namespace std;

FILE* open(string filename){
	FILE* file;
	if(filename.empty()){
		printf("Error: No input file specified!\n");
		return false;
	}

	fopen_s(&file, filename.c_str() , "a");
	
	if(file == NULL){
		printf("Unable to open file %s\n",filename.c_str());
		exit(1);
	}
	return file;
}

void close(FILE* file){
	fclose(file);
}

void writeSolidElementInfo(Building building, FILE* file){
	int totalCells = 0;
	int pointCount=1;
	for(vector<BuildingCell>::iterator cell = building.buildingCells.begin(); cell!=building.buildingCells.end();cell++){
		if(cell->roofType!=NoSegments){
			totalCells++;
		}
	}

	if(totalCells > 1){
		fprintf(file,"1, 1008, %d,\n",totalCells);
	}

	for(vector<BuildingCell>::iterator cell = building.buildingCells.begin(); cell!=building.buildingCells.end();cell++){
		if(cell->roofType != NoSegments){
			int nrFaces = 0;
			nrFaces += cell->roof->getNumGeometries();
			nrFaces += cell->walls->getNumGeometries();
			nrFaces += 1; //floor
			if(pointCount != 1) fprintf(file,",\n");
			fprintf(file,"%d, 1007, 1,\n",pointCount);
			fprintf(file,"%d, 1006, %d,\n",pointCount,nrFaces);
			fprintf(file,"%d, 1003, 1",pointCount);
			pointCount+=((OGRPolygon*)cell->roof->getGeometryRef(0))->getExteriorRing()->getNumPoints()*3;
			for(int roofPart=1;roofPart<cell->roof->getNumGeometries();roofPart++){ //roof
				fprintf(file,",\n%d, 1003, 1",pointCount);
				pointCount+=((OGRPolygon*)cell->roof->getGeometryRef(roofPart))->getExteriorRing()->getNumPoints()*3;
			}
			for(int wallPart=0;wallPart<cell->walls->getNumGeometries();wallPart++){ //walls
				fprintf(file,",\n%d, 1003, 1",pointCount);
				pointCount+=((OGRPolygon*)cell->walls->getGeometryRef(wallPart))->getExteriorRing()->getNumPoints()*3;
			}
			fprintf(file,",\n%d, 1003, 1",pointCount);
			pointCount+=cell->floor->getExteriorRing()->getNumPoints()*3; //floor
			
		}
	}
	fprintf(file,")\n,");
}

void writeSolidOrdinate(Building building, FILE* file){
	bool first = true;
	for(vector<BuildingCell>::iterator cell = building.buildingCells.begin(); cell!=building.buildingCells.end();cell++){
		if(cell->roofType != NoSegments){
			if(!first) fprintf(file,",\n");
			for(int roofPart=0;roofPart<cell->roof->getNumGeometries();roofPart++){ //roof
				OGRLinearRing* ring = ((OGRPolygon*)cell->roof->getGeometryRef(roofPart))->getExteriorRing(); 
				for(int pointnr=0;pointnr<ring->getNumPoints();pointnr++){
					OGRPoint point;
					ring->getPoint(pointnr, &point);
					if(pointnr==0)
						fprintf(file,"%lf, %lf, %lf",point.getX(),point.getY(),point.getZ());
					else
						fprintf(file,",   %lf, %lf, %lf",point.getX(),point.getY(),point.getZ());
				}
				fprintf(file,",  --Roofpart");
				fprintf(file,"\n");
				
			}
			for(int wallPart=0;wallPart<cell->walls->getNumGeometries();wallPart++){ //walls
				OGRLinearRing* ring = ((OGRPolygon*)cell->walls->getGeometryRef(wallPart))->getExteriorRing(); 
				for(int pointnr=0;pointnr<ring->getNumPoints();pointnr++){
					OGRPoint point;
					ring->getPoint(pointnr, &point);
					if(pointnr==0)
						fprintf(file,"%lf, %lf, %lf",point.getX(),point.getY(),point.getZ());
					else
						fprintf(file,",   %lf, %lf, %lf",point.getX(),point.getY(),point.getZ());
				}
				fprintf(file,",  --Wallpart");
				fprintf(file,"\n");
				
			}
			OGRLinearRing* ring = cell->floor->getExteriorRing(); //floor
			for(int pointnr=0;pointnr<ring->getNumPoints();pointnr++){
				OGRPoint point;
				ring->getPoint(pointnr, &point);
				if(pointnr==0)
					fprintf(file,"%lf, %lf, %lf",point.getX(),point.getY(),point.getZ());
				else
					fprintf(file,",   %lf, %lf, %lf",point.getX(),point.getY(),point.getZ());
			}
			//fprintf(file,",  --Bottom");
			first=false;
		}
	}
	fprintf(file,")\n");
}

void writeBuildingToOCI(Building building){
	int totalCells = 0;
	for(vector<BuildingCell>::iterator cell = building.buildingCells.begin(); cell!=building.buildingCells.end();cell++){
		if(cell->roofType!=NoSegments){
			totalCells++;
		}
	}
	if(totalCells>0){
		FILE* file = open("OCI_SDO_SOLIDS.txt");
		fprintf(file,"insert into SOLID3d values (%ld, sdo_geometry(3008, NULL, NULL,\n",building.id);
		fprintf(file,"\tsdo_elem_info_array(\n--\n");
		writeSolidElementInfo(building,file);

		fprintf(file,"\nsdo_ordinate_array(\n");
		writeSolidOrdinate(building,file);
		fprintf(file,"));\n\n");

		close(file);
	}
}

void writeBlockToOCI(map<long,Building> block){
	for(map<long,Building>::iterator building=block.begin();building!=block.end();building++){
		writeBuildingToOCI(building->second);
	}
}