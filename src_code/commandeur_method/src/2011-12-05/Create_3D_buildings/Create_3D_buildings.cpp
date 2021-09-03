#include <vld.h>

//Memory usage
#include <Windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib") //added

#include "ogr_api.h"
#include "SHPReader.h"
#include "SHPWriter.h"
#include "OCIWriter.h"
#include "BuildingBlock.h"
#include <set>

using namespace std;

void writeConvexHull(SHPWriter* convexHullWriter, map<int,Segment> segments, int iBuilding){		
	for(map<int,Segment>::iterator segment = segments.begin();segment!=segments.end();++segment){
		OGRGeometry* geom = segment->second.points->ConvexHull();
			
		if(((OGRPolygon*) geom)->get_Area() > 2){
			OGRFeature *feature = OGRFeature::CreateFeature(convexHullWriter->Layer->GetLayerDefn());
			feature->SetGeometry(segment->second.convexHull);
			feature->SetField("Polygon_id", iBuilding);
			feature->SetField("Segment_id", segment->second.id);
			feature->SetField("Normal_X", segment->second.normal[0]);
			feature->SetField("Normal_Y", segment->second.normal[1]);
			feature->SetField("Normal_Z", segment->second.normal[2]);
			feature->SetField("Nr_points", segment->second.points->getNumGeometries());
			feature->SetField("Area", segment->second.convexHull->get_Area());
			convexHullWriter->writeFeature(feature);
		}	
	}
}

void writeRoofLOD1(SHPWriter* writer, BuildingCell cell, int buildingNr){
	OGRFeature *feature = OGRFeature::CreateFeature(writer->Layer->GetLayerDefn());
	feature->SetGeometry(cell.roof);
	if(cell.roofType==Flat){
		feature->SetField("RoofType", "Flat");
		feature->SetField("RoofHeight", cell.roofHeight);
	}else if(cell.roofType==Shed){
		feature->SetField("RoofType", "Shed");
		feature->SetField("Eaves", cell.eavesHeight);
		feature->SetField("Ridge", cell.ridgeHeight);
		feature->SetField("DirectionX", cell.roofDirection[0]);
		feature->SetField("DirectionY", cell.roofDirection[1]);
	}else if(cell.roofType==NoSegments){
		feature->SetField("RoofType", "NoSegments");
	}else if(cell.roofType==Gabled){
		//write gabled in 2 polygons
		feature->SetGeometry(cell.roof->getGeometryRef(0));
		feature->SetField("RoofType", "Gabled");
		feature->SetField("Eaves", cell.eavesHeight);
		feature->SetField("Ridge", cell.ridgeHeight);
		feature->SetField("DirectionX", cell.roofDirection[0]);
		feature->SetField("DirectionY", cell.roofDirection[1]);
	}else{
		feature->SetField("RoofType", "Unknown");
		feature->SetField("RoofHeight", cell.roofHeight);
	}
	feature->SetField("SurfType", "Roof");
	feature->SetField("RoofAngle", cell.roofAngle);
	feature->SetField("Polygon_id", buildingNr);
	writer->writeFeature(feature);
}

void writeRoof(SHPWriter* writer, BuildingCell cell, int buildingNr){
	OGRFeature *feature = OGRFeature::CreateFeature(writer->Layer->GetLayerDefn());
	feature->SetGeometry(cell.roof);
	if(cell.roofType==Flat){
		feature->SetField("RoofType", "Flat");
		feature->SetField("RoofHeight", cell.roofHeight);
	}else if(cell.roofType==Shed){
		feature->SetField("RoofType", "Shed");
		feature->SetField("Eaves", cell.eavesHeight);
		feature->SetField("Ridge", cell.ridgeHeight);
		feature->SetField("DirectionX", cell.roofDirection[0]);
		feature->SetField("DirectionY", cell.roofDirection[1]);
	}else if(cell.roofType==NoSegments){
		feature->SetField("RoofType", "NoSegments");
	}else if(cell.roofType==Gabled){
		//write gabled in 2 polygons
		feature->SetGeometry(cell.roof->getGeometryRef(0));
		feature->SetField("RoofType", "Gabled");
		feature->SetField("Eaves", cell.eavesHeight);
		feature->SetField("Ridge", cell.ridgeHeight);
		feature->SetField("SurfType", "Roof");
		feature->SetField("RoofAngle", cell.roofAngle);
		feature->SetField("Polygon_id", buildingNr);
		feature->SetField("DirectionX", cell.roofDirection[0]);
		feature->SetField("DirectionY", cell.roofDirection[1]);
		writer->writeFeature(feature);
		feature = OGRFeature::CreateFeature(writer->Layer->GetLayerDefn());
		feature->SetGeometry(cell.roof->getGeometryRef(1));
		//write second roof polygon
		feature->SetField("RoofType", "Gabled");
		feature->SetField("Eaves", cell.eavesHeight);
		feature->SetField("Ridge", cell.ridgeHeight);
		feature->SetField("DirectionX", -cell.roofDirection[0]); //Oposite direction
		feature->SetField("DirectionY", -cell.roofDirection[1]); //Oposite direction
	}else{
		feature->SetField("RoofType", "Unknown");
		feature->SetField("RoofHeight", cell.roofHeight);
	}
	feature->SetField("SurfType", "Roof");
	feature->SetField("RoofAngle", cell.roofAngle);
	feature->SetField("Polygon_id", buildingNr);
	writer->writeFeature(feature);
}

void writeWalls(SHPWriter* writer, OGRMultiPolygon* walls, int buildingNr){
	for(int wall=0;wall<walls->getNumGeometries();wall++){
		OGRFeature *feature = OGRFeature::CreateFeature(writer->Layer->GetLayerDefn());
		feature->SetGeometry(walls->getGeometryRef(wall));
		feature->SetField("SurfType", "Wall");
		feature->SetField("Polygon_id", buildingNr);
		writer->writeFeature(feature);
	}
}

void writeFloor(SHPWriter* writer, OGRPolygon* floor, int buildingNr){
	OGRFeature *feature = OGRFeature::CreateFeature(writer->Layer->GetLayerDefn());
	feature->SetGeometry(floor);
	feature->SetField("SurfType", "Floor");
	feature->SetField("Polygon_id", buildingNr);
	writer->writeFeature(feature);

}

void writeBuildingLOD1(SHPWriter* writer, Building building){
	for(vector<BuildingCell>::iterator buildingPart = building.buildingCells.begin();buildingPart!=building.buildingCells.end();++buildingPart){
		if(buildingPart->roofType==NoSegments){
			writeRoof(writer,*buildingPart,building.id);
		}else{
			writeRoofLOD1(writer,*buildingPart,building.id);
			writeWalls(writer,buildingPart->walls,building.id);
			writeFloor(writer,buildingPart->floor,building.id);
		}
	}
}

void writeBuildingLOD2(SHPWriter* writer, Building building){
	for(vector<BuildingCell>::iterator buildingPart = building.buildingCells.begin();buildingPart!=building.buildingCells.end();++buildingPart){
		if(buildingPart->roofType==NoSegments){
			writeRoof(writer,*buildingPart,building.id);
		}else{
			writeRoof(writer,*buildingPart,building.id);
			writeWalls(writer,buildingPart->walls,building.id);
			writeFloor(writer,buildingPart->floor,building.id);
		}
	}
}

void writeBuildingBlockLOD1(SHPWriter* roofWriter, map<long,Building> block){
	for(map<long,Building>::iterator building = block.begin();building!=block.end();building++){
		writeBuildingLOD1(roofWriter,building->second);
	}
	roofWriter->Layer->SyncToDisk();
}

void writeBuildingBlockLOD2(SHPWriter* roofWriter, map<long,Building> block){
	for(map<long,Building>::iterator building = block.begin();building!=block.end();building++){
		writeBuildingLOD2(roofWriter,building->second);
	}
	roofWriter->Layer->SyncToDisk();
}

bool addBuildingToBlock(map<long,Building> &block, SHPReader* reader, long iBuilding, vector<bool> &buildingsProcessed){
	Building building;
	OGRMultiPolygon* buildingGeom = reader->readBuilding(iBuilding);
	vector<long> adjacentBuildings = reader->getAdjacentBuildings(buildingGeom);

	building.adjacentBuildings.insert(building.adjacentBuildings.end(),adjacentBuildings.begin(),adjacentBuildings.end());

	building.id = iBuilding;
	building.buildingGeom = buildingGeom;
	block[iBuilding] = building;
	return 1;
}

map<long,Building> readBuildingBlock(SHPReader* reader, long iBuilding, vector<bool> &buildingsProcessed){
	map<long,Building> block;
	set<long> processingList;
	processingList.insert(iBuilding);

	while(!processingList.empty()){
		long nextBuilding = *(processingList.begin());
		processingList.erase(nextBuilding);
		if(block.find(nextBuilding) == block.end()){
			if(addBuildingToBlock(block,reader,nextBuilding, buildingsProcessed)){
				processingList.insert(block[nextBuilding].adjacentBuildings.begin(),block[nextBuilding].adjacentBuildings.end());
			}
		}
	}
	return block;
}

void helpMessage(){
	printf("\nOptions:		Description (default)\n");
	printf(" -ip arg		Input polygon shape file\n");
	printf(" -is arg		Input segmented point shape file\n");
	printf(" -o arg			Output 3d shapefile\n");
	printf(" -average arg	Averaging roof height\n");
	printf(" -LOD1			Create an LOD1 model\n");
	printf(" -overwrite		Overwrite existing files\n");
	exit(1);
}

int main(int argc, char* argv[]){
	string input_decomp_filename;
	string input_segm_filename;
	string output_filename;
	bool overwriteFiles = false;
	bool writeOCI = false;
	bool LOD1 = false;
	double average_height = 0.25;

	for(int a = 1; a < argc; a++){
		if(strcmp(argv[a],"-ip") == 0){
			a++;
			input_decomp_filename = string(argv[a]);
		}
		else if(strcmp(argv[a],"-is") == 0){
			a++;
			input_segm_filename = string(argv[a]);
		}
		else if(strcmp(argv[a],"-o") == 0){
			a++;
			output_filename = string(argv[a]);
		}
		else if(strcmp(argv[a],"-overwrite") == 0){
			overwriteFiles = true;
		}
		else if(strcmp(argv[a],"-average") == 0){
			a++;
			average_height = atof(argv[a]);
		}else if(strcmp(argv[a],"-OCI") == 0 || strcmp(argv[a],"-oci") == 0){
			writeOCI = true;
		}else if(strcmp(argv[a],"-LOD1") == 0 || strcmp(argv[a],"-lod1") == 0){
			LOD1 = true;
		}
		else helpMessage();
	}
	if(input_decomp_filename.empty() || input_segm_filename.empty() || output_filename.empty()){
		helpMessage();
	}

	if(overwriteFiles){
		string syntax = "del ";
		syntax.append(output_filename.substr(0,output_filename.find_last_of("."))).append(".*");
		system(syntax.c_str());
		//system("del convexHull.*");
		system("del OCI_SDO_SOLIDS.txt");
	}

	RegisterOGRShape();

	//read shp with decomposed polygons
	SHPReader* decomposedReader = new SHPReader(input_decomp_filename);

	//read shp with points
	SHPReader* segmentReader = new SHPReader(input_segm_filename);

	//SHPWriter* convexHullWriter = new SHPWriter("convexHull.shp",wkbPolygon);
	//convexHullWriter->addField("Polygon_id",OFTInteger);
	//convexHullWriter->addField("Segment_id",OFTInteger);
	//convexHullWriter->addField("Normal_X",OFTInteger);
	//convexHullWriter->addField("Normal_Y",OFTInteger);
	//convexHullWriter->addField("Normal_Z",OFTInteger);
	//convexHullWriter->addField("Nr_points",OFTInteger);
	//convexHullWriter->addField("Area",OFTReal);

	SHPWriter* roofWriter = new SHPWriter(output_filename,wkbPolygon25D);
	roofWriter->addField("Polygon_id",OFTInteger);
	roofWriter->addField("RoofType",OFTString);
	roofWriter->addField("SurfType",OFTString);
	roofWriter->addField("RoofAngle",OFTReal);
	roofWriter->addField("RoofHeight",OFTReal);
	roofWriter->addField("Eaves",OFTReal);
	roofWriter->addField("Ridge",OFTReal);
	roofWriter->addField("DirectionX",OFTInteger);
	roofWriter->addField("DirectionY",OFTInteger);

	PROCESS_MEMORY_COUNTERS pmc;
	clock_t overall_start = clock();

	long totalBlocks = 0;
	long start = 0;
	long end   = decomposedReader->Layer->GetFeatureCount();
	
	long nextPolygonNr = start;
	vector<bool> buildingsProcessed;
	for(int poly=0;poly<end;poly++){
		buildingsProcessed.insert(buildingsProcessed.begin()+poly,false);
	}
	long processedBuildings=0;

	while(nextPolygonNr<end){
		clock_t timer_start = clock();
		try{
			printf("\n============ Building nr: %d ============\n",nextPolygonNr);
			map<long,Building> block = readBuildingBlock(decomposedReader, nextPolygonNr,buildingsProcessed);
			processedBuildings += block.size();
			for(map<long,Building>::iterator building = block.begin();building!=block.end();building++){
				buildingsProcessed[building->second.id]=true;
		
				building->second.segments = findSegmentsInBuilding(segmentReader->Layer,building->second.buildingGeom);
		
				//writeConvexHull(convexHullWriter,building->second.segments,building->second.id);

				//Calculate Building RoofTypes
				processBuilding(building->second,average_height);
			}

			if(block.size() > 1){
				adjustBuildingBlock(block,average_height);
			}

			//For each building; add heights to cells and create walls
			if(LOD1){
				createBuildingBlockGeometryLOD1(block);
				writeBuildingBlockLOD1(roofWriter,block);
			}else{
				createBuildingBlockGeometryLOD2(block);
				writeBuildingBlockLOD2(roofWriter,block);
			}

			if(writeOCI){
				writeBlockToOCI(block);
			}

			while(nextPolygonNr!=end && buildingsProcessed.at(nextPolygonNr)==true){
				nextPolygonNr++;
			}
			totalBlocks++;
		}catch(char* err){
			printf("Building processing failure: %s",err);
		}
		GetProcessMemoryInfo(GetCurrentProcess(), &pmc,sizeof(pmc));
		printf("Total time needed to process building block is %lf seconds\n",(double)(clock()-timer_start)/CLOCKS_PER_SEC);
		printf("Memory footprint: %5lu MB \n", pmc.WorkingSetSize/1048576);
	}
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc,sizeof(pmc));
	printf("Total number of buildings processed: %d\n",processedBuildings);
	printf("Total time needed to create buildings is %lf seconds\n",(double)(clock()-overall_start)/CLOCKS_PER_SEC);
	printf("Average processing time per building: %g seconds\n", ((double)(clock()-overall_start)/CLOCKS_PER_SEC)/processedBuildings);
	printf("Average buildings per block: %ld\n", processedBuildings/totalBlocks);
	printf("Average processing time per block: %g seconds\n", ((double)(clock()-overall_start)/CLOCKS_PER_SEC)/totalBlocks);
	printf("Maximum memory footprint: %5lu MB \n", pmc.PeakWorkingSetSize/1048576);

	decomposedReader->close();
	delete decomposedReader;
	segmentReader->close();
	delete segmentReader;
	//convexHullWriter->close();
	//delete convexHullWriter;
	roofWriter->close();
	delete roofWriter;
}