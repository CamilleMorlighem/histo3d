#include "SHPReader.h"
#include "SHPWriter.h"
#include <ctime>

using namespace std;

OGRMultiPoint* clipBuilding(OGRLayer* inputLayer, OGRGeometry* buildingGeometry){
	OGRMultiPoint* points = new OGRMultiPoint;
	points->setCoordinateDimension(3);
	
	if(buildingGeometry != NULL){
		inputLayer->SetSpatialFilter(buildingGeometry->UnionCascaded());

		OGRFeature* feat = OGRFeature::CreateFeature(inputLayer->GetLayerDefn());
		while( (feat = inputLayer->GetNextFeature()) != NULL){
			points->addGeometry(feat->GetGeometryRef());
		}
	}
	return points;
}

bool spatialFileExists(string filename){
	filename.erase(filename.end()-4,filename.end());
	filename.append(".qix");
	FILE *istream;
	if ((fopen_s(&istream, filename.c_str(),"r")) == NULL){
		fclose(istream);
		return true;
	}
	return false;
}

void usage(){
	printf("\nOptions:		Description\n");
	printf(" -i arg			Input filename (point SHP)\n");
	printf(" -o arg			Output filename\n");
	printf(" -clipsrc arg		Clipping shapefile (Polygon SHP)\n");
	printf(" -spat			Create spatial index\n");
}

int main(int argc, char* argv[]){
	//-i out_points_test.shp -shp -o small_clipped.shp -clipsrc out_decomposedpoly_BAG_test.shp
	string input_filename;
	string output_filename;
	string clip_filename;
	bool use_spatial_index = false;
	int spatial_depth = 12;
	clock_t timer_start;

	for(int a = 1; a < argc; a++){
		if(strcmp(argv[a],"-i") == 0){
			a++;
			input_filename = string(argv[a]);
		}
		else if(strcmp(argv[a],"-o") == 0){
			a++;
			output_filename = string(argv[a]);
		}
		else if(strcmp(argv[a],"-clipsrc") == 0){
			a++;
			clip_filename = string(argv[a]);
		}
		else if(strcmp(argv[a],"-spat") == 0){
			use_spatial_index = true;
			a++;
			if(atoi(argv[a]) > 0 && atoi(argv[a]) <= 12){
				spatial_depth = atoi(argv[a]);
			}else a--;
		}
	}

	if(input_filename.empty() || output_filename.empty() || clip_filename.empty()){
		usage();
		exit(1);
	}

	//OGRRegisterAll();
	RegisterOGRShape();

	SHPWriter* pointReader = new SHPWriter();
	pointReader->openExisting(input_filename);

	if(use_spatial_index && !spatialFileExists(input_filename)){
		printf("Creating spatial index with depth %d\n",spatial_depth);
		timer_start = clock();

		string pszSQLCommand = "CREATE SPATIAL INDEX ON \"";
		pszSQLCommand.append(pointReader->Layer->GetName()).append("\" %d",spatial_depth);
		pointReader->ExecuteSQL(pszSQLCommand.c_str(),NULL,"OGRSQL");
		//writer->CreateSpatialIndex(spatial_depth);
		printf("Total time needed to create spatial index is %lf seconds\n",(double)(clock()-timer_start)/CLOCKS_PER_SEC);
	}

	timer_start = clock();
	printf("Clipping points\n");

	SHPReader* buildingReader = new SHPReader(clip_filename);

	//Set extent of the buildings to the extent of the points
	OGREnvelope envelope;
	pointReader->Layer->GetExtent(&envelope); //writer is the shapefile with points
	OGRLinearRing* ring = new OGRLinearRing;
	ring->addPoint(envelope.MinX,envelope.MinY);
	ring->addPoint(envelope.MinX,envelope.MaxY);
	ring->addPoint(envelope.MaxX,envelope.MaxY);
	ring->addPoint(envelope.MaxX,envelope.MinY);
	ring->closeRings();
	buildingReader->Layer->SetSpatialFilter(ring);

	SHPWriter* clippedWriter = new SHPWriter(output_filename.c_str(),wkbPoint25D);

	int totalBuildings = buildingReader->Layer->GetFeatureCount();
	int building_nr = 0;
	int totalSegments = 0;
	//testing single building to ascii
	while(buildingReader->readPolygon()){
		printf("Processing polygon %5d of %5d\r",building_nr,totalBuildings);
		OGRMultiPoint* points = clipBuilding(pointReader->Layer,buildingReader->Geometry);
		if(points->getNumGeometries() > 0){
			for(int i=0;i<points->getNumGeometries();i++){
				OGRPoint* point = (OGRPoint*)points->getGeometryRef(i);
				OGRFeature *feature = OGRFeature::CreateFeature(clippedWriter->Layer->GetLayerDefn());
				feature->SetGeometry(point);
				clippedWriter->writeFeature(feature);
			}
		}
		if(building_nr%50==0){
			clippedWriter->Layer->SyncToDisk();
		}
		building_nr++;
	}

	clippedWriter->close();
	pointReader->close();
	printf("Total number of buildings processed: %d\n",building_nr);
	printf("Total time needed to clip points is %lf seconds\n",(double)(clock()-timer_start)/CLOCKS_PER_SEC);
	printf("Average processing time per building: %g seconds\n", ((double)(clock()-timer_start)/CLOCKS_PER_SEC)/building_nr);
}