#include <vld.h>

//Memory usage
#include <Windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib") //added

#include "ASCIIReader.h"
#include "SegmentReader.h"
#include "SHPReader.h"
#include "SHPWriter.h"
#include <ctime>
#include <vector>

using namespace std;

const double PI = 4*atan(1.0);

SHPWriter* readFromASCII(const char* input_filename){
	ASCIIReader* reader = new ASCIIReader(input_filename);
	SHPWriter* writer = new SHPWriter("out_points.shp",wkbPoint25D);

	OGRPoint point;

	printf("Reading ascii from %s\n",input_filename);
	clock_t timer_start = clock();

	while(reader->readLine()){
		point.setX(reader->point[0]);
		point.setY(reader->point[1]);
		point.setZ(reader->point[2]);

		OGRFeature *feature = OGRFeature::CreateFeature(writer->Layer->GetLayerDefn());
		feature->SetGeometry(&point);
		writer->writeFeature(feature);
			
		if(reader->pointcount%100000==0){
			writer->Layer->SyncToDisk();
			printf("Writting points upto %ld\r",reader->pointcount);
		}
	}
	reader->close();
	printf("Total time needed to process %ld points is %lf seconds\n",reader->pointcount,(double)(clock()-timer_start)/CLOCKS_PER_SEC);

	return writer;
}

OGRMultiPoint* clipBuilding(OGRLayer* inputLayer, OGRGeometry* buildingGeometry, OGRMultiPoint* points){
	//OGRMultiPoint* points = new OGRMultiPoint;
	points->setCoordinateDimension(3);
	
	printf("Clipping building\n");
	if(buildingGeometry != NULL){
		inputLayer->SetSpatialFilter(buildingGeometry->UnionCascaded());

		OGRFeature* feat = OGRFeature::CreateFeature(inputLayer->GetLayerDefn());
		while( (feat = inputLayer->GetNextFeature()) != NULL){
			points->addGeometry(feat->GetGeometryRef());
		}
		OGRFeature::DestroyFeature(feat);
	}
	return points;
}

OGRMultiPoint* clipBuilding(OGRLayer* inputLayer, OGRGeometry* buildingGeometry){
	OGRMultiPoint* points = new OGRMultiPoint;
	points->setCoordinateDimension(3);
	
	printf("Clipping building\n");
	if(buildingGeometry != NULL){
		inputLayer->SetSpatialFilter(buildingGeometry->UnionCascaded());

		OGRFeature* feat = OGRFeature::CreateFeature(inputLayer->GetLayerDefn());
		while( (feat = inputLayer->GetNextFeature()) != NULL){
			points->addGeometry(feat->GetGeometryRef());
		}
		OGRFeature::DestroyFeature(feat);
	}
	return points;
}

double angleWithHorizontalVector(vector<int> vector1){
	//vector<int> vector2;
	//vector2.push_back(0);
	//vector2.push_back(0);
	//vector2.push_back(1);

	//double lenA = sqrt((double)vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);
	//double lenB = sqrt((double)vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);

	//double dotAB = vector1[0]*vector2[0] + vector1[1]*vector2[1] + vector1[2]*vector2[2];

	//return abs(acos(dotAB / (lenA*lenB)) * (180/PI));

	double lenA = sqrt((double)vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);

	return abs(acos(vector1[2] / (lenA)) * (180/PI));
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
	printf(" -i arg			Input filename (ASCII/SHP)\n");
	printf(" -o arg			Output filename\n");
	printf(" -shp			Input type is shapefile\n");
	printf(" -clipsrc arg		Clipping shapefile\n");
	printf(" -spat			Create spatial index\n");
}

int main(int argc, char* argv[]){
	//-i out_points_test.shp -shp -o small_clipped.shp -clipsrc out_decomposedpoly_BAG_test.shp
	string input_filename;
	string output_filename;
	string clip_filename;
	bool ascii = true;
	bool use_spatial_index = false;
	bool write_ascii = false;
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
		else if(strcmp(argv[a],"-shp") == 0){
			ascii = false;
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
		else if(strcmp(argv[a],"-writeascii") == 0){
			write_ascii = true;
		}
	}

	if(input_filename.empty() || output_filename.empty()){
		usage();
		exit(1);
	}

	RegisterOGRShape();

	SHPWriter* writer;
	if(ascii){
		writer = readFromASCII(input_filename.c_str());
	}else{
		writer = new SHPWriter();
		writer->openExisting(input_filename);
	}

	if(use_spatial_index && !spatialFileExists(input_filename)){
		printf("Creating spatial index with depth %d\n",spatial_depth);
		timer_start = clock();

		string pszSQLCommand = "CREATE SPATIAL INDEX ON \"";
		pszSQLCommand.append(writer->Layer->GetName()).append("\" %d",spatial_depth);
		writer->ExecuteSQL(pszSQLCommand.c_str(),NULL,"OGRSQL");
		//writer->CreateSpatialIndex(spatial_depth);
		printf("Total time needed to create spatial index is %lf seconds\n",(double)(clock()-timer_start)/CLOCKS_PER_SEC);
	}

	timer_start = clock();
	printf("Clipping points\n");

	SHPReader* buildingReader = new SHPReader(clip_filename);

	//Set extent of the buildings to the extent of the points
	OGREnvelope envelope;
	writer->Layer->GetExtent(&envelope); //writer is the shapefile with points
	OGRLinearRing ring;
	ring.addPoint(envelope.MinX,envelope.MinY);
	ring.addPoint(envelope.MinX,envelope.MaxY);
	ring.addPoint(envelope.MaxX,envelope.MaxY);
	ring.addPoint(envelope.MaxX,envelope.MinY);
	ring.closeRings();
	buildingReader->Layer->SetSpatialFilter(&ring);

	SHPWriter* segmentWriter = new SHPWriter(output_filename.c_str(),wkbMultiPoint25D);
	segmentWriter->addField("Polygon_id",OFTInteger);
	segmentWriter->addField("Segment_id",OFTInteger);
	segmentWriter->addField("Normal_X",OFTInteger);
	segmentWriter->addField("Normal_Y",OFTInteger);
	segmentWriter->addField("Normal_Z",OFTInteger);
	segmentWriter->addField("Nr_points",OFTInteger);

	PROCESS_MEMORY_COUNTERS pmc;
	int totalBuildings = buildingReader->Layer->GetFeatureCount();
	int building_nr = 0;
	int totalSegments = 0;

	while(buildingReader->readNextGeometry()){
		printf("----------------------------------\nProcessing polygon %ld of %d\n----------------------------------\n", 
					building_nr,totalBuildings);
		OGRMultiPoint points;
		clipBuilding(writer->Layer,buildingReader->Geometry,&points);
		if(points.getNumGeometries() > 0){
			FILE* outasciiFile;
			fopen_s(&outasciiFile, "out_points.xyz" , "w");
			for(int i=0;i<points.getNumGeometries();i++){
				fprintf(outasciiFile,"%.3lf,%.3lf,%.3lf\n", ((OGRPoint*)points.getGeometryRef(i))->getX(),
															((OGRPoint*)points.getGeometryRef(i))->getY(),
															((OGRPoint*)points.getGeometryRef(i))->getZ());
			}
			fclose(outasciiFile);
			delete outasciiFile;

			//Run segmentation
			string segmentation_arguments = "segmentation_complete.exe -i out_points.xyz -o out_points_segm.xyz -minsegsize 20 -maxslope 70 -dim 2";
			if( system(segmentation_arguments.c_str())==0 ){
				//Reading segmented points from ASCII
				SegmentReader* segmentReader = new SegmentReader("out_points_segm.xyz");

				//Writing segments to SHP
				while(segmentReader->createSegment()){
					vector<int> normal;
					normal.push_back(segmentReader->last_normal_vector[0]);
					normal.push_back(segmentReader->last_normal_vector[1]);
					normal.push_back(segmentReader->last_normal_vector[2]);
					if(segmentReader->segmentPoints.getNumGeometries() >= 20
						&& angleWithHorizontalVector(normal) < 70){
						OGRFeature *feature = OGRFeature::CreateFeature(segmentWriter->Layer->GetLayerDefn());
						feature->SetGeometry(&segmentReader->segmentPoints);
						
						feature->SetField("Polygon_id", buildingReader->polygonID);
						feature->SetField("Segment_id", segmentReader->last_segment_id);
						feature->SetField("Normal_X", segmentReader->last_normal_vector[0]);
						feature->SetField("Normal_Y", segmentReader->last_normal_vector[1]);
						feature->SetField("Normal_Z", segmentReader->last_normal_vector[2]);
						feature->SetField("Nr_points", segmentReader->segmentPoints.getNumGeometries());

						segmentWriter->writeFeature(feature);

						if(write_ascii){
							output_filename.erase(output_filename.end()-4,output_filename.end());
							output_filename.append(".xyz");
							FILE* outascii;
							fopen_s(&outascii, output_filename.c_str() , "a");
							for(int i=0;i<segmentReader->segmentPoints.getNumGeometries();i++){
								fprintf(outascii,"%.3lf,%.3lf,%.3lf,%d\n",
																((OGRPoint*)segmentReader->segmentPoints.getGeometryRef(i))->getX(),
																((OGRPoint*)segmentReader->segmentPoints.getGeometryRef(i))->getY(),
																((OGRPoint*)segmentReader->segmentPoints.getGeometryRef(i))->getZ(),totalSegments);
							}
							fclose(outascii);
							delete outascii;

							totalSegments++;
						}
					}
				}
				segmentReader->close();
				delete segmentReader;
			}
			//Remove out_points.xyz
			segmentWriter->Layer->SyncToDisk();
			system("del out_points.xyz");
			system("del out_points_segm.xyz");
		}
		//OGRGeometryFactory::destroyGeometry(points);
		building_nr++;
		printf("----------------------------------\n");
		GetProcessMemoryInfo(GetCurrentProcess(), &pmc,sizeof(pmc));
		printf("Memory footprint: %5lu MB \n", pmc.WorkingSetSize/1048576);
		printf("----------------------------------\n\n");
	}
	buildingReader->close();
	delete buildingReader;
	segmentWriter->close();
	delete segmentWriter;
	writer->close();
	delete writer;
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc,sizeof(pmc));
	printf("Total number of buildings processed: %d\n",building_nr);
	printf("Total time needed to clip points is %lf seconds\n",(double)(clock()-timer_start)/CLOCKS_PER_SEC);
	printf("Average processing time per building: %g seconds\n", ((double)(clock()-timer_start)/CLOCKS_PER_SEC)/building_nr-1);
}