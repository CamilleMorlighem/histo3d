#include "SegmentReader.h"

using namespace std;

SegmentReader::SegmentReader(string filename){
	last_seg = false;
	point[0] = point[1] = point[2] = 0;
	normal_vector[0] = normal_vector[1] = normal_vector[2] = 0;
	segment_id = last_segment_id = 0;
	pointcount = 0;
	input_file = 0;
	line = 0;
	//segmentPoints = new OGRMultiPoint;

	open(filename);
	readLine();
}

SegmentReader::~SegmentReader(){
	//OGRGeometryFactory::destroyGeometry(segmentPoints);
}

bool SegmentReader::open(string filename){
	if(filename.empty()){
		cout << "Error: No input file specified!\n";
		return false;
	}

	fopen_s(&input_file, filename.c_str() , "r");
	
	if(input_file == NULL){
		cout << "Unable to open file...\n";
		exit(1);
	}
	if (line == 0){
		line = (char*)malloc(sizeof(char)*256);
	}
	return true;
}

void SegmentReader::close(){
	fclose(input_file);
	free(line);
}

bool SegmentReader::readLine(){
	if(fgets(line, 256, input_file) != NULL){
		if(sscanf_s(line,"%lf,%lf,%lf,%d,%d,%d,%d",&point[0],&point[1],&point[2],
				&normal_vector[0],&normal_vector[1],&normal_vector[2],&segment_id)==7){
			pointcount++;
			return true;
		}else{
			printf("Could not process line...\n");
		}
	}
	//printf("Could not read line...\n");
	return false;
}

bool SegmentReader::createSegment(){
	this->segmentPoints.empty();
	OGRPoint newpoint;
	int segment=segment_id;
	
	newpoint.setX(this->point[0]);
	newpoint.setY(this->point[1]);
	newpoint.setZ(this->point[2]);
	this->segmentPoints.addGeometry(&newpoint);
	last_segment_id = segment_id;
	last_normal_vector[0] = normal_vector[0];
	last_normal_vector[1] = normal_vector[1];
	last_normal_vector[2] = normal_vector[2];

	while(readLine()){
		if(segment==segment_id){
			newpoint.setX(this->point[0]);
			newpoint.setY(this->point[1]);
			newpoint.setZ(this->point[2]);
			this->segmentPoints.addGeometry(&newpoint);
		}else{
			return true; //Segment is done
		}
	}
	if(!last_seg){
		last_seg = true;
		return true; //EOF
	}
	return false;
/*
	if(readLine()){
		OGRPoint* newpoint = new OGRPoint;
		
		int segment=segment_id;
		while(this->segment_id==segment){
			newpoint->setX(this->point[0]);
			newpoint->setY(this->point[1]);
			newpoint->setZ(this->point[2]);
			this->segmentPoints->addGeometry(newpoint);
			readLine();
		}
		return true;
	}
	return false;*/
}