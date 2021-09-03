#include "ASCIIReader.h"

using namespace std;

ASCIIReader::ASCIIReader(string filename){
	point[0] = point[1] = point[2] = 0;
	input_file = 0;
	line = 0;
	pointcount = 0;

	open(filename);
}

bool ASCIIReader::open(string filename){
	if(filename.empty()){
		cout << "Error: No input file specified!\n";
		return false;
	}

	fopen_s(&input_file, filename.c_str() , "r");
	
	if(input_file == NULL){
		printf("Unable to open file %s\n",filename.c_str());
		exit(1);
	}
	if (line == 0){
		line = (char*)malloc(sizeof(char)*256);
	}
	return true;
}

void ASCIIReader::close(){
	fclose(input_file);
}

bool ASCIIReader::readLine(){
	if(fgets(line, 256, input_file) != NULL){
		if(sscanf_s(line,"%lf,%lf,%lf",&point[0],&point[1],&point[2])==3){
			pointcount++;
			return true;
		}else{
			printf("Could not process line...\n");
		}
	}
	//printf("Could not read line...\n");
	return false;
}