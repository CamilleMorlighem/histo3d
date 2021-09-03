#include <iostream>

class ASCIIReader{
	public:
		std::string filename;
		FILE* input_file;
		char* line;
		double point[3];
		long pointcount;

		ASCIIReader();
		ASCIIReader(std::string);

		bool ASCIIReader::open(std::string);
		void ASCIIReader::close();
		bool ASCIIReader::readLine();
};