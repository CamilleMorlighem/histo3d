#include <iostream>
#include "ogrsf_frmts.h"

class SegmentReader{
	public:
		std::string filename;
		FILE* input_file;
		bool last_seg;
		char* line;
		double point[3];
		int normal_vector[3],last_normal_vector[3],segment_id,last_segment_id;
		long pointcount;
		OGRMultiPoint segmentPoints;

		virtual ~SegmentReader();
		SegmentReader(std::string);

		bool SegmentReader::open(std::string);
		void SegmentReader::close();
		bool SegmentReader::readLine();
		bool SegmentReader::createSegment();
};