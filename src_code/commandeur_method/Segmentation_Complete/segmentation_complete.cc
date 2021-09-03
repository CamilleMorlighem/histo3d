/*
--------------------------------------------------------------------------------
                               Include files
--------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LaserUnit.h"

/*
--------------------------------------------------------------------------------
                      The main laser2ascii function
--------------------------------------------------------------------------------
*/


/*
    Copyright 2010 University of Twente and Delft University of Technology
 
       This file is part of the Mapping libraries and tools, developed
  for research, education and projects in photogrammetry and laser scanning.

  The Mapping libraries and tools are free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation, either version 3 of the License,
                   or (at your option) any later version.

 The Mapping libraries and tools are distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
                GNU General Public License for more details.

      You should have received a copy of the GNU General Public License
          along with the Mapping libraries and tools.  If not, see
                      <http://www.gnu.org/licenses/>.

----------------------------------------------------------------------------*/


/*
--------------------------------------------------------------------------------
 Segmentation of laser altimetry data.
 Files can be specified explicitly or by a file filter.
 The segmented points are written to file. Optionally, a meta data file is
 generated.

 Initial creation:
 Author : George Vosselman
 Date   : 19-12-2006

*/
/*
--------------------------------------------------------------------------------
                               Include files
--------------------------------------------------------------------------------
*/

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <strings.h>
//#include <malloc.h>
//#include <time.h>
#include "LaserBlock.h"
//#include "TINEdges.h"
//#include "BNF_io.h"
//#include <iostream>
//#include <ANN.h>
//#include "KNNFinder.h"
//#include "LaserPoints.h"
//#include <math.h>
//#include "LaserDataTypes.h"

/*
--------------------------------------------------------------------------------
                         Declaration of C functions
--------------------------------------------------------------------------------
*/

extern "C" void parse_filter(char **, char *);
extern "C" char *get_full_filename(const char *, const char *, int *);

/*
--------------------------------------------------------------------------------
                         The main segmentlaser function
--------------------------------------------------------------------------------
*/

void timer_start(clock_t *time1)
{
  *time1 = clock();
}

void timer_end(clock_t time1, char *string)
{
  clock_t time2;
  time2 = clock();
  printf("Time used for %s: %5.2f minutes\n\n", string,
         (double) (time2 - time1) / (60 * CLOCKS_PER_SEC));
}

LaserBlock readASCII(char *ascii_file, int column_x, int column_y, int column_z, int header_lines, int remove_double_points){
  LaserBlock           block;
  LaserSubUnit         pointset;
  LaserPoint           point, previous_point;
  FILE                 *ascii;
  char                 line[2048], *comma;
  double               value[21];
  int                  i, double_points, total_points;
  
  // START Reading from ascii
    // ====================================================================
    /* Open the input file */
    ascii = Open_Compressed_File(ascii_file, "r");
    if (!ascii) {
    fprintf(stderr, "Error opening input file %s\n", ascii_file);
    exit(0);
    }
    
    /* Skip the header records */
    
    for (i=0; i<header_lines; i++) fgets(line, 2048, ascii);
    if (feof(ascii)) {
    fprintf(stderr, "Error: end of file reached after reading header lines.\n");
    exit(0);
    }
    
    // Set the point type
    pointset.Scanner().SetPointType(NormalPoint);

    /* Process all records */
    total_points = double_points = 0;
    do {
    if (fgets(line, 1024, ascii)) {
    
    // Remove the comma's
    
      while ((comma = strchr(line, ',')) != NULL) *comma = ' ';
    
    // Read the next line
    
      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
             value+ 1, value+ 2, value+ 3, value+ 4, value+ 5,
             value+ 6, value+ 7, value+ 8, value+ 9, value+10,
             value+11, value+12, value+13, value+14, value+15,
             value+16, value+17, value+18, value+19, value+20);
    
    // Copy the data to the laser point
    
      point.X() = value[column_x];
      point.Y() = value[column_y];
      point.Z() = value[column_z];
    
    // Put the laser point into the data set
    
      if (point != previous_point || !remove_double_points) {
        pointset.push_back(point);
        previous_point = point;
      }
      else double_points++;
        
    // Counter
    
      if ((pointset.size() / 10000) * 10000 == pointset.size()) {
        printf(" %d / %d\r", pointset.size(), double_points);
        fflush(stdout);
      }
    }
    
    // Output of subunit if maximum size or EOF is reached
    
    if (feof(ascii)) {      
      total_points += pointset.size();
      block.push_back(pointset);   
    }
    } while (!feof(ascii));
    Close_Compressed_File(ascii);
    
    printf("Total of %d records read.\n\n", total_points + double_points);
    if (double_points)
      printf("%d duplicate records encountered.\n", double_points);
           
  return block;
}

void writeASCII(LaserPoints laser_out, char *outfile){
  LaserPoints::iterator  out_point;
  FILE   *ascii;
  int    store_x = true;
  int    store_y = true;
  int    store_z = true;
  int    store_normx = true;
  int    store_normy = true;
  int    store_normz = true;
  int    store_sn = true;
  
  // Open the output file
  ascii = fopen(outfile, "w");
  if (!ascii) {
    fprintf(stderr, "Error opening output file %s\n", outfile);
    exit(0);
  }

  for (out_point=laser_out.begin(); out_point!=laser_out.end(); out_point++) {
    if (store_x) fprintf(ascii, "%.3f, ", out_point->X());
    if (store_y) fprintf(ascii, "%.3f, ", out_point->Y());
    if (store_z) fprintf(ascii, "%.3f, ", out_point->Z());
    if (store_normx) fprintf(ascii, "%d, ", out_point->Reflectance());
    if (store_normy) fprintf(ascii, "%d, ", out_point->Label());
    if (store_normz) fprintf(ascii, "%d, ", out_point->PulseLength());
    if (store_sn) fprintf(ascii, "%d", out_point->SegmentNumber());
    fprintf(ascii, "\n");
  }
  
  printf("Written %s with %d points.\n", outfile,laser_out.size());
  
  fclose(ascii);
}

void segmentation_complete(char *infile, char *outfile,
                  int column_x, int column_y, int column_z,
                  int remove_double_points, int header_lines,
                  const SegmentationParameters &parameters){
    
    //Segmentation
    LaserBlock           block;
    LaserBlock::iterator unitptr;
    LaserUnit::iterator  subunitptr;
    clock_t              start;
    
    //Normal vectors
    LaserPoints          laser_points, sel_laser_points, laser_out;
    vector <int>         seg_nums;
    vector <int>::iterator seg_id;
    
    int                             dir_x, dir_y, dir_z, index2;
    double                          PI=3.14159; 
    FILE                            *ascii;
    Plane                           plane;
    PointNumberList                 pnl_seg;
    
    // START Reading from ascii
    printf("Begin reading from %s\n",infile);
    block = readASCII(infile, column_x, column_y, column_z, header_lines, remove_double_points);
    // END Reading from ascii
    
    printf("Begin of segmentation\n");

    // Loop over all units and subunits
    for (unitptr=block.begin(); unitptr!=block.end(); unitptr++) {
      for (subunitptr=unitptr->begin(); subunitptr!=unitptr->end();
	       subunitptr++) {

          //Sort points
          subunitptr->SortOnCoordinates();
          
          // START segmentation of the data
          timer_start(&start);
          subunitptr->SurfaceGrowing(parameters);
          timer_end(start, "segmenting");
          // END segmentation of the data
          

          // START Adding normal vectors to the segments
          printf("Begin calculation of normal vectors\n");
          //get segment numbers
          seg_nums = subunitptr->AttributeValues(SegmentNumberTag);
            
          if (seg_nums.size() <= 1){
            cout<<"WARNING: input data set contains only one segment number or is not segmented."<<endl;
          }
          printf("Number of segments found: %d\n", seg_nums.size());
               
          //iterate over segment numbers
          for (seg_id = seg_nums.begin(), index2=0; seg_id!=seg_nums.end(); seg_id++, index2++){
            printf("%7d  %5.1f% \r", index2, 100.0 * index2 / seg_nums.size());
            //get points from segment
            sel_laser_points = subunitptr->SelectTagValue(SegmentNumberTag, *seg_id);
            //get point number list for segment
            pnl_seg = sel_laser_points.TaggedPointNumberList(SegmentNumberTag,*seg_id);
                
            // fit plane
            plane = sel_laser_points.FitPlane(*seg_id, *seg_id, SegmentNumberTag);
                
            // assign normal vector to laser points
            dir_x = int(plane.Normal()[0]*180/PI); // divide by 2 for making all angles less than 180 (needed for 8 bit representation)
            dir_y = int(plane.Normal()[1]*180/PI);
            dir_z = int(plane.Normal()[2]*180/PI);
                
            //Make sure that the Z component is positive
            if (dir_z < 0){
               dir_x = -dir_x;
               dir_y = -dir_y;
               dir_z = -dir_z;
            }
                  
            //Make the vector unitary (and multiply by 100 so it can be stored in an integer -required by parameter type-)
            double magnitude = sqrt((dir_x*dir_x)+(dir_z*dir_z)+(dir_z*dir_z));
            dir_x = int((dir_x/magnitude)*100);
            dir_y = int((dir_y/magnitude)*100);
            dir_z = int((dir_z/magnitude)*100);
                
            sel_laser_points.SetAttribute(ReflectanceTag, dir_x);
            sel_laser_points.SetAttribute(LabelTag, dir_y);
            sel_laser_points.SetAttribute(PulseLengthTag , dir_z);

            //erase points
            laser_out.AddPoints(sel_laser_points);
            sel_laser_points.ErasePoints();
          }
          printf("Points not belonging to segments %d\n",subunitptr->size()-laser_out.size());
          // END Adding normal vectors to the segments
                
                   
          // START Remove small segments
             //TODO?
          // END Remove small segments
          
          
          // START Write the point data to the ASCII file
             writeASCII(laser_out,outfile);
          // END Write the point data to the ASCII file
          
          // Erase the points
          laser_out.ErasePoints();
          subunitptr->ErasePoints();
      } //END for subunitptr
    } //END for unitptr
} //END segmentation_complete
