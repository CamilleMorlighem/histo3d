#ifndef _BUILDINGBLOCK_H_INCLUDED
#define _BUILDINGBLOCK_H_INCLUDED
#include "BuildingBlock.h"
#endif
#include <set>

using namespace std;

void createBuildingBlockGeometryLOD1(map<long,Building> &block){
	for(map<long,Building>::iterator building = block.begin();building!=block.end();building++){
		createBuildingGeometryLOD1(building->second.buildingCells);
	}
}

void createBuildingBlockGeometryLOD2(map<long,Building> &block){
	for(map<long,Building>::iterator building = block.begin();building!=block.end();building++){
		createBuildingGeometryLOD2(building->second.buildingCells);
	}
}

double angleDifferenceVectors(vector<int> vector1, vector<int> vector2){
	if(vector1==vector2){
		return 0;
	}else{
		double lenA = sqrt((double)vector1[0]*vector1[0] + vector1[1]*vector1[1] + vector1[2]*vector1[2]);
		double lenB = sqrt((double)vector2[0]*vector2[0] + vector2[1]*vector2[1] + vector2[2]*vector2[2]);

		double dotAB = (vector1[0]*vector2[0]) + (vector1[1]*vector2[1]) + (vector1[2]*vector2[2]);
		double bla = abs(acos(dotAB / (lenA*lenB)) *(180/PI));
		double temp = abs(acos(dotAB / (lenA*lenB)) *(180/PI));
		return abs(acos(dotAB / (lenA*lenB)) *(180/PI));
	}
}

void adjustBuildingsWithSameRoof(map<long,Building> &block, set<long> sameRoof){
	if(block[*(sameRoof.begin())].buildingCells.at(0).roofType == Flat
			|| block[*(sameRoof.begin())].buildingCells.at(0).roofType == UnknownType){
		double height = 0;
		int roofCount = 0;
		for(set<long>::iterator bld = sameRoof.begin();bld!=sameRoof.end();bld++){
			height += block[*bld].buildingCells.at(0).roofHeight;
			roofCount++;
		}
		double averageHeight = ceil((height/roofCount)*100)/100;
		for(set<long>::iterator bld = sameRoof.begin();bld!=sameRoof.end();bld++){
			block[*bld].buildingCells.at(0).roofHeight = averageHeight;
		}
	}else if(block[*(sameRoof.begin())].buildingCells.at(0).roofType == Shed
			|| block[*(sameRoof.begin())].buildingCells.at(0).roofType == Gabled){
		double minEavesHeight = block[*(sameRoof.begin())].buildingCells.at(0).eavesHeight;
		double minRidgeHeight = block[*(sameRoof.begin())].buildingCells.at(0).ridgeHeight;
		for(set<long>::iterator bld = sameRoof.begin();bld!=sameRoof.end();bld++){
			if(block[*bld].buildingCells.at(0).eavesHeight < minEavesHeight)
				minEavesHeight = block[*bld].buildingCells.at(0).eavesHeight;
			if(block[*bld].buildingCells.at(0).ridgeHeight < minRidgeHeight)
				minRidgeHeight = block[*bld].buildingCells.at(0).ridgeHeight;
		}
		minEavesHeight = ceil(minEavesHeight*100)/100;
		minRidgeHeight = ceil(minRidgeHeight*100)/100;
		for(set<long>::iterator bld = sameRoof.begin();bld!=sameRoof.end();bld++){
			block[*bld].buildingCells.at(0).eavesHeight = minEavesHeight;
			block[*bld].buildingCells.at(0).ridgeHeight = minRidgeHeight;
		}
	}
}

set<long> findBuildingsWithSameRoof(map<long,Building> &block,long buildingNr){
	vector<long> processList;
	set<long> sameRoof;
	//vector<long> differentRoof;
	
	processList.push_back(buildingNr);

	while(!processList.empty()){
		long nextBld = processList.back();
		processList.pop_back();
		
		sameRoof.insert(nextBld);
		//Building has one rooftype
		if(block[nextBld].buildingCells.size() == 1){
			for(vector<long>::iterator adjBuilding = block[nextBld].adjacentBuildings.begin();adjBuilding!=block[nextBld].adjacentBuildings.end();adjBuilding++){
				//Adjacent building has one rooftype, same rooftype and same direction
				if(sameRoof.find(*adjBuilding) == sameRoof.end()
						&& block[*adjBuilding].buildingCells.size()==1 
						&& block[nextBld].buildingCells.at(0).roofType == block[*adjBuilding].buildingCells.at(0).roofType){
					if(block[nextBld].buildingCells.at(0).roofType == Gabled){
						if(angleDifferenceVectors(block[nextBld].buildingCells.at(0).roofDirection,block[*adjBuilding].buildingCells.at(0).roofDirection) < 5 ){
							sameRoof.insert(*adjBuilding);
							processList.push_back(*adjBuilding);
						}
					}else{
						sameRoof.insert(*adjBuilding);
						processList.push_back(*adjBuilding);
					}
				}
			}
		}
	}
	adjustBuildingsWithSameRoof(block,sameRoof);
	return sameRoof;
}

void adjustBuildingBlock(map<long,Building> &block){
	set<long> processList;
	for(map<long,Building>::iterator building = block.begin();building!=block.end();building++){
		processList.insert(building->first);
	}

	while(!processList.empty()){
		long nextBld = *(processList.begin());
		set<long> processedList = findBuildingsWithSameRoof(block, nextBld);
		for(set<long>::iterator processed=processedList.begin();processed!=processedList.end();processed++){
			processList.erase(*processed);
		}
	}
}