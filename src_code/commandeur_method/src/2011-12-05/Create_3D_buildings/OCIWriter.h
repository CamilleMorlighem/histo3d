#ifndef _BUILDING_H_INCLUDED
#define _BUILDING_H_INCLUDED
#include "Building.h"
#endif

#include <map>

void writeBuildingToOCI(Building building);
void writeBlockToOCI(std::map<long,Building> block);