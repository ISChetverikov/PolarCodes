#pragma once
#include <string>

#include <unordered_map>

enum simulatorType { MC, UnknownSimulation };

simulatorType simulatorFromString(std::string str);
std::string simulationTypeToString(simulatorType);