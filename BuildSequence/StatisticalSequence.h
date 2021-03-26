#pragma once

#include <string>

void BuiltSequenceStatistically(std::string folder, int m, int k, int maxTestsCount, int maxRejectionsCount, double snr);
double PickUpSnr(int m, int k);