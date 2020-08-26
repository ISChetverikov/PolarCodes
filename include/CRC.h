#pragma once 

#include <vector>

class CRC {
private:
	std::vector<int> _poly;

public:
	CRC(std::vector<int> poly);
	std::vector<int> GetSum(std::vector<int> bits);
	bool Check(std::vector<int> bytes, std::vector<int> crc);
};