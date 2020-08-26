#include "../include/CRC.h"
#include <vector>
#include <iostream>

int main1() {
	std::vector<int> poly = { 1, 1, 1, 0, 0, 0, 0, 0, 1};
	//std::vector<int> input = { 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<int> input = { 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0 };

	auto c = CRC(poly);
	std::vector<int> hash = c.Calculate(input);
	int length = hash.size();
	for (int i = length - 1; i >= 0; i--)
	{
		std::cout << hash[i];
	}
	std::cout << std::endl;

	return 0;
}