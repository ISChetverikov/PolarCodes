#include <vector>
#include <fstream>
#include <string>

#include "../include/ScDecoderTreeMaker.h"
#include "../include/PolarCode.h"


void DumpInfo1(std::string filename, std::string info) {
	std::ofstream resultsFileStream;
	resultsFileStream.open(filename, std::fstream::out | std::fstream::app);
	resultsFileStream << info << std::endl;

	resultsFileStream.close();
}

std::string VecToStr1(std::vector<double> vec) {
	std::string result;

	for (size_t i = 0; i < vec.size(); i++)
	{
		result += std::to_string(vec[i]) + ", ";
	}

	return result;
}


int main5() {
	std::vector<double> output = { 0.070364, 0.884861, 0.831962, 0.013296, 0.987174, 0.584293, 0.677645, 0.134929 };
	auto code = new PolarCode(3, 4, { 0, 1, 2, 4, 3,5,6,7 }, {});
	auto decoderPtr = new ScDecoderTreeMaker(code, 1.0);
	decoderPtr->Decode(output);
	auto fullTree = decoderPtr->GetPathInfo();
	auto filename = "fullTree.debug";
	DumpInfo1(filename, VecToStr1(output));
	DumpInfo1(filename, fullTree);
	return 0;
}