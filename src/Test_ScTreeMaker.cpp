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


int main3() {
	std::vector<double> output = { 0.010119, 0.269376, 0.714051, 0.025632, 0.421688, 0.062596, 0.928681, 0.009483, 0.929548, 0.969867, 0.948404, 0.005403, 0.266173, 0.957171, 0.011451, 0.715712, };
	auto code = new PolarCode(4, 7, { 0, 1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14, 15 }, {});
	auto decoderPtr = new ScDecoderTreeMaker(code, 1.0);
	decoderPtr->Decode(output);
	auto fullTree = decoderPtr->GetPathInfo();
	auto filename = "fullTree.debug";
	DumpInfo1(filename, VecToStr1(output));
	DumpInfo1(filename, fullTree);
	return 0;
}