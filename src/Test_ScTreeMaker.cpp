#include <vector>
#include <fstream>
#include <string>

#include "../include/ScDecoderTreeMaker.h"
#include "../include/PolarCode.h"
#include "../include/CommonTransformations.h"


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


int main0() {
	std::vector<double> output = { 0.91975, 0.00149736, 0.998306, 0.999511, 0.958257, 0.954386, 0.296819, 0.934197, 0.0126016, 0.432695, 0.236357, 0.332338, 7.6258e-05, 0.137163, 0.253909, 0.897426 };
	/*for (size_t i = 0; i < output.size(); i++)
	{
		output[i] = LlrToP1(output[i]);
	}*/
	auto code = new PolarCode(4, 8, { 0, 1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14, 15 }, {});
	auto decoderPtr = new ScDecoderTreeMaker(code, 1.0);
	decoderPtr->Decode(output);
	auto fullTree = decoderPtr->GetPathInfo();
	auto filename = "fullTree.debug";
	DumpInfo1(filename, VecToStr1(output));
	DumpInfo1(filename, fullTree);
	return 0;
}