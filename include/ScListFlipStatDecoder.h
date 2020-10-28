#pragma once

#include "ScListDecoder.h"

class ScListFlipStatDecoder : public ScListDecoder {
private:
	std::vector<int> _singleFlips;
	std::vector<std::tuple<int, int>> _doubleFlips;

protected:
	void DecodeFlipListInternal(std::vector<double> inLlr, std::vector<int> flipPositions);
	std::vector<int> TakeListStatResult(bool& isError);
public:
	ScListFlipStatDecoder(PolarCode * code, int L);

	std::vector<int> Decode(std::vector<double> llr) override;
	~ScListFlipStatDecoder() {};
};