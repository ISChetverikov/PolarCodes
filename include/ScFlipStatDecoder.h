#pragma once

#include "ScCrcAidedDecoder.h"

class ScFlipStatDecoder : public ScCrcAidedDecoder {
private:
	std::vector<int> _unfrozenBits;
	std::vector<int> _singleFlipStatistic;
	std::vector<std::vector<int>> _doubleFlipStatistic;

public:
	ScFlipStatDecoder(PolarCode * code);

	std::string GetStatistic() override;
	void ClearStatistic() override;

	std::vector<int> Decode(std::vector<double> llr) override;
	~ScFlipStatDecoder() {};
};