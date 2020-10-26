#pragma once

#include <map>
#include "ScListDecoder.h"

class ScListFlipOracleStatDecoder : public ScListDecoder {
private:
	std::map<int, int> _singleOracleFlipsStat;
	std::map<std::tuple<int, int>, int> _doubleOracleFlipsStat;
	std::map<std::tuple<int, int, int>, int> _tripleOracleFlipsStat;

	int _count = 0;
	int _countErrorneous = 0;
	int _singleSuccessfulFlips = 0;
	int _doubleSuccessfulFlips = 0;
	int _tripleSuccessfulFlips = 0;
	
protected:

	// Returns exact flips by comparision with _codeword at each step of list decoding
	std::vector<int> DecodeFlipOracleListInternal(std::vector<double> inLlr);

public:
	ScListFlipOracleStatDecoder(PolarCode * code, int L);

	std::string GetStatistic() override;
	void ClearStatistic() override;

	std::vector<int> Decode(std::vector<double> llr) override;
	~ScListFlipOracleStatDecoder() {};
};