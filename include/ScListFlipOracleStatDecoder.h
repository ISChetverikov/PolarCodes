#pragma once

#include <map>
#include "ScListFlipStatDecoder.h"

class ScListFlipOracleStatDecoder : public ScListFlipStatDecoder {
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
	int GetFirstErrorPositionMetric(std::vector<std::vector<int>> candidates, std::vector<double> metrics, std::vector<int> originalCodeword);
	int GetFirstErrorPositionMaxDistance(std::vector<std::vector<int>> candidates, std::vector<int> originalCodeword);

public:
	ScListFlipOracleStatDecoder(PolarCode * code, int L);

	std::string GetStatistic() override;
	void ClearStatistic() override;

	std::vector<int> Decode(std::vector<double> llr) override;
	~ScListFlipOracleStatDecoder() {};
};