#pragma once

#include "BaseDecoder.h"
#include "Domain.h"
#include "CRC.h"

class ScFlipProgDecoder : public BaseDecoder {

private:
	std::vector<std::vector<double>> _beliefTree;
	std::vector<std::vector<int>> _uhatTree;
	size_t _treeHeight;
	std::vector<int> _maskWithCrc;
	std::vector<int> _x;
	CRC * _crcPtr;
	std::vector<double> _subchannelsMeansGa;

	double f(double llr1, double llr2);
	double g(double llr1, double llr2, int u1);
	int L(double llr1);
	size_t log2(size_t n);

	void FillLeftMessageInTree(std::vector<double>::iterator leftIt,
		std::vector<double>::iterator rightIt,
		std::vector<double>::iterator outIt,
		size_t n);
	void FillRightMessageInTree(std::vector<double>::iterator leftIt,
		std::vector<double>::iterator rightIt,
		std::vector<int>::iterator uhatIt,
		std::vector<double>::iterator outIt,
		size_t n);

	void PassDown(size_t iter);
	void PassUp(size_t iter);

	void DecodeFrom(int position);
	bool IsCrcPassed(std::vector<int> codeword);

	std::vector<int> GetCriticalSet(std::vector<int> mask, int position);
	std::vector<int> SortCriticalBits(std::vector<int> criticalSet, std::vector<double> llrs);

public:
	ScFlipProgDecoder(PolarCode * code);

	std::vector<int> Decode(std::vector<double> llr) override;
	void SetSigma(double sigma) override;
	domain GetDomain() override;

	~ScFlipProgDecoder() {};
};