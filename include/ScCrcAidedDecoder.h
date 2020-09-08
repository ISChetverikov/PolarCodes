#pragma once

#pragma once

#include "ScDecoder.h"
#include "Domain.h"
#include "CRC.h"

// logically it is abstract class
class ScCrcAidedDecoder : public ScDecoder {

protected:
	
	CRC * _crcPtr;

	void DecodeFrom(int position);
	void DecodeFromTo(int position, int endPosition);
	bool IsCrcPassed(std::vector<int> codeword);

public:
	ScCrcAidedDecoder(PolarCode * code, domain domain, bool isMinSum);
	~ScCrcAidedDecoder() {};
};