#include "../include//DecoderType.h"

decoderType decoderTypeFromString(std::string str) {
	std::unordered_map<std::string, decoderType> decoderTypeResolver = {
		{"SC", decoderType::SC},
		{"SCFano", decoderType::SCFano}
	};

	if (decoderTypeResolver.count(str) > 0)
		return decoderTypeResolver[str];

	return UnknownDecoder;
}

std::string decoderTypeToString(decoderType type) {

	std::unordered_map<decoderType, std::string> decoderTypeStringResolver = {
		{decoderType::SC, "SC"},
		{decoderType::SCFano, "SCFano"},
	};

	return decoderTypeStringResolver[type];
}