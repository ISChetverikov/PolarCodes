#pragma once

#include <vector>
#include <tuple>

// for (1024, 512) code!!!!
class FlipStatistic {
private:
	static const std::vector<int> _singles;
	static const std::vector<std::tuple<int, int>> _pairs;

public:
	static std::vector<int> GetSingles() {
		return _singles;
	}

	static std::vector<std::tuple<int, int>> GetPairs() {
		return _pairs;
	}

};

const std::vector<int> FlipStatistic::_singles = { 335, 343, 543, 347, 559, 364, 426, 372, 451, 376 };//, 238, 318, 453, 454, 480, 349, 662, 237, 661, 428, 567, 370, 433, 434, 317, 659, 243, 425, 675, 399, 350, 618, 359, 665, 666, 614, 245, 617, 457, 458, 834, 654, 315, 436, 620, 836, 571, 412, 422, 235, 440, 833, 613, 246, 840, 779, 677, 898, 625, 775, 668, 460, 808, 591, 781, 897, 787, 816, 678, 363, 848, 626, 782, 720, 604, 573, 365, 900, 628, 712, 736, 222, 804, 407, 250, 423, 371, 249, 466, 574, 904, 465, 366, 864, 789, 632, 319, 481, 681, 223, 599, 252, 790, 413, 411, 682, 468, 482, 912, 373, 655, 713, 472, 615, 689, 414, 806, 684, 805, 710, 427, 803, 374, 794, 377, 707, 690, 835, 484, 928, 793, 605, 603, 796, 809, 714, 663, 488, 239, 960, 709, 692, 496, 455, 837, 606, 619, 721, 810, 722, 696, 817, 812, 838, 378, 429, 716, 435, 667, 627, 437, 461, 724, 467, 438, 737, 738, 899, 849 };
const std::vector<std::tuple<int, int>> FlipStatistic::_pairs = {
{335, 343},
{343, 347},
{347, 349},
{238, 237},
{343, 349},
{451, 453},
{335, 347},
{543, 559},
{335, 364},
{347, 364},
//{453, 454},
//{343, 372},
//{343, 318},
//{543, 453},
//{426, 480},
//{662, 661},
//{543, 454},
//{543, 480},
//{347, 376},
//{372, 376},
//{454, 480},
//{426, 428},
//{335, 238},
//{343, 454},
//{347, 372},
//{559, 318},
//{559, 454},
//{364, 372},
//{426, 451},
//{451, 454},
//{335, 426},
//{335, 376},
//{335, 318},
//{335, 661},
//{343, 238},
//{343, 480},
//{343, 428},
//{347, 451},
//{559, 451},
//{364, 376},
//{426, 376},
//{426, 318},
//{376, 318},
//{453, 480},
//{335, 543},
//{335, 559},
//{335, 372},
//{335, 480},
//{335, 349},
//{335, 428},
//{343, 543},
//{343, 559},
//{343, 451},
//{343, 237},
//{343, 661},
//{543, 426},
//{543, 349},
//{543, 237},
//{347, 318},
//{559, 453},
//{559, 662},
//{559, 237},
//{559, 428},
//{364, 451},
//{364, 453},
//{364, 349},
//{426, 454},
//{426, 349},
//{372, 454},
//{451, 318},
//{451, 480},
//{376, 453},
//{376, 349},
//{238, 662},
//{318, 428},
//{453, 428},
//{335, 451},
//{335, 237},
//{343, 426},
//{343, 453},
//{543, 451},
//{543, 318},
//{347, 426},
//{347, 238},
//{347, 453},
//{559, 364},
//{364, 318},
//{364, 480},
//{364, 662},
//{364, 237},
//{426, 372},
//{426, 453},
//{426, 662},
//{426, 237},
//{372, 451},
//{372, 318},
//{372, 349},
//{372, 662},
//{372, 237},
//{451, 376},
//{238, 318},
//{238, 349},
//{318, 237},
//{454, 237},
//{480, 428},
//{661, 428},
//{335, 453},
//{335, 662},
//{343, 364},
//{343, 376},
//{543, 347},
//{543, 364},
//{543, 372},
//{543, 238},
//{543, 662},
//{543, 661},
//{347, 559},
//{347, 454},
//{347, 480},
//{347, 237},
//{347, 428},
//{559, 426},
//{559, 238},
//{559, 480},
//{559, 349},
//{364, 426},
//{364, 238},
//{364, 454},
//{364, 428},
//{372, 661},
//{372, 428},
//{451, 238},
//{451, 349},
//{451, 662},
//{451, 428},
//{376, 480},
//{376, 662},
//{376, 237},
//{238, 453},
//{238, 480},
//{318, 453},
//{318, 661},
//{453, 662},
//{454, 428},
//{480, 349},
//{480, 662},
//{349, 661},
//{349, 428},
//{662, 237},
//{237, 428}
};
