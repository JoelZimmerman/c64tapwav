#ifndef _SYNTH_H
#define _SYNTH_H 1

#include <vector>

struct pulse {
	// all values in samples
	double start, end;
};

std::vector<float> synth(const std::vector<pulse> &pulses);

#endif // !defined(_SYNTH_H)
