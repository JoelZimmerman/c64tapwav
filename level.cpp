// A small program to try to even out the audio levels within a file
// (essentially a compressor with infinite lookahead).

#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>

#include "filter.h"

// A final scalar to get the audio within approximately the right range.
// Increase to _lower_ overall volume.
#define DAMPENING_FACTOR 5.0

// 6dB/oct per round.
#define FILTER_DEPTH 4

std::vector<float> level_samples(const std::vector<float> &pcm, float min_level, float lpfilter_freq, int sample_rate)
{
	// filter forwards, then backwards (perfect phase filtering)
	std::vector<float> filtered_samples, refiltered_samples, leveled_samples;
	filtered_samples.resize(pcm.size());
	refiltered_samples.resize(pcm.size());
	leveled_samples.resize(pcm.size());

	Filter filter = Filter::lpf(2.0 * M_PI * lpfilter_freq / sample_rate);
	for (unsigned i = 0; i < pcm.size(); ++i) {
		filtered_samples[i] = filter.update(fabs(pcm[i]));
	}
	filter.reset();
	for (unsigned i = pcm.size(); i --> 0; ) {
		refiltered_samples[i] = filter.update(filtered_samples[i]);
	}

	for (int i = 1; i < FILTER_DEPTH; ++i) {
		filter.reset();
		for (unsigned i = 0; i < pcm.size(); ++i) {
			filtered_samples[i] = filter.update(refiltered_samples[i]);
		}
		filter.reset();
		for (unsigned i = pcm.size(); i --> 0; ) {
			refiltered_samples[i] = filter.update(filtered_samples[i]);
		}
	}

	for (unsigned i = 0; i < pcm.size(); ++i) {
		float f = DAMPENING_FACTOR * std::max<float>(refiltered_samples[i], min_level);
		leveled_samples[i] = pcm[i] / f;
	}

	return leveled_samples;
}
