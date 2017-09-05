#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>

#include "common.h"

using namespace std;

#define LANCZOS_RADIUS 10
#define RESOLUTION 256

// Cutoff frequency of low-pass filter (must be max. WAVE_FREQ / 2 or you
// will get aliasing)
#define LPFILTER_FREQ 22050

// Cutoff frequency of high-pass filter
#define HPFILTER_FREQ 40

#define NORMALIZED_LPFILTER_FREQ (float(LPFILTER_FREQ) / float(WAVE_FREQ))
#define LANCZOS_EFFECTIVE_RADIUS (LANCZOS_RADIUS / (NORMALIZED_LPFILTER_FREQ * 2.0))

#define SIZEOF_HALF_TABLE ((int)((LANCZOS_EFFECTIVE_RADIUS * RESOLUTION)))
static double *integrated_window;

double sinc(double x)
{
	if (fabs(x) < 1e-6) {
		return 1.0f - fabs(x);
	} else {
		return sin(x) / x;
	}
}

double weight(double x)
{
	if (fabs(x) > LANCZOS_EFFECTIVE_RADIUS) {
		return 0.0f;
	}
	float t = 2.0 * M_PI * NORMALIZED_LPFILTER_FREQ * x;
	return sinc(t) * sinc(t / LANCZOS_RADIUS);
}

void make_table()
{
	integrated_window = new double[2 * SIZEOF_HALF_TABLE + 1];

	double sum = 0.0;
	for (int i = 0; i <= 2 * SIZEOF_HALF_TABLE; ++i) {
		float t = (i - SIZEOF_HALF_TABLE) / (float)RESOLUTION;
		integrated_window[i] = sum;
		sum += weight(t) * NORMALIZED_LPFILTER_FREQ * 2.0 / (float)RESOLUTION;
		//printf("%f %f %f\n", t, weight(t), sum);
	}
//	exit(0);
}

// integral from -inf to window
double window_one_sided_integral(double to)
{
	double array_pos = to * RESOLUTION + SIZEOF_HALF_TABLE;

	int whole = int(floor(array_pos));
	double frac = array_pos - whole;
	if (whole < 0) {
		return 0.0;
	}
	if (whole >= 2 * SIZEOF_HALF_TABLE) {
		return 1.0;
	}
	return integrated_window[whole] + frac * (integrated_window[whole + 1] - integrated_window[whole]);
}

double window_integral(double from, double to)
{
	return window_one_sided_integral(to) - window_one_sided_integral(from);
}

struct pulse {
	// all values in samples
	double start, end;
};
vector<pulse> pulses;

static float a1, a2, b0, b1, b2;
static float d0, d1;

static void filter_init(float cutoff_radians)
{
	float resonance = 1.0f / sqrt(2.0f);
	float sn, cs;
	sincosf(cutoff_radians, &sn, &cs);

	float alpha = float(sn / (2 * resonance));

	// coefficients for highpass filter
        float a0 = 1 + alpha;
	b0 = (1 + cs) * 0.5f;
	b1 = -(1 + cs);
	b2 = b0;
        a1 = -2 * cs;
        a2 = 1 - alpha;

	float invA0 = 1.0f / a0;
	b0 *= invA0;
	b1 *= invA0;
	b2 *= invA0;
	a1 *= invA0;
	a2 *= invA0;

	// reset filter delays
	d0 = d1 = 0.0f;
}

static float filter_update(float in)
{
	float out = b0*in + d0;
	d0 = b1 * in - a1 * out + d1;
	d1 = b2 * in - a2 * out;
	return out;
}

vector<float> synth(const vector<pulse> &pulses)
{
	make_table();

	int len_samples = int(ceil(pulses.back().start + (pulses.back().end - pulses.back().start) * 2));
	
	fprintf(stderr, "%d pulses, total %.2f seconds (%d samples)\n", int(pulses.size()), len_samples / float(WAVE_FREQ), len_samples);

	int pulse_begin = 0;

	vector<float> samples;
	samples.reserve(len_samples);

	for (int i = 0; i < len_samples; ++i) {
//	for (int i = 50000000; i < 51000000; ++i) {
//		double t = i * 0.01;
		double t = i;
		double sample = -0.5;
		for (unsigned j = pulse_begin; j < pulses.size(); ++j) {
			if (t - pulses[j].end > LANCZOS_EFFECTIVE_RADIUS) {
				++pulse_begin;
				continue;
			}
			if (pulses[j].start - t > LANCZOS_EFFECTIVE_RADIUS) {
				break;
			}
			float contribution = window_integral(pulses[j].start - t, pulses[j].end - t);
			sample += contribution;
		}
		samples.push_back(sample);
	}

	// filter forwards, then backwards (perfect phase filtering)
	vector<float> filtered_samples, refiltered_samples;
	filtered_samples.resize(len_samples);
	refiltered_samples.resize(len_samples);

	filter_init(M_PI * HPFILTER_FREQ / WAVE_FREQ);
	for (int i = 0; i < len_samples; ++i) {
		filtered_samples[i] = filter_update(samples[i]);
	}
	filter_init(M_PI * HPFILTER_FREQ / WAVE_FREQ);
	for (int i = len_samples; i --> 0; ) {
		refiltered_samples[i] = filter_update(filtered_samples[i]);
	}

	return refiltered_samples;
}
