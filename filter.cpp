#include <math.h>

#include "filter.h"

void Filter::reset()
{
	d0 = d1 = 0.0f;
}

Filter Filter::lpf(float cutoff_radians)
{
	float resonance = 1.0f / sqrt(2.0f);
	float sn = sin(cutoff_radians), cs = cos(cutoff_radians);
	float alpha = float(sn / (2 * resonance));

	// coefficients for lowpass filter
        float a0 = 1 + alpha;
	float b0 = (1 - cs) * 0.5f;
	float b1 = 1 - cs;
	float b2 = b0;
        float a1 = -2 * cs;
        float a2 = 1 - alpha;

	return Filter(a0, a1, a2, b0, b1, b2);
}

Filter Filter::hpf(float cutoff_radians)
{
	float resonance = 1.0f / sqrt(2.0f);
	float sn = sin(cutoff_radians), cs = cos(cutoff_radians);
	float alpha = float(sn / (2 * resonance));

	// coefficients for highpass filter
	float a0 = 1 + alpha;
	float b0 = (1 + cs) * 0.5f;
	float b1 = -(1 + cs);
	float b2 = b0;
	float a1 = -2 * cs;
	float a2 = 1 - alpha;

	return Filter(a0, a1, a2, b0, b1, b2);
}
