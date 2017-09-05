#ifndef _FILTER_H
#define _FILTER_H 1

class Filter {
public:
	Filter(float a0, float a1, float a2, float b0, float b1, float b2)
		: a1(a1 / a0), a2(a2 / a0), b0(b0 / a0), b1(b1 / a0), b2(b2 / a0), d0(0.0f), d1(0.0f) {}

	inline float update(float in)
	{
		float out = b0*in + d0;
		d0 = b1 * in - a1 * out + d1;
		d1 = b2 * in - a2 * out;
		return out;
	}

	void reset();

	static Filter lpf(float cutoff_radians);
	static Filter hpf(float cutoff_radians);

private:
	float a1, a2, b0, b1, b2;
	float d0, d1;
};

#endif  // !defined(_FILTER_H)
