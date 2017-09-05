#include <stdio.h>
#include <math.h>

#include "common.h"
#include "synth.h"

using namespace std;

int main(int argc, char **argv)
{
	FILE *fp = fopen(argv[1], "rb");
	fseek(fp, 14, SEEK_SET);

	int x = 0;
	
	vector<pulse> pulses;
	while (!feof(fp)) {
		int len = getc(fp);
		int cycles;
		if (len == 0) {
			int a = getc(fp);
			int b = getc(fp);
			int c = getc(fp);
			cycles = a | (b << 8) | (c << 16);
		} else {
			cycles = len * 8;
		}
		pulse p;
		p.start = float(x) * WAVE_FREQ / C64_FREQ;
		p.end = (float(x) + cycles * 0.5) * WAVE_FREQ / C64_FREQ;
		pulses.push_back(p);
		x += cycles;
	}

	vector<float> samples = synth(pulses);

	for (unsigned i = 0; i < samples.size(); ++i) {
		//printf("%f %f\n", samples[i], refiltered_samples[i]);
		short s = lrintf(samples[i] * 16384.0f);
		fwrite(&s, 2, 1, stdout);
	}
}
