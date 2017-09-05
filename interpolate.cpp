#include "interpolate.h"

double lanczos_table[LANCZOS_RADIUS * LANCZOS_RESOLUTION];

void make_lanczos_weight_table()
{
	for (int i = 0; i < LANCZOS_RADIUS * LANCZOS_RESOLUTION; ++i) {
		float x = double(i) / LANCZOS_RESOLUTION;
		lanczos_table[i] = lanczos_weight(x);
	}
}
