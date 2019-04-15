#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct Bspline Bspline;

void free_spline(Bspline *ptr);

Bspline *new_spline(const float *x, const float *y, size_t len, uint8_t itype);

float sample_spline(Bspline *spline, float t);
