#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct BSpline BSpline;

typedef struct Bspline Bspline;

void free_spline(Bspline *ptr);

const float *knot_domain(BSpline *spline);

BSpline *new_bspline(const float *points, const float *knots, size_t lpts, size_t lkts);

Bspline *new_spline(const float *x, const float *y, size_t len, uint8_t itype);

float sample_spline(Bspline *spline, float t);
