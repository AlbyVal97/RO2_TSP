#ifndef CONVEX_HULL_H_

#define CONVEX_HULL_H_

#include <stdbool.h>

typedef struct tPoint {
	int id;
	double x, y;
} Point;

bool ccw(const Point* a, const Point* b, const Point* c);

int comparePoints(const void* lhs, const void* rhs);

void fatal(const char* message);

void* xmalloc(size_t n);

void* xrealloc(void* p, size_t n);

void printPoints(const Point* points, int len);

Point* convexHull(Point p[], int len, int* hsize);

#endif   /* CONVEX_HULL_H_ */