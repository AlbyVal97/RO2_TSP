#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "convex_hull.h"

bool ccw(const Point* a, const Point* b, const Point* c) {
    return (b->x - a->x) * (c->y - a->y)
         > (b->y - a->y) * (c->x - a->x);
}

int comparePoints(const void* lhs, const void* rhs) {
    const Point* lp = lhs;
    const Point* rp = rhs;
    if (lp->x < rp->x)
        return -1;
    if (rp->x < lp->x)
        return 1;
    if (lp->y < rp->y)
        return -1;
    if (rp->y < lp->y)
        return 1;
    return 0;
}

void fatal(const char* message) {
    fprintf(stderr, "%s\n", message);
    exit(1);
}

void* xmalloc(size_t n) {
    void* ptr = malloc(n);
    if (ptr == NULL)
        fatal("Out of memory");
    return ptr;
}

void* xrealloc(void* p, size_t n) {
    void* ptr = realloc(p, n);
    if (ptr == NULL)
        fatal("Out of memory");
    return ptr;
}

void printPoints(const Point* points, int len) {
    printf("[");
    if (len > 0) {
        const Point* ptr = points;
        const Point* end = points + len;
        printf("(#%d: %f, %f)", ptr->id, ptr->x, ptr->y);
        ++ptr;
        for (; ptr < end; ++ptr)
            printf(", (#%d: %f, %f)", ptr->id, ptr->x, ptr->y);
    }
    printf("]");
}

Point* convexHull(Point p[], int len, int* hsize) {
    if (len == 0) {
        *hsize = 0;
        return NULL;
    }

    int i, size = 0, capacity = 4;
    Point* hull = xmalloc(capacity * sizeof(Point));

    qsort(p, len, sizeof(Point), comparePoints);

    /* lower hull */
    for (i = 0; i < len; ++i) {
        while (size >= 2 && !ccw(&hull[size - 2], &hull[size - 1], &p[i]))
            --size;
        if (size == capacity) {
            capacity *= 2;
            hull = xrealloc(hull, capacity * sizeof(Point));
        }
        assert(size >= 0 && size < capacity);
        hull[size++] = p[i];
    }

    /* upper hull */
    int t = size + 1;
    for (i = len - 1; i >= 0; i--) {
        while (size >= t && !ccw(&hull[size - 2], &hull[size - 1], &p[i]))
            --size;
        if (size == capacity) {
            capacity *= 2;
            hull = xrealloc(hull, capacity * sizeof(Point));
        }
        assert(size >= 0 && size < capacity);
        hull[size++] = p[i];
    }
    --size;
    assert(size >= 0);
    hull = xrealloc(hull, size * sizeof(Point));
    *hsize = size;
    return hull;
}