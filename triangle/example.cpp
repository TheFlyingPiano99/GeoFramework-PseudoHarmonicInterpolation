#define ANSI_DECLARATORS
#define REAL double
#define VOID void
extern "C" {
#include <triangle.h>
}

// Example: domain_ contains the points of the closed polygon

// Input points
size_t n = domain_.size();
DoubleVector points; points.reserve(2 * n);
for (const auto& p : domain_) {
    points.push_back(p[0]);
    points.push_back(p[1]);
}

// Input segments : just a closed polygon
std::vector<int> segments; segments.reserve(2 * n);
for (size_t i = 0; i < n; ++i) {
    segments.push_back(i);
    segments.push_back(i + 1);
}
segments.back() = 0;

// Setup output data structure
struct triangulateio in, out;
in.pointlist = &points[0];
in.numberofpoints = n;
in.numberofpointattributes = 0;
in.pointmarkerlist = nullptr;
in.segmentlist = &segments[0];
in.numberofsegments = n;
in.segmentmarkerlist = nullptr;
in.numberofholes = 0;
in.numberofregions = 0;

// Setup output data structure
out.pointlist = nullptr;
out.pointattributelist = nullptr;
out.pointmarkerlist = nullptr;
out.trianglelist = nullptr;
out.triangleattributelist = nullptr;
out.segmentlist = nullptr;
out.segmentmarkerlist = nullptr;

// Call the library function [with maximum triangle area = resolution]
std::ostringstream cmd;
cmd << "pqa" << std::fixed << resolution << "DBPzQ";
triangulate(const_cast<char*>(cmd.str().c_str()), &in, &out, (struct triangulateio*)nullptr);

// Now the result is:
// - out.numberofpoints
// - out.pointlist (stored in the format x0,y0,x1,y1,x2,y2,...)
// - out.numberoftriangles
// - out.trianglelist (stored in the format T0a,T0b,T0c,T1a,T1b,T1c,T2a,T2b,T2c,...)
