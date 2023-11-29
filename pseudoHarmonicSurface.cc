#include <fstream>

#include "pseudoHarmonicSurface.hh"
#include <QDebug>
#define ANSI_DECLARATORS
#define REAL double
#define VOID void
extern "C" {
#include <triangle/triangle.h>
}


PseudoHarmonicSurface::PseudoHarmonicSurface(
    const std::function<BaseTraits::Point(double t)>& curve,
    const std::function<double(double x, double y)>& height
)
    : Object(""), curve(curve), height(height) {
    this->reload();
    qInfo() << "Finished\n";
}

PseudoHarmonicSurface::~PseudoHarmonicSurface() {
}

void PseudoHarmonicSurface::draw(const Visualization &vis) const {
  Object::draw(vis);
  // No extra elements are drawn.
}

void PseudoHarmonicSurface::drawWithNames(const Visualization &vis) const {
  // No movable control points are drawn
}


void PseudoHarmonicSurface::updateBaseMesh() {

    size_t resolution = 60;

    // Clear previous data:
    mesh.clear();
    std::vector<BaseMesh::VertexHandle> handles, tri;

    // Create helper surface:
    std::function<Geometry::Point2D(double)> f_c = [this](double t) {
        auto p = this->curve(t);
        return Geometry::Point2D(p[0], p[1]);
    };
    std::function<double(Geometry::Point2D)> f_h = [this](Geometry::Point2D x) {
        return this->height(x[0], x[1]);
        };
    auto helperSurface = Geometry::ModifiedGordonWixomSurface(f_c, f_h);
    Geometry::Point2D min = helperSurface.getBoundingRectangleMin();
    Geometry::Point2D max = helperSurface.getBoundingRectangleMax();

    const auto& discretizedCurve = helperSurface.getDiscretizedCurve();
    // Input points
    size_t n = discretizedCurve.size();
    Geometry::DoubleVector points; points.reserve(2 * n);
    for (const auto& p : discretizedCurve) {
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
    in.edgelist = nullptr;
    in.edgemarkerlist = nullptr;
    in.holelist = nullptr;
    in.neighborlist = nullptr;
    in.normlist = nullptr;
    in.numberofcorners = 0;
    in.numberofedges = 0;
    in.pointlist = &points[0];
    in.numberofpoints = n;
    in.numberofpointattributes = 0;
    in.pointmarkerlist = nullptr;
    in.segmentlist = &segments[0];
    in.numberofsegments = n;
    in.segmentmarkerlist = nullptr;
    in.numberofholes = 0;
    in.numberofregions = 0;
    in.numberoftriangleattributes = 0;
    in.numberoftriangles = 0;
    in.trianglelist = nullptr;

    // Setup output data structure
    out.edgelist = nullptr;
    out.edgemarkerlist = nullptr;
    out.holelist = nullptr;
    out.neighborlist = nullptr;
    out.normlist = nullptr;
    out.numberofcorners = 0;
    out.numberofedges = 0;
    out.numberofholes = 0;
    out.numberofregions = 0;
    out.numberofsegments = 0;
    out.numberofpointattributes = 0;
    out.numberofpoints = 0;
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

    qInfo() << "Triangulation started.";
    triangulate(const_cast<char*>(cmd.str().c_str()), &in, &out, (struct triangulateio*)nullptr);
    qInfo() << "Triangulation finished.";

    // Now the result is:
    // - out.numberofpoints
    // - out.pointlist (stored in the format x0,y0,x1,y1,x2,y2,...)
    // - out.numberoftriangles
    // - out.trianglelist (stored in the format T0a,T0b,T0c,T1a,T1b,T1c,T2a,T2b,T2c,...)


    // Sample surface:
    for (size_t i = 0; i < out.numberofpoints; i++) {
        auto x = Geometry::Point2D(out.pointlist[i * 2], out.pointlist[i * 2 + 1]);
        double height = helperSurface.eval(x);
        auto vhd = mesh.add_vertex(BaseTraits::Point(x[0], x[1], height));
        handles.push_back(vhd);
    }

    // Build triangles:
    for (size_t i = 0; i < out.numberoftriangles; i++) {
            tri.clear();
            tri.push_back(handles[out.trianglelist[i * 3]]);
            tri.push_back(handles[out.trianglelist[i * 3 + 1]]);
            tri.push_back(handles[out.trianglelist[i * 3 + 2]]);
            mesh.add_face(tri);
    }
    Object::updateBaseMesh(false, false);
    qInfo() << "Finished initialization.\n";
}

bool PseudoHarmonicSurface::reload() {
    this->updateBaseMesh();
    return true;
}

void PseudoHarmonicSurface::setCurve(const std::function<BaseTraits::Point (double)> &func)
{
    curve = func;
}
