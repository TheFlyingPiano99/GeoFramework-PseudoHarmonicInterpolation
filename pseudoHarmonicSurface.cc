#include <fstream>

#include "pseudoHarmonicSurface.hh"
#include <QDebug>

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

    size_t resolution = 50;

    // Clear previous data:
    mesh.clear();
    std::vector<BaseMesh::VertexHandle> handles, tri;

    // Create helper surface:
    std::function<Geometry::Point2D(double)> f1 = [this](double t) {
        auto p = this->curve(t);
        return Geometry::Point2D(p[0], p[1]);
    };
    std::function<double(Geometry::Point2D)> f2 = [this](Geometry::Point2D x) {
        return this->height(x[0], x[1]);
        };
    auto helperSurface = Geometry::ModifiedGordonWixomSurface(f1, f2);

    // Sample surface:
    for (size_t i = 0; i < resolution; ++i) {
        double u = (double)i / (double)(resolution - 1);
        for (size_t j = 0; j < resolution; ++j) {
            double v = (double)j / (double)(resolution - 1);
            auto x = Geometry::Point2D(u, v);
            double height = helperSurface.eval(x);
            auto vhd = mesh.add_vertex(BaseTraits::Point(x[0], x[1], height));
            handles.push_back(vhd);
        }
    }
    // Build triangles:
    for (size_t i = 0; i < resolution - 1; ++i)
        for (size_t j = 0; j < resolution - 1; ++j) {
            tri.clear();
            tri.push_back(handles[i * resolution + j]);
            tri.push_back(handles[i * resolution + j + 1]);
            tri.push_back(handles[(i + 1) * resolution + j]);
            mesh.add_face(tri);
            tri.clear();
            tri.push_back(handles[(i + 1) * resolution + j]);
            tri.push_back(handles[i * resolution + j + 1]);
            tri.push_back(handles[(i + 1) * resolution + j + 1]);
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
