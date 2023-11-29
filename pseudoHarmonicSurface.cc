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
    // Sample surface:
    for (size_t i = 0; i < resolution; ++i) {
        double u = (double)i / (double)(resolution - 1);
        auto x = Geometry::Point2D((1.0 - u) * min[0] + u * max[0], 0);
        qInfo() << x[0] << ", " << x[1];
        auto intersections = helperSurface.findLineCurveIntersections(x, Geometry::Vector2D(std::cos(M_PI * 0.499), std::sin(M_PI * 0.499)));
        qInfo() << intersections.size();
        if (intersections.size() % 2 != 0) {
            continue;
        }
        for (int j = 0; j < (int)intersections.size() - 1; j += 2) {
            double sub_min_y = intersections[j][1];
            double sub_max_y = intersections[j + 1][1];
            int sub_resolution = resolution / intersections.size() * 2;
            for (size_t k = 0; k < sub_resolution; ++k) {
                double v = (double)k / (double)(sub_resolution - 1);
                x[1] = (1.0 - v) * sub_min_y + v * sub_max_y;
                double height = helperSurface.eval(x);
                auto vhd = mesh.add_vertex(BaseTraits::Point(x[0], x[1], height));
                handles.push_back(vhd);
            }
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
