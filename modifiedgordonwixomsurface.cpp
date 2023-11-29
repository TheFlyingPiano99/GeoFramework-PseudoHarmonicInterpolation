#include "modifiedgordonwixomsurface.h"
#include <math.h>

Geometry::ModifiedGordonWixomSurface::ModifiedGordonWixomSurface(const std::function<Point2D(double)>& _curve,
                                                                 const std::function<double(Point2D)>& _height)
    : curve(_curve), height(_height)
{
    discretizeCurve();
}

double Geometry::ModifiedGordonWixomSurface::eval(const Point2D &x)
{
    constexpr int n = 100;
    constexpr double delta_theta = 2.0 * M_PI / n;
    double integral_den = 0.0;
    double integral_div = 0.0;
    for (int i = 0; i < n; i++) {
        Vector2D direction(std::cos(i * delta_theta), std::sin(i * delta_theta));
        auto intersections = findLineCurveIntersections(x, direction);
        if (2 > intersections.size()) {
            continue;
        }
        double a = 0.0;
        double b = 0.0;
        double c = 1.0 / (intersections[0] - x).length();
        double d = 0.0;
        for (int i = 0; i < intersections.size(); i++) {
            double distance = (intersections[i] - x).length();
            a += ((i == 0 || i % 2 == 1)? 1.0 : -1.0) * height(intersections[i]) / distance;
            b += ((i == 0 || i % 2 == 1)? 1.0 : -1.0) / distance;
            if (0 < i) {
                d += ((i == 0 || i % 2 == 1)? 1.0 : -1.0) / distance;
            }
        }
        c *= d;
        integral_den += a / b * c * delta_theta;
        integral_div += c * delta_theta;
    }
    double u = integral_den / integral_div;
    return u;
}

void Geometry::ModifiedGordonWixomSurface::setCurve(const std::function<Point2D(double)> &_curve)
{
    curve = _curve;
    discretizeCurve();
}

void Geometry::ModifiedGordonWixomSurface::setHeight(const std::function<double(Point2D)>&_height)
{
    height = _height;
}

void Geometry::ModifiedGordonWixomSurface::discretizeCurve()
{
    discretizedCurve.clear();
    constexpr int n = 128;
    for (int i = 0; i < n; i++) {
        discretizedCurve.push_back(curve(i / (double)n));
    }
}

std::vector<Geometry::Point2D> Geometry::ModifiedGordonWixomSurface::findLineCurveIntersections(
    const Point2D& x, const Vector2D& direction
)
{
    std::vector<Point2D> intersection_points;
    for (int i = 0; i < discretizedCurve.size(); i++){
        Point2D p0 = discretizedCurve[i];
        Point2D p1 = (i == discretizedCurve.size() - 1)? discretizedCurve[0] : discretizedCurve[i + 1];
        Point2D sectionDiff = p1 - p0;
        double sectionLength = sectionDiff.length();
        if (sectionLength < std::numeric_limits<double>::min()) {
            continue;
        }
        Vector2D sectionDir = sectionDiff / sectionLength;
        double t = x[1] / sectionDir[1] + direction[1] / (direction[0] * sectionDir[1])
                                              * (p0[0] - x[0]) - p0[1] / sectionDir[1];
        if (t >= 0 && t < sectionLength) {
            intersection_points.push_back(p0 + sectionDir * t);
        }
    }
    return intersection_points;
}

