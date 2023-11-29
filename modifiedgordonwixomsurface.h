#pragma once

#include "geometry.hh"
#include <functional>

namespace Geometry {

  class ModifiedGordonWixomSurface
  {
  public:
    /*
     * Receives a function: t in [0, 1] -> R^3 describing a closed curve
     * The surface will interpolated inside the closed curve
    */
    ModifiedGordonWixomSurface(const std::function<Point2D(double)>& curve, const std::function<double(Point2D)>& height);

    double eval(const Point2D& x);

    void setCurve(const std::function<Point2D(double)>& curve);

    void setHeight(const std::function<double(Point2D)>& curve);

    /*
     * Returns an array of intersection points
    */
    std::vector<Point2D> findLineCurveIntersections(const Point2D& x, const Vector2D& direction);

    Point2D getBoundingRectangleMin() const;

    Point2D getBoundingRectangleMax() const;

private:
    void discretizeCurve();
    Point2D boundingRectangleMin;
    Point2D boundingRectangleMax;
    std::function<Geometry::Point2D(double)> curve;
    std::function<double(Point2D)> height;
    std::vector<Point2D> discretizedCurve;
  };
}
