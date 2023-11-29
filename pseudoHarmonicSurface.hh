#pragma once

#include "object.hh"
#include <functional>
#include "modifiedgordonwixomsurface.h"


class PseudoHarmonicSurface : public Object {
public:
    PseudoHarmonicSurface(
        const std::function<BaseTraits::Point(double t)>& curve,
        const std::function<double(double x, double y)>& height
    );

  virtual ~PseudoHarmonicSurface();

  // Used only to draw extra elements, such as control points
  virtual void draw(const Visualization &vis) const override;

  // Used to draw clickable things
  virtual void drawWithNames(const Visualization &vis) const override;

  // Creates mesh representation of the surface
  virtual void updateBaseMesh() override;

  // Used to load the surface from file
  virtual bool reload() override;

  // Returns the position of a selected object
  virtual Vector postSelection(int selected) override {}

  // Sets the position of the selected object
  virtual void movement(int selected, const Vector &pos) override {}

  void setCurve(const std::function<BaseTraits::Point(double)>& func);


private:

  std::function<BaseTraits::Point(double t)> curve; // t in [0, 1] -> R^2 (z coordinate is ignored)
  std::function<double(double x, double y)> height;
};
