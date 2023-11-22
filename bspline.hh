#pragma once

#include "object.hh"

/*
 * Default B-spline surface
 * Loaded from '*.dbs' files
*/
class BSpline : public Object {
public:
  BSpline(std::string filename);
  virtual ~BSpline();

  // Used only to draw extra elements, such as control points
  virtual void draw(const Visualization &vis) const override;

  // Used to draw clickable things
  virtual void drawWithNames(const Visualization &vis) const override;

  // Returns the position of a selected object
  virtual Vector postSelection(int selected) override;

  // Sets the position of the selected object
  virtual void movement(int selected, const Vector &pos) override;

  // Creates mesh representation of the surface
  virtual void updateBaseMesh() override;

  // Used to load the surface from file
  virtual bool reload() override;

private:
  size_t degree[2];
  size_t no_of_control_points[2];
  std::vector<float> knots[2];
  std::vector<Vector> control_points;

  size_t calculate_cp_index(size_t a, size_t b, size_t nv);
};
