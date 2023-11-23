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

  double getFullness() const;

  void setFullness(const double f);

private:
  size_t degree[2];
  size_t no_of_control_points[2];
  std::vector<double> knots[2];
  std::vector<Vector> control_points;
  double fullness = 0.5;
  size_t cp_index(size_t a, size_t b) const;

  void calculateInnerControlPoints();

  double gamma(unsigned int u_or_v, unsigned int idx) const;

  /*
   * sign in {0, 1}
   * 0 ... negative
   * 1 ... positive
  */
  double delta(unsigned int i, unsigned int j, int sign, unsigned int u_or_v) const;
};
