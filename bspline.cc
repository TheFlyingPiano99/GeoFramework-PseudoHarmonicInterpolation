#include <fstream>

#include "bspline.hh"
#include "bsplineHelper.hh"

BSpline::BSpline(std::string filename) : Object(filename) {
  reload();
}

BSpline::~BSpline() {
}

void BSpline::draw(const Visualization &vis) const {

  Object::draw(vis);
  if (vis.show_control_points) {
    glDisable(GL_LIGHTING);
    glLineWidth(3.0);
    glColor3d(0.3, 0.3, 1.0);
    size_t m = degree[1] + 1;
    for (size_t k = 0; k < 2; ++k)
      for (size_t i = 0; i <= degree[k]; ++i) {
        glBegin(GL_LINE_STRIP);
        for (size_t j = 0; j <= degree[1-k]; ++j) {
          size_t const index = k ? j * m + i : i * m + j;
          const auto &p = control_points[index];
          glVertex3dv(p.data());
        }
        glEnd();
      }
    glLineWidth(1.0);
    glPointSize(8.0);
    glColor3d(1.0, 0.0, 1.0);
    glBegin(GL_POINTS);
    for (const auto &p : control_points)
      glVertex3dv(p.data());
    glEnd();
    glPointSize(1.0);
    glEnable(GL_LIGHTING);
  }
}

void BSpline::drawWithNames(const Visualization &vis) const {
  if (!vis.show_control_points)
    return;
  for (size_t i = 0; i < control_points.size(); ++i) {
    const auto &p = control_points[i];
    glPushName(i);
    glRasterPos3dv(p.data());
    glPopName();
  }
}

Vector BSpline::postSelection(int selected) {
  return control_points[selected];
}

void BSpline::movement(int selected, const Vector &pos) {
  control_points[selected] = pos;
}

static void bernstein(size_t n, double u, std::vector<double> &coeff) {
  coeff.clear(); coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
}

void BSpline::updateBaseMesh() {
  // TODO: rewrite for B-spline

  size_t resolution = 50;

  mesh.clear();
  std::vector<BaseMesh::VertexHandle> handles, tri;

  std::vector<Geometry::Vector3D> cps;
  for (int i = 0; i < control_points.size(); i++) {
    auto p = control_points[i];
    cps.push_back(Geometry::Vector3D(p[0], p[1], p[2]));
  }
  Geometry::BSSurface helperSurface(degree[0], degree[1], cps);

  // Sample surface:
  for (size_t i = 0; i < resolution; ++i) {
    double u = (double)i / (double)(resolution - 1);
    for (size_t j = 0; j < resolution; ++j) {
      double v = (double)j / (double)(resolution - 1);
      auto p = helperSurface.eval(u, v);
      BaseTraits::Point mp;
      mp[0] = p[0];
      mp[1] = p[1];
      mp[2] = p[2];
      handles.push_back(mesh.add_vertex(mp));
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
}

size_t BSpline::calculate_cp_index(size_t u, size_t v, size_t nv) {
  return u * nv + v;
}

bool BSpline::reload() {
  try {
    std::ifstream f(filename);
    f.exceptions(std::ios::failbit | std::ios::badbit);
    f >> degree[0] >> degree[1];
    std::cout << "U degree: " << degree[0] << " V degree: " << degree[1] << std::endl;
    f >> no_of_control_points[0] >> no_of_control_points[1];
    std::cout << "U No. of cps: " << no_of_control_points[0] << " V No. of cps: " << no_of_control_points[1] << std::endl;
    control_points.resize(no_of_control_points[0] * no_of_control_points[1]);
    for (size_t i = 0; i <= no_of_control_points[0] + degree[0]; i++) {
      float knot;
      f >> knot;
      knots[0].push_back(knot);
      std::cout << "U knot: " << knot << " ";
    }
    for (size_t i = 0; i <= no_of_control_points[1] + degree[1]; i++) {
      float knot;
      f >> knot;
      knots[1].push_back(knot);
      std::cout << "V knot: " << knot << " ";
    }

    control_points.resize(
        no_of_control_points[0] + no_of_control_points[1]
        - 1 + no_of_control_points[0]
        - 2 + no_of_control_points[1] - 3
    );
    std::cout << "\nAfter Resize" << std::endl;
    for (size_t u = 0; u < no_of_control_points[0]; u++) {
      size_t cp_index = calculate_cp_index(u, 0, no_of_control_points[1]);
      f >> control_points[cp_index][0] >> control_points[cp_index][1] >> control_points[cp_index][2];
      std::cout << control_points[cp_index][0] << " " <<  control_points[cp_index][1] << " " << control_points[cp_index][2] << std::endl;
    }
    std::cout << "After Resize" << std::endl;
    for (size_t v = 1; v < no_of_control_points[1]; v++) {
      size_t cp_index = calculate_cp_index(no_of_control_points[0] - 1, v, no_of_control_points[1]);
      f >> control_points[cp_index][0] >> control_points[cp_index][1] >> control_points[cp_index][2];
      std::cout << control_points[cp_index][0] << " " <<  control_points[cp_index][1] << " " << control_points[cp_index][2] << std::endl;
    }
    std::cout << "After Resize" << std::endl;
    for (size_t u = no_of_control_points[0] - 2; u >= 0; u--) {
      size_t cp_index = calculate_cp_index(u, no_of_control_points[1] - 1, no_of_control_points[1]);
      f >> control_points[cp_index][0] >> control_points[cp_index][1] >> control_points[cp_index][2];
      std::cout << control_points[cp_index][0] << " " <<  control_points[cp_index][1] << " " << control_points[cp_index][2] << std::endl;
    }
    std::cout << "After Resize" << std::endl;
    for (size_t v = no_of_control_points[1] - 2; v >= 1; v--) {
      size_t cp_index = calculate_cp_index(0, v, no_of_control_points[1]);
      f >> control_points[cp_index][0] >> control_points[cp_index][1] >> control_points[cp_index][2];
    }
    std::cout << "After Resize" << std::endl;

  } catch(std::ifstream::failure &) {
    return false;
  }
  updateBaseMesh();
  return true;
}
