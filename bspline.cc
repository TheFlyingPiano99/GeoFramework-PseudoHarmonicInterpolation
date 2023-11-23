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

    for (size_t i = 0; i < no_of_control_points[0]; i++) {
      glBegin(GL_LINE_STRIP);
      for (size_t j = 0; j < no_of_control_points[1]; j++) {
          const auto &p = control_points[cp_index(i, j)];
        glVertex3dv(p.data());
      }
      glEnd();
    }
    for (size_t j = 0; j < no_of_control_points[1]; j++) {
      glBegin(GL_LINE_STRIP);
      for (size_t i = 0; i < no_of_control_points[0]; i++) {
        const auto &p = control_points[cp_index(i, j)];
        glVertex3dv(p.data());
      }
      glEnd();
    }

    glPointSize(8.0);
    glColor3d(1.0, 0.0, 1.0);
    glBegin(GL_POINTS);
    for (int u = 0; u < no_of_control_points[0]; u++) {
      const auto &p = control_points[cp_index(u, 0)];
      glVertex3dv(p.data());
    }
    for (int v = 1; v < no_of_control_points[1] - 1; v++) {
      const auto &p = control_points[cp_index(0, v)];
      glVertex3dv(p.data());
    }
    for (int v = 1; v < no_of_control_points[1] - 1; v++) {
      const auto &p = control_points[cp_index(no_of_control_points[0] - 1, v)];
      glVertex3dv(p.data());
    }
    for (int u = 0; u < no_of_control_points[0]; u++) {
      const auto &p = control_points[cp_index(u, no_of_control_points[1] - 1)];
      glVertex3dv(p.data());
    }
    glEnd();

    glLineWidth(1.0);
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
  Geometry::BSSurface helperSurface(degree[0], degree[1], knots[0], knots[1], cps);

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

size_t BSpline::cp_index(size_t u, size_t v) const {
  return u * no_of_control_points[1] + v;
}

void BSpline::calculateInnerControlPoints() {
    // Initial position:
    for (int i = 1; i < no_of_control_points[0] - 1; i++) {
      for (int j = 1; j < no_of_control_points[1] - 1; j++) {
      control_points[cp_index(i, j)] =
          0.5 * ((1.0 - gamma(0, i)) * control_points[cp_index(0, j)]
                                              + gamma(0, i) * control_points[cp_index(no_of_control_points[0] - 1, j)]
                                              + (1.0 - gamma(1, j)) * control_points[cp_index(i, 0)]
                                              + gamma(1, j) * control_points[cp_index(i, no_of_control_points[1] - 1)]);
      }
    }

    // Iterative rearrangement:
    double epsilon = 0.01;
    double max_change = 1.0;
    while (max_change > epsilon) {
      std::vector<Vector> prev_cps = control_points;
      max_change = -1.0;
      for (int i = 1; i < no_of_control_points[0] - 1; i++) {
        for (int j = 1; j < no_of_control_points[1] - 1; j++) {
          control_points[cp_index(i, j)]
            = fullness * (
                  delta(i, j, 1, 0) * prev_cps[cp_index(i - 1, j)]
                        + delta(i, j, 0, 0) * prev_cps[cp_index(i + 1, j)]
                          + delta(i, j, 1, 1) * prev_cps[cp_index(i, j - 1)]
                          + delta(i, j, 0, 1) * prev_cps[cp_index(i, j + 1)]
              )
              + (1.0 - 2.0 * fullness) * (
                    delta(i, j, 1, 0) * delta(i, j, 1, 1) * prev_cps[cp_index(i - 1, j - 1)]
                                          + delta(i, j, 1, 0) * delta(i, j, 0, 1) * prev_cps[cp_index(i - 1, j + 1)]
                                          + delta(i, j , 0, 0) * delta(i, j, 1, 1) * prev_cps[cp_index(i + 1, j - 1)]
                                          + delta(i, j, 0, 0) * delta(i, j, 0, 1) * prev_cps[cp_index(i + 1, j + 1)]
              );
            // Find maximal change in current iteration:
            auto diff = control_points[cp_index(i, j)] - prev_cps[cp_index(i, j)];
            double change = diff.length();
            if (max_change < change) {
              max_change = change;
            }
        }
      }
      std::cout << "Max change: " << max_change << std::endl;
    }
}

double BSpline::gamma(unsigned int u_or_v, unsigned int idx) const {
    double sum = 0.0;
    for (int k = 1; k <= degree[u_or_v]; k++) {
      sum += knots[u_or_v][idx + k];
    }
    return sum / (double)degree[u_or_v];
}

double BSpline::delta(unsigned int i, unsigned int j, int sign, unsigned int u_or_v) const {
    return (gamma(u_or_v, ((u_or_v == 0)? i:j) + sign) - gamma(u_or_v, ((u_or_v == 0)? i:j) + sign - 1))
           / (gamma(u_or_v, ((u_or_v == 0)? i:j) + 1) - gamma(u_or_v, ((u_or_v == 0)? i:j) - 1));
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
    for (size_t u = 0; u < no_of_control_points[0]; u++) {
      size_t idx = cp_index(u, 0);
      std::cout << "U = " << u << " Idx = "<< idx << std::endl;
      f >> control_points[idx][0] >> control_points[idx][1] >> control_points[idx][2];
      std::cout << control_points[idx][0] << " " <<  control_points[idx][1] << " " << control_points[idx][2] << std::endl;
    }
    for (size_t v = 1; v < no_of_control_points[1]; v++) {
      size_t idx = cp_index(no_of_control_points[0] - 1, v);
      std::cout << "V = " << v << " Idx = "<< idx << std::endl;
      f >> control_points[idx][0] >> control_points[idx][1] >> control_points[idx][2];
      std::cout << control_points[idx][0] << " " <<  control_points[idx][1] << " " << control_points[idx][2] << std::endl;
    }
    for (int u = no_of_control_points[0] - 2; u >= 0; u--) {
      size_t idx = cp_index(u, no_of_control_points[1] - 1);
      std::cout << "U = " << u << " Idx = "<< idx << std::endl;
      f >> control_points[idx][0] >> control_points[idx][1] >> control_points[idx][2];
      std::cout << control_points[idx][0] << " " <<  control_points[idx][1] << " " << control_points[idx][2] << std::endl;
    }
    for (int v = no_of_control_points[1] - 2; v >= 1; v--) {
      size_t idx = cp_index(0, v);
      std::cout << "V = " << v << " Idx = "<< idx << std::endl;
      f >> control_points[idx][0] >> control_points[idx][1] >> control_points[idx][2];
      std::cout << control_points[idx][0] << " " <<  control_points[idx][1] << " " << control_points[idx][2] << std::endl;
    }

    calculateInnerControlPoints();

  } catch(std::ifstream::failure &) {
    return false;
  }
  updateBaseMesh();
  return true;
}

double BSpline::getFullness() const {
  return fullness;
}

void BSpline::setFullness(const double f) {
  fullness = f;
}
