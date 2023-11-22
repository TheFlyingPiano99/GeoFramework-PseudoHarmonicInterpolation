class BSSurface {
public:
    // Constructors
    BSSurface();
    BSSurface(const BSSurface &) = default;
    BSSurface(size_t deg_u, size_t deg_v, const PointVector &cpts);
    BSSurface(size_t deg_u, size_t deg_v, const DoubleVector &knots_u, const DoubleVector &knots_v,
              const PointVector &cpts);
    BSSurface &operator=(const BSSurface &) = default;

    // Evaluation
    Point3D eval(double u, double v) const;
    Point3D eval(double u, double v, size_t nr_der, VectorMatrix &der) const;

    // Coordinates
    std::array<size_t, 2> numControlPoints() const;
    Point3D controlPoint(size_t i, size_t j) const;
    Point3D &controlPoint(size_t i, size_t j);
    const PointVector &controlPoints() const;
    PointVector &controlPoints();

    // Parameterization
    const BSBasis &basisU() const;
    const BSBasis &basisV() const;
    void swapUV();
    void reverseU();
    void reverseV();
    void normalize();

    // Algorithms
    BSSurface insertKnotU(double u, size_t r) const;
    BSSurface insertKnotV(double v, size_t r) const;

private:
    BSSurface insertKnotU(double u, size_t k, size_t s, size_t r) const;
    BSSurface insertKnotV(double v, size_t k, size_t s, size_t r) const;

    size_t n_u_, n_v_;
    BSBasis basis_u_, basis_v_;
    PointVector cp_;
};
