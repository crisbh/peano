#include "../../apparent-horizon-finder.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace exahype2::ccz4::ahf::tests::spacetimes {

  namespace {

    bool approxEqual(double a, double b, double tol) {
      return std::abs(a - b) <= tol;
    }

    // A scatter of off-axis Cartesian probe points away from the
    // punctures/poles.
    const std::vector<Vector3> probes{
      {1.0, 0.5, 0.3},
      {0.7, -0.4, 0.9},
      {-0.6, 0.8, -0.5},
      {2.0, 0.1, -1.2},
      {0.4, 1.3, 0.6}};

    template <class Spacetime>
    void checkMetricInverseIdentity(
      const Spacetime& spacetime,
      const char*      label
    ) {
      for (const Vector3& X : probes) {
        const InterpolatedGeometryData data = MetricData(spacetime, X);
        const Tensor3                  ginv = inverse(data.spatialMetric);
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            double product = 0.0;
            for (int k = 0; k < 3; ++k) {
              product += data.spatialMetric[i][k] * ginv[k][j];
            }
            assert(approxEqual(product, i == j ? 1.0 : 0.0, 1.0e-10));
          }
        }
      }
      std::cout << label << ": gamma gamma^{-1} = I" << std::endl;
    }

  } // namespace

  // Complex-step MetricData yields a metric whose closed-form inverse recovers
  // the identity for every spacetime.
  void testMetricInverseIdentity() {
    checkMetricInverseIdentity(
      SchwarzschildIsotropic{1.0},
      "Schwarzschild (isotropic)"
    );
    checkMetricInverseIdentity(
      SchwarzschildKerrSchild{1.0},
      "Schwarzschild (Kerr-Schild)"
    );
    checkMetricInverseIdentity(KerrKerrSchild{1.0, 0.6}, "Kerr (Kerr-Schild)");
    checkMetricInverseIdentity(
      BrillLindquist::single(1.0),
      "Brill-Lindquist (single)"
    );
  }

  // The two time-symmetric slices have K_ij = 0 exactly through the
  // complex-step path.
  void testTimeSymmetricExtrinsicCurvatureVanishes() {
    const SchwarzschildIsotropic isotropic{1.0};
    const BrillLindquist         puncture = BrillLindquist::single(1.4);
    for (const Vector3& X : probes) {
      for (const InterpolatedGeometryData& data :
           {MetricData(isotropic, X), MetricData(puncture, X)}) {
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            assert(approxEqual(data.extrinsicCurvature[i][j], 0.0, 1.0e-10));
          }
        }
        assert(approxEqual(data.traceK, 0.0, 1.0e-10));
      }
    }
    std::cout << "Time-symmetric slices give K_ij = 0 and trK = 0" << std::endl;
  }

  // The Kerr-Schild slices carry a non-trivial, finite, symmetric K_ij.
  void testKerrSchildExtrinsicCurvatureFiniteSymmetric() {
    const SchwarzschildKerrSchild schwarzschild{1.0};
    const KerrKerrSchild          kerr{1.0, 0.6};
    bool                          sawNonZero = false;
    for (const Vector3& X : probes) {
      for (const InterpolatedGeometryData& data :
           {MetricData(schwarzschild, X), MetricData(kerr, X)}) {
        assert(std::isfinite(data.traceK));
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            assert(std::isfinite(data.extrinsicCurvature[i][j]));
            assert(approxEqual(
              data.extrinsicCurvature[i][j],
              data.extrinsicCurvature[j][i],
              1.0e-12
            ));
            if (std::abs(data.extrinsicCurvature[i][j]) > 1.0e-6) {
              sawNonZero = true;
            }
          }
        }
      }
    }
    assert(sawNonZero);
    std::cout
      << "Kerr-Schild slices give finite, symmetric, non-zero K_ij"
      << std::endl;
  }

  // A coordinate sphere at the known horizon radius makes Theta vanish for the
  // two spherically symmetric spacetimes (the horizon is a coordinate sphere
  // there).
  void testExpansionVanishesOnSphericalHorizon() {
    const ReferencePoint origin{0.0, 0.0, 0.0};
    const SpectralBasis  basis = BuildSpectralBasis(4);

    struct Case {
      const char* label;
      double      radius;
    };
    const SchwarzschildIsotropic  isotropic{1.0};
    const SchwarzschildKerrSchild kerrSchild{1.0};

    {
      const std::vector<double>
        a = SphereCoefficients(basis, isotropic.horizonRadius());
      const InterpolatedSurface surface = BuildAnalyticSurface(
        basis,
        a,
        isotropic,
        origin
      );
      double maxTheta = 0.0;
      for (const InterpolatedSurfacePoint& point : surface.points) {
        maxTheta = std::max(maxTheta, std::abs(computeExpansion(point)));
      }
      assert(maxTheta < 1.0e-9);
      std::cout
        << "Schwarzschild isotropic horizon (r = M/2): max|Theta| = "
        << maxTheta << std::endl;
    }
    {
      const std::vector<double>
        a = SphereCoefficients(basis, kerrSchild.horizonRadius());
      const InterpolatedSurface surface = BuildAnalyticSurface(
        basis,
        a,
        kerrSchild,
        origin
      );
      double maxTheta = 0.0;
      for (const InterpolatedSurfacePoint& point : surface.points) {
        maxTheta = std::max(maxTheta, std::abs(computeExpansion(point)));
      }
      assert(maxTheta < 1.0e-9);
      std::cout
        << "Schwarzschild Kerr-Schild horizon (r = 2M): max|Theta| = "
        << maxTheta << std::endl;
    }
  }

} // namespace exahype2::ccz4::ahf::tests::spacetimes

int main() {
  using namespace exahype2::ccz4::ahf::tests::spacetimes;

  std::cout << "Apparent Horizon Finder spacetime tests begin..." << std::endl;

  testMetricInverseIdentity();                        // complex-step differentiation test
  testTimeSymmetricExtrinsicCurvatureVanishes();      // complext-step differentiation test
  testKerrSchildExtrinsicCurvatureFiniteSymmetric();
  testExpansionVanishesOnSphericalHorizon();

  std::cout << "All spacetime tests completed." << std::endl;
  return 0;
}
