#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class ConePlacer {
  public:
    ConePlacer(ManifoldSurfaceMesh& mesh_, IntrinsicGeometryInterface& geom);

    std::vector<std::pair<Vertex, double>> optimalCones();

    Vector<double> DRminimize(const Vector<double>& initialS);

    Vector<double> stepU(const Vector<double>& s);
    Vector<double> stepV(const Vector<double>& s, const Vector<double>& u);
    Vector<double> stepS(const Vector<double>& s, const Vector<double>& u,
                         const Vector<double>& v);

  public:
    ManifoldSurfaceMesh& mesh;
    IntrinsicGeometryInterface& geom;

    Vector<double> weights;

    Vector<double> k_orig, k_orig_i, k_orig_b;
    SparseMatrix<double> L, L_ii, L_ib;
    SparseMatrix<double> rA, rA_ii, rA_inv, rA_inv_ii; // square root area
    std::unique_ptr<PositiveDefiniteSolver<double>> L_solver, L_ii_solver;
    BlockDecompositionResult<double> Ldecomp;

    VertexData<size_t> vIdx;
    Vector<bool> isInterior;
    size_t nInterior, nBoundary, nVertices;

    double sigma       = 5;
    double alpha       = 1;
    double minAlpha    = 2e-10;
    double eta         = 0;
    double rho         = 0.5;
    double xi          = 1;
    double epsilon     = 1e-3;
    double epsilon_min = 1e-10;
};
