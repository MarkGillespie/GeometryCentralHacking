#pragma once

#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"

#include "geometrycentral/numerical/linear_solvers.h"

#include "polyscope/surface_mesh.h"

#include "homology_basis.h"
#include "utils.h"

#include <queue>

using namespace geometrycentral;
using namespace geometrycentral::surface;

//== Hodge Decomposition
class HodgeDecomposition {
  public:
    ManifoldSurfaceMesh& mesh;
    IntrinsicGeometryInterface& geom;

    const SparseMatrix<double>&d0, &d1, &star0, &star0inv, &star1, &star1inv,
        &star2, &star2inv;
    std::unique_ptr<PositiveDefiniteSolver<double>> L0FormSolver, L2FormSolver;

    HodgeDecomposition(ManifoldSurfaceMesh& mesh,
                       IntrinsicGeometryInterface& geom);
    ~HodgeDecomposition();

    //== Computes primal 0-form alpha, 1-form gamma, and 2-form beta such that
    // omega = d0 * alpha + star1 * d1.transpose() * star2inv * beta + gamma
    Vector<double>
    computePrimalExactPotential(const Vector<double>& omega) const;
    Vector<double>
    computePrimalCoexactPotential(const Vector<double>& omega) const;
    Vector<double>
    computePrimalHarmonicComponent(const Vector<double>& omega) const;
    // returned in order α, β, γ, where ω = dα + δβ + γ
    std::tuple<Vector<double>, Vector<double>, Vector<double>>
    computePrimalHodgeDecomposition(const Vector<double>& omega) const;

    //== Computes dual 0-form alpha, 1-form gamma, and 2-form beta such that
    // omega = d1.transpose() * alpha + star1 * d0 * star0inv * beta + gamma
    Vector<double> computeDualExactPotential(const Vector<double>& omega) const;
    Vector<double>
    computeDualCoexactPotential(const Vector<double>& omega) const;
    Vector<double>
    computeDualHarmonicComponent(const Vector<double>& omega) const;
    // returned in order α, β, γ, where ω = dα + δβ + γ
    std::tuple<Vector<double>, Vector<double>, Vector<double>>
    computeDualHodgeDecomposition(const Vector<double>& omega) const;

    //== Compute Hodge decompositions directly on MeshData
    // WARNING: only works properly on compressed meshes
    VertexData<double>
    computePrimalExactPotential(const EdgeData<double>& omega) const;
    FaceData<double>
    computePrimalCoexactPotential(const EdgeData<double>& omega) const;
    EdgeData<double>
    computePrimalHarmonicComponent(const EdgeData<double>& omega) const;
    // returned in order α, β, γ, where ω = dα + δβ + γ
    std::tuple<VertexData<double>, FaceData<double>, EdgeData<double>>
    computePrimalHodgeDecomposition(const EdgeData<double>& omega) const;

    FaceData<double>
    computeDualExactPotential(const EdgeData<double>& omega) const;
    VertexData<double>
    computeDualCoexactPotential(const EdgeData<double>& omega) const;
    EdgeData<double>
    computeDualHarmonicComponent(const EdgeData<double>& omega) const;
    // returned in order α, β, γ, where ω = dα + δβ + γ
    std::tuple<FaceData<double>, VertexData<double>, EdgeData<double>>
    computeDualHodgeDecomposition(const EdgeData<double>& omega) const;
};

//== Compute cohomology basis
std::vector<Vector<double>>
computePrimalCohomologyBasis(ManifoldSurfaceMesh& mesh,
                             IntrinsicGeometryInterface& geom);

std::vector<Vector<double>> computePrimalCohomologyBasis(
    ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
    const std::vector<std::vector<Halfedge>>& dualHomologyBasis);

std::vector<Vector<double>>
computePrimalCohomologyBasis(ManifoldSurfaceMesh& mesh,
                             IntrinsicGeometryInterface& geom,
                             const std::vector<std::vector<Halfedge>>& basis,
                             const HodgeDecomposition& hodgeDecomposition);

std::vector<Vector<double>>
computeDualCohomologyBasis(ManifoldSurfaceMesh& mesh,
                           IntrinsicGeometryInterface& geom);
std::vector<Vector<double>> computeDualCohomologyBasis(
    ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
    const std::vector<std::vector<Halfedge>>& primalHomologyBasis);

std::vector<Vector<double>> computeDualCohomologyBasis(
    ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
    const std::vector<std::vector<Halfedge>>& primalHomologyBasis,
    const HodgeDecomposition& hodgeDecomposition);

class TrivialConnection {
  public:
    ManifoldSurfaceMesh& mesh;
    IntrinsicGeometryInterface& geom;
    std::vector<std::vector<Halfedge>> primalHomologyBasis, dualHomologyBasis;
    std::vector<Vector<double>> primalCohomologyBasis, dualCohomologyBasis;
    HodgeDecomposition hodgeDecomposition;
    DenseMatrix<double> primalPeriodMatrix, dualPeriodMatrix;
    Eigen::ColPivHouseholderQR<DenseMatrix<double>> primalPeriodSolver,
        dualPeriodSolver;
    HalfedgeData<double> angleAlongHalfedge, angleAcrossHalfedge;
    HalfedgeData<Vector2> connectionTransportVectorAlongHalfedge,
        connectionTransportVectorAcrossHalfedge;
    size_t N = 1;

    TrivialConnection(ManifoldSurfaceMesh& mesh,
                      IntrinsicGeometryInterface& geom, size_t N = 1);
    void ensureHavePrimalData();
    void ensureHaveDualData();

    DenseMatrix<double> buildPrimalPeriodMatrix() const;
    DenseMatrix<double> buildDualPeriodMatrix() const;

    Vector<double>
    computePrimalCoexactComponent(const Vector<double>& singularity) const;
    Vector<double>
    computeDualCoexactComponent(const Vector<double>& singularity) const;

    Vector<double>
    computePrimalHarmonicComponent(const Vector<double>& delBeta) const;
    Vector<double>
    computeDualHarmonicComponent(const Vector<double>& delBeta) const;

    const HalfedgeData<Vector2>&
    computeVertexConnection(const Vector<double>& singularity);
    const HalfedgeData<Vector2>&
    computeFaceConnection(const Vector<double>& singularity);
    const HalfedgeData<Vector2>&
    computeVertexConnection(const FaceData<double>& singularity);
    const HalfedgeData<Vector2>&
    computeFaceConnection(const VertexData<double>& singularity);

    VertexData<Vector2> parallelTransport(Vertex v, Vector2 vec) const;
    FaceData<Vector2> parallelTransport(Face f, Vector2 vec) const;
};
