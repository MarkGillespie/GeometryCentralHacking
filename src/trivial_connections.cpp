#include "trivial_connections.h"

HodgeDecomposition::HodgeDecomposition(ManifoldSurfaceMesh& mesh_,
                                       IntrinsicGeometryInterface& geom_)
    : mesh(mesh_), geom(geom_), d0(geom.d0), d1(geom.d1), star0(geom.hodge0),
      star0inv(geom.hodge0Inverse), star1(geom.hodge1),
      star1inv(geom.hodge1Inverse), star2(geom.hodge2),
      star2inv(geom.hodge2Inverse) {

    geom.requireDECOperators();

    SparseMatrix<double> L0Form = d0.transpose() * star1 * d0;
    SparseMatrix<double> L2Form = d1 * star1inv * d1.transpose();

    shiftDiagonal<double>(L0Form, 1e-8);
    shiftDiagonal<double>(L2Form, 1e-8);

    L0FormSolver.reset(new PositiveDefiniteSolver<double>(L0Form));
    L2FormSolver.reset(new PositiveDefiniteSolver<double>(L2Form));
}

HodgeDecomposition::~HodgeDecomposition() { geom.unrequireDECOperators(); }

Vector<double> HodgeDecomposition::computePrimalExactPotential(
    const Vector<double>& omega) const {
    Vector<double> rhs   = d0.transpose() * star1 * omega;
    Vector<double> alpha = L0FormSolver->solve(rhs);
    return alpha;
}

Vector<double> HodgeDecomposition::computePrimalCoexactPotential(
    const Vector<double>& omega) const {
    Vector<double> starBeta = L2FormSolver->solve(d1 * omega);

    return star2inv * starBeta;
}

std::tuple<Vector<double>, Vector<double>, Vector<double>>
HodgeDecomposition::computePrimalHodgeDecomposition(
    const Vector<double>& omega) const {
    Vector<double> alpha   = computePrimalExactPotential(omega);
    Vector<double> beta    = computePrimalCoexactPotential(omega);
    Vector<double> dAlpha  = d0 * alpha;
    Vector<double> delBeta = star1inv * d1.transpose() * star2 * beta;
    Vector<double> gamma   = omega - dAlpha - delBeta;
    return std::make_tuple(alpha, beta, gamma);
}

Vector<double> HodgeDecomposition::computePrimalHarmonicComponent(
    const Vector<double>& omega) const {
    Vector<double> alpha   = computePrimalExactPotential(omega);
    Vector<double> beta    = computePrimalCoexactPotential(omega);
    Vector<double> dAlpha  = d0 * alpha;
    Vector<double> delBeta = star1inv * d1.transpose() * star2 * beta;
    return omega - dAlpha - delBeta;
}

Vector<double> HodgeDecomposition::computeDualExactPotential(
    const Vector<double>& omega) const {
    Vector<double> rhs   = d1 * star1inv * omega;
    Vector<double> alpha = L2FormSolver->solve(rhs);
    return alpha;
}

Vector<double> HodgeDecomposition::computeDualCoexactPotential(
    const Vector<double>& omega) const {
    Vector<double> rhs      = d0.transpose() * omega;
    Vector<double> starBeta = L0FormSolver->solve(rhs);

    return star0 * starBeta;
}

std::tuple<Vector<double>, Vector<double>, Vector<double>>
HodgeDecomposition::computeDualHodgeDecomposition(
    const Vector<double>& omega) const {
    Vector<double> alpha   = computeDualExactPotential(omega);
    Vector<double> beta    = computeDualCoexactPotential(omega);
    Vector<double> dAlpha  = d1.transpose() * alpha;
    Vector<double> delBeta = star1 * d0 * star0inv * beta;
    Vector<double> gamma   = omega - dAlpha - delBeta;
    return std::make_tuple(alpha, beta, gamma);
}

Vector<double> HodgeDecomposition::computeDualHarmonicComponent(
    const Vector<double>& omega) const {
    Vector<double> alpha   = computeDualExactPotential(omega);
    Vector<double> beta    = computeDualCoexactPotential(omega);
    Vector<double> dAlpha  = d1.transpose() * alpha;
    Vector<double> delBeta = star1 * d0 * star0inv * beta;
    return omega - dAlpha - delBeta;
}

VertexData<double> HodgeDecomposition::computePrimalExactPotential(
    const EdgeData<double>& omega) const {
    return VertexData<double>(mesh, computePrimalExactPotential(omega.raw()));
}

FaceData<double> HodgeDecomposition::computePrimalCoexactPotential(
    const EdgeData<double>& omega) const {
    return FaceData<double>(mesh, computePrimalCoexactPotential(omega.raw()));
}

EdgeData<double> HodgeDecomposition::computePrimalHarmonicComponent(
    const EdgeData<double>& omega) const {
    return EdgeData<double>(mesh, computePrimalHarmonicComponent(omega.raw()));
}

std::tuple<VertexData<double>, FaceData<double>, EdgeData<double>>
HodgeDecomposition::computePrimalHodgeDecomposition(
    const EdgeData<double>& omega) const {
    std::tuple<Vector<double>, Vector<double>, Vector<double>> decomp =
        computePrimalHodgeDecomposition(omega.raw());
    return std::make_tuple(VertexData<double>(mesh, std::get<0>(decomp)),
                           FaceData<double>(mesh, std::get<1>(decomp)),
                           EdgeData<double>(mesh, std::get<2>(decomp)));
}

FaceData<double> HodgeDecomposition::computeDualExactPotential(
    const EdgeData<double>& omega) const {
    return FaceData<double>(mesh, computeDualExactPotential(omega.raw()));
}
VertexData<double> HodgeDecomposition::computeDualCoexactPotential(
    const EdgeData<double>& omega) const {
    return VertexData<double>(mesh, computeDualCoexactPotential(omega.raw()));
}
EdgeData<double> HodgeDecomposition::computeDualHarmonicComponent(
    const EdgeData<double>& omega) const {
    return EdgeData<double>(mesh, computeDualHarmonicComponent(omega.raw()));
}
std::tuple<FaceData<double>, VertexData<double>, EdgeData<double>>
HodgeDecomposition::computeDualHodgeDecomposition(
    const EdgeData<double>& omega) const {
    std::tuple<Vector<double>, Vector<double>, Vector<double>> decomp =
        computeDualHodgeDecomposition(omega.raw());
    return std::make_tuple(FaceData<double>(mesh, std::get<0>(decomp)),
                           VertexData<double>(mesh, std::get<0>(decomp)),
                           EdgeData<double>(mesh, std::get<0>(decomp)));
}

std::vector<Vector<double>>
computePrimalCohomologyBasis(ManifoldSurfaceMesh& mesh,
                             IntrinsicGeometryInterface& geom) {
    return computePrimalCohomologyBasis(mesh, geom,
                                        computeDualHomologyBasis(mesh));
}

std::vector<Vector<double>> computePrimalCohomologyBasis(
    ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
    const std::vector<std::vector<Halfedge>>& dualHomologyBasis) {
    HodgeDecomposition hodgeDecomposition(mesh, geom);
    return computePrimalCohomologyBasis(mesh, geom, dualHomologyBasis,
                                        hodgeDecomposition);
}

std::vector<Vector<double>> computePrimalCohomologyBasis(
    ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
    const std::vector<std::vector<Halfedge>>& dualHomologyBasis,
    const HodgeDecomposition& hodgeDecomposition) {

    geom.requireEdgeIndices();
    const EdgeData<size_t>& eIdx = geom.edgeIndices;

    std::vector<Vector<double>> cohomologyBasis;
    for (const std::vector<Halfedge>& generator : dualHomologyBasis) {
        // build closed primal 1-form
        Vector<double> omega = Vector<double>::Zero(mesh.nEdges());
        for (Halfedge he : generator) {
            omega(eIdx[he.edge()]) += he.orientation() ? 1 : -1;
        }

        cohomologyBasis.push_back(
            hodgeDecomposition.computePrimalHarmonicComponent(omega));
    }
    geom.unrequireEdgeIndices();

    return cohomologyBasis;
}

std::vector<Vector<double>>
computeDualCohomologyBasis(ManifoldSurfaceMesh& mesh,
                           IntrinsicGeometryInterface& geom) {
    return computeDualCohomologyBasis(mesh, geom,
                                      computePrimalHomologyBasis(mesh));
}
std::vector<Vector<double>> computeDualCohomologyBasis(
    ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
    const std::vector<std::vector<Halfedge>>& primalHomologyBasis) {
    HodgeDecomposition hodgeDecomposition(mesh, geom);
    return computeDualCohomologyBasis(mesh, geom, primalHomologyBasis,
                                      hodgeDecomposition);
}

std::vector<Vector<double>> computeDualCohomologyBasis(
    ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom,
    const std::vector<std::vector<Halfedge>>& primalHomologyBasis,
    const HodgeDecomposition& hodgeDecomposition) {

    geom.requireEdgeIndices();
    const EdgeData<size_t>& eIdx = geom.edgeIndices;

    std::vector<Vector<double>> cohomologyBasis;
    for (const std::vector<Halfedge>& generator : primalHomologyBasis) {
        // build closed dual 1-form
        Vector<double> omega = Vector<double>::Zero(mesh.nEdges());
        for (Halfedge he : generator) {
            omega(eIdx[he.edge()]) += he.orientation() ? 1 : -1;
        }

        cohomologyBasis.push_back(
            hodgeDecomposition.computeDualHarmonicComponent(omega));
    }
    geom.unrequireEdgeIndices();

    return cohomologyBasis;
}

TrivialConnection::TrivialConnection(ManifoldSurfaceMesh& mesh_,
                                     IntrinsicGeometryInterface& geom_,
                                     size_t N_)
    : mesh(mesh_), geom(geom_), hodgeDecomposition(mesh, geom), N(N_) {}

void TrivialConnection::ensureHavePrimalData() {
    if (!primalHomologyBasis.empty()) return;

    primalHomologyBasis = computePrimalHomologyBasis(mesh);

    // If we're a topological sphere, return early
    if (primalHomologyBasis.empty()) return;

    // We build the dual cohomology basis and manually convert to the basis of
    // primal 1-forms (computing the dual basis is easier since we need the
    // primal homology generators later)
    primalCohomologyBasis = computeDualCohomologyBasis(
        mesh, geom, primalHomologyBasis, hodgeDecomposition);
    const SparseMatrix<double>& star1inv = geom.hodge1Inverse;
    for (Vector<double>& harmonicForm : primalCohomologyBasis) {
        harmonicForm = star1inv * harmonicForm;
    }

    primalPeriodMatrix = buildPrimalPeriodMatrix();

    primalPeriodSolver =
        Eigen::ColPivHouseholderQR<DenseMatrix<double>>(primalPeriodMatrix);
}

void TrivialConnection::ensureHaveDualData() {
    if (!dualHomologyBasis.empty()) return;

    dualHomologyBasis = computeDualHomologyBasis(mesh);

    // If we're a topological sphere, return early
    if (dualHomologyBasis.empty()) return;

    // We build the primal cohomology basis and manually convert to the basis of
    // dual 1-forms (computing the primal basis is easier since we need the dual
    // homology generators later)
    dualCohomologyBasis = computePrimalCohomologyBasis(
        mesh, geom, dualHomologyBasis, hodgeDecomposition);
    const SparseMatrix<double>& star1 = geom.hodge1;
    for (Vector<double>& harmonicForm : dualCohomologyBasis) {
        harmonicForm = star1 * harmonicForm;
    }

    dualPeriodMatrix = buildDualPeriodMatrix();

    dualPeriodSolver =
        Eigen::ColPivHouseholderQR<DenseMatrix<double>>(dualPeriodMatrix);
}

DenseMatrix<double> TrivialConnection::buildPrimalPeriodMatrix() const {
    geom.requireEdgeIndices();
    const EdgeData<size_t>& eIdx = geom.edgeIndices;

    DenseMatrix<double> periodMatrix(primalHomologyBasis.size(),
                                     primalCohomologyBasis.size());

    for (size_t i = 0; i < primalHomologyBasis.size(); i++) {
        const std::vector<Halfedge>& l_i = primalHomologyBasis[i];
        for (size_t j = 0; j < primalCohomologyBasis.size(); j++) {
            const Vector<double>& xi_j = primalCohomologyBasis[j];

            double integral = 0;
            for (Halfedge he : l_i) {
                double sign = he.orientation() ? 1 : -1;
                integral += sign * xi_j(eIdx[he.edge()]);
            }
            periodMatrix(i, j) = integral;
        }
    }
    geom.unrequireEdgeIndices();

    return periodMatrix;
}

DenseMatrix<double> TrivialConnection::buildDualPeriodMatrix() const {
    geom.requireEdgeIndices();
    const EdgeData<size_t>& eIdx = geom.edgeIndices;

    DenseMatrix<double> periodMatrix(dualHomologyBasis.size(),
                                     dualCohomologyBasis.size());

    for (size_t i = 0; i < dualHomologyBasis.size(); i++) {
        const std::vector<Halfedge>& l_i = dualHomologyBasis[i];
        for (size_t j = 0; j < dualCohomologyBasis.size(); j++) {
            const Vector<double>& xi_j = dualCohomologyBasis[j];

            double integral = 0;
            for (Halfedge he : l_i) {
                double sign = he.orientation() ? 1 : -1;
                integral += sign * xi_j(eIdx[he.edge()]);
            }
            periodMatrix(i, j) = integral;
        }
    }
    geom.unrequireEdgeIndices();

    return periodMatrix;
}

Vector<double> TrivialConnection::computePrimalCoexactComponent(
    const Vector<double>& singularity) const {

    FaceData<double> faceGaussianCurvatures(mesh, -M_PI);
    geom.requireCornerScaledAngles();
    for (Corner c : mesh.corners()) {
        faceGaussianCurvatures[c.face()] += geom.cornerScaledAngles[c];
    }
    Vector<double> rhs(mesh.nFaces());
    geom.requireFaceIndices();
    const FaceData<size_t>& fIdx = geom.faceIndices;
    for (Face f : mesh.faces()) {
        rhs(fIdx[f]) = -faceGaussianCurvatures[f] +
                       2 * M_PI * singularity(fIdx[f]) / ((double)N);
    }
    geom.unrequireFaceIndices();

    Vector<double> starBeta = hodgeDecomposition.L2FormSolver->solve(rhs);

    const SparseMatrix<double>& star1inv = geom.hodge1Inverse;
    const SparseMatrix<double>& d1       = geom.d1;

    return star1inv * d1.transpose() * starBeta;
}

Vector<double> TrivialConnection::computeDualCoexactComponent(
    const Vector<double>& singularity) const {
    Vector<double> rhs(mesh.nVertices());
    geom.requireVertexIndices();
    geom.requireVertexGaussianCurvatures();
    const VertexData<size_t>& vIdx = geom.vertexIndices;
    const VertexData<double>& K    = geom.vertexGaussianCurvatures;
    for (Vertex v : mesh.vertices()) {
        rhs(vIdx[v]) = -K[v] + 2 * M_PI * singularity(vIdx[v]) / ((double)N);
    }
    geom.unrequireVertexGaussianCurvatures();
    geom.unrequireVertexIndices();

    Vector<double> starBeta = hodgeDecomposition.L0FormSolver->solve(rhs);

    const SparseMatrix<double>& star1 = geom.hodge1;
    const SparseMatrix<double>& d0    = geom.d0;

    return star1 * d0 * starBeta;
}

Vector<double> TrivialConnection::computePrimalHarmonicComponent(
    const Vector<double>& delBeta) const {

    // If we're a topological sphere, harmonic component is zero
    if (primalHomologyBasis.empty())
        return Vector<double>::Zero(delBeta.rows());

    geom.requireTransportVectorsAlongHalfedge();
    geom.requireEdgeIndices();

    const HalfedgeData<Vector2>& leviCivita =
        geom.transportVectorsAlongHalfedge;
    const EdgeData<size_t>& eIdx = geom.edgeIndices;

    Vector<double> signature(primalHomologyBasis.size());
    for (size_t i = 0; i < primalHomologyBasis.size(); i++) {
        const std::vector<Halfedge>& li = primalHomologyBasis[i];

        // represent holonomy as a complex number to reduce mod 2π
        Vector2 holonomy = Vector2{1, 0};
        for (Halfedge he : li) {
            double sign          = he.orientation() ? 1 : -1;
            double signedDelBeta = sign * delBeta(eIdx[he.edge()]);
            holonomy *= leviCivita[he] * Vector2::fromAngle(signedDelBeta);
        }
        holonomy     = holonomy.pow(N);
        signature(i) = -holonomy.arg() / ((double)N);
    }
    geom.unrequireEdgeIndices();
    geom.unrequireTransportVectorsAlongHalfedge();

    Vector<double> z = primalPeriodSolver.solve(signature);

    Vector<double> gamma = Vector<double>::Zero(mesh.nEdges());
    for (size_t i = 0; i < primalCohomologyBasis.size(); i++) {
        gamma += z(i) * primalCohomologyBasis[i];
    }

    return gamma;
}

Vector<double> TrivialConnection::computeDualHarmonicComponent(
    const Vector<double>& delBeta) const {

    // If we're a topological sphere, harmonic component is zero
    if (dualHomologyBasis.empty()) return Vector<double>::Zero(delBeta.rows());

    geom.requireTransportVectorsAcrossHalfedge();
    geom.requireEdgeIndices();

    const HalfedgeData<Vector2>& leviCivita =
        geom.transportVectorsAcrossHalfedge;
    const EdgeData<size_t>& eIdx = geom.edgeIndices;

    Vector<double> signature(dualHomologyBasis.size());
    for (size_t i = 0; i < dualHomologyBasis.size(); i++) {
        const std::vector<Halfedge>& li = dualHomologyBasis[i];

        // represent holonomy as a complex number to reduce mod 2π
        Vector2 holonomy = Vector2{1, 0};
        for (Halfedge he : li) {
            double sign          = he.orientation() ? 1 : -1;
            double signedDelBeta = sign * delBeta(eIdx[he.edge()]);
            holonomy *= leviCivita[he] * Vector2::fromAngle(signedDelBeta);
        }
        holonomy     = holonomy.pow(N);
        signature(i) = -holonomy.arg() / ((double)N);
    }
    geom.unrequireEdgeIndices();
    geom.unrequireTransportVectorsAcrossHalfedge();

    Vector<double> z     = dualPeriodSolver.solve(signature);
    Vector<double> gamma = Vector<double>::Zero(mesh.nEdges());
    for (size_t i = 0; i < dualCohomologyBasis.size(); i++) {
        gamma += z(i) * dualCohomologyBasis[i];
    }

    return gamma;
}

const HalfedgeData<Vector2>&
TrivialConnection::computeVertexConnection(const Vector<double>& singularity) {
    ensureHavePrimalData();

    bool satisfyGaussBonnet =
        abs(singularity.sum() - ((int)N) * mesh.eulerCharacteristic()) < 1e-8;
    if (!satisfyGaussBonnet) {
        std::cerr << "Singularities sum to " << singularity.sum()
                  << " whereas N times the Euler characteristic is "
                  << ((int)N) * mesh.eulerCharacteristic() << " (N=" << N
                  << ", Euler characteristic=" << mesh.eulerCharacteristic()
                  << " )" << std::endl;
        throw_verbose_runtime_error(
            "invalid singularity input to trivial connections");
        return connectionTransportVectorAlongHalfedge;
    }

    Vector<double> delBeta = computePrimalCoexactComponent(singularity);
    Vector<double> gamma   = computePrimalHarmonicComponent(delBeta);
    Vector<double> phi     = delBeta + gamma;

    angleAlongHalfedge                     = HalfedgeData<double>(mesh);
    connectionTransportVectorAlongHalfedge = HalfedgeData<Vector2>(mesh);
    geom.requireEdgeIndices();
    geom.requireTransportVectorsAlongHalfedge();

    const EdgeData<size_t>& eIdx = geom.edgeIndices;
    const HalfedgeData<Vector2>& leviCivita =
        geom.transportVectorsAlongHalfedge;

    for (Halfedge he : mesh.halfedges()) {
        double sign      = he.orientation() ? 1 : -1;
        double signedPhi = sign * phi(eIdx[he.edge()]);

        angleAlongHalfedge[he] = leviCivita[he].arg() + signedPhi;
        connectionTransportVectorAlongHalfedge[he] =
            leviCivita[he] * Vector2::fromAngle(signedPhi);
    }
    geom.unrequireTransportVectorsAlongHalfedge();
    geom.unrequireEdgeIndices();

    return connectionTransportVectorAlongHalfedge;
}

const HalfedgeData<Vector2>&
TrivialConnection::computeFaceConnection(const Vector<double>& singularity) {
    ensureHaveDualData();

    bool satisfyGaussBonnet =
        abs(singularity.sum() - ((int)N) * mesh.eulerCharacteristic()) < 1e-8;
    if (!satisfyGaussBonnet) {
        std::cerr << "Singularities sum to " << singularity.sum()
                  << " whereas N times the Euler characteristic is "
                  << ((int)N) * mesh.eulerCharacteristic() << " (N=" << N
                  << ", Euler characteristic=" << mesh.eulerCharacteristic()
                  << " )" << std::endl;
        throw_verbose_runtime_error(
            "invalid singularity input to trivial connections");
        return connectionTransportVectorAcrossHalfedge;
    }

    Vector<double> delBeta = computeDualCoexactComponent(singularity);
    Vector<double> gamma   = computeDualHarmonicComponent(delBeta);
    Vector<double> phi     = delBeta + gamma;

    connectionTransportVectorAcrossHalfedge = HalfedgeData<Vector2>(mesh);
    angleAcrossHalfedge                     = HalfedgeData<double>(mesh);
    geom.requireEdgeIndices();
    geom.requireTransportVectorsAcrossHalfedge();

    const EdgeData<size_t>& eIdx = geom.edgeIndices;
    const HalfedgeData<Vector2>& leviCivita =
        geom.transportVectorsAcrossHalfedge;

    for (Halfedge he : mesh.halfedges()) {
        if (he.edge().isBoundary()) {
            angleAcrossHalfedge[he] = std::numeric_limits<double>::quiet_NaN();
            connectionTransportVectorAcrossHalfedge[he] = Vector2::undefined();
        } else {
            double sign      = he.orientation() ? 1 : -1;
            double signedPhi = sign * phi(eIdx[he.edge()]);

            angleAcrossHalfedge[he] = leviCivita[he].arg() + signedPhi;
            connectionTransportVectorAcrossHalfedge[he] =
                leviCivita[he] * Vector2::fromAngle(signedPhi);
        }
    }
    geom.unrequireTransportVectorsAcrossHalfedge();
    geom.unrequireEdgeIndices();

    return connectionTransportVectorAcrossHalfedge;
}

const HalfedgeData<Vector2>& TrivialConnection::computeVertexConnection(
    const FaceData<double>& singularity) {
    return computeVertexConnection(singularity.raw());
}

const HalfedgeData<Vector2>& TrivialConnection::computeFaceConnection(
    const VertexData<double>& singularity) {
    return computeFaceConnection(singularity.raw());
}

VertexData<Vector2>
TrivialConnection::parallelTransport(Vertex vStart, Vector2 vecStart) const {
    VertexData<Vector2> field(mesh, Vector2::undefined());

    std::queue<Vertex> toVisit;

    field[vStart] = vecStart;
    toVisit.push(vStart);

    while (!toVisit.empty()) {
        Vertex v = toVisit.front();
        toVisit.pop();

        verbose_assert(!std::isnan(field[v].x), "??? Nan");

        for (Halfedge he : v.outgoingHalfedges()) {
            Vertex w = he.tipVertex();
            if (std::isnan(field[w].x)) {
                field[w] =
                    connectionTransportVectorAlongHalfedge[he] * field[v];
                toVisit.push(w);
            }
        }
    }

    return field;
}

FaceData<Vector2> TrivialConnection::parallelTransport(Face fStart,
                                                       Vector2 vecStart) const {
    FaceData<Vector2> field(mesh, Vector2::undefined());

    std::queue<Face> toVisit;

    field[fStart] = vecStart;
    toVisit.push(fStart);

    while (!toVisit.empty()) {
        Face f = toVisit.front();
        toVisit.pop();

        if (f.isBoundaryLoop()) continue;

        verbose_assert(!std::isnan(field[f].x), "??? Nan");

        for (Halfedge he : f.adjacentHalfedges()) {
            Face g = he.twin().face();
            if (std::isnan(field[g].x)) {
                field[g] =
                    connectionTransportVectorAcrossHalfedge[he] * field[f];
                toVisit.push(g);
            }
        }
    }

    return field;
}
