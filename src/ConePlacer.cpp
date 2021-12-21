#include "ConePlacer.h"

ConePlacer::ConePlacer(ManifoldSurfaceMesh& mesh_,
                       IntrinsicGeometryInterface& geom_)
    : mesh(mesh_), geom(geom_) {

    vIdx = mesh.getVertexIndices();

    isInterior = Vector<bool>(mesh.nVertices());
    for (Vertex v : mesh.vertices()) {
        isInterior(vIdx[v]) = !v.isBoundary();
    }
    nVertices = mesh.nVertices();
    nInterior = mesh.nInteriorVertices();
    nBoundary = nVertices - nInterior;

    geom.requireCotanLaplacian();
    SparseMatrix<double> L = geom.cotanLaplacian;
    shiftDiagonal(L, 1e-12);
    geom.unrequireCotanLaplacian();

    Ldecomp = blockDecomposeSquare(L, isInterior);

    Lii = Ldecomp.AA;
    Lib = Ldecomp.AB;


    geom.requireVertexDualAreas();
    double surfaceArea = 0;
    for (Vertex v : mesh.vertices()) surfaceArea += geom.vertexDualAreas[v];
    std::vector<Eigen::Triplet<double>> TA, TAinv;
    for (Vertex v : mesh.vertices()) {
        size_t iV = vIdx[v];
        TA.emplace_back(iV, iV, sqrt(geom.vertexDualAreas[v] / surfaceArea));
        TAinv.emplace_back(iV, iV, sqrt(surfaceArea / geom.vertexDualAreas[v]));
    }
    rA = SparseMatrix<double>(nVertices, nVertices);
    rA.setFromTriplets(std::begin(TA), std::end(TA));
    rA_inv = SparseMatrix<double>(nVertices, nVertices);
    rA_inv.setFromTriplets(std::begin(TA), std::end(TA));

    BlockDecompositionResult<double> rAdecomp =
        blockDecomposeSquare(rA, isInterior);
    rA_ii = rAdecomp.AA;

    BlockDecompositionResult<double> rA_inv_decomp =
        blockDecomposeSquare(rA_inv, isInterior);
    rA_inv_ii = rA_inv_decomp.AA;
    geom.unrequireVertexDualAreas();

    geom.requireVertexAngleSums();
    k_orig = Vector<double>(nVertices);
    for (Vertex v : mesh.vertices())
        k_orig[v] =
            (v.isBoundary() ? M_PI : 2 * M_PI) - geom.vertexAngleSums[v];
    geom.unrequireVertexAngleSums();
    decomposeVector(Ldecomp, k_orig, k_orig_i, k_orig_b);
}

std::vector<std::pair<Vertex, double>> ConePlacer::optimalCones() {
    size_t iter = 0;

    weights                   = Vector<double>::Ones(nInteriorVertices);
    double curvatureChange    = 1;
    Vector<double> curvatures = k_orig;
    while (curvatureChange > 1e-8 && iter++ < 25) {
        Vector<double> newCurvatures = DRminimize(curvatures);

        curvatureChange = (newCurvatures - curvatures).norm();
        curvatures      = newCurvatures;
    }

    Vector<double> boundaryZeros = Vector<double>::Zero(nBoundaryVertices);
    Vector<double> allCurvatures =
        reassembleVector(Ldecomp, curvatures, boundaryZeros);

    std::vector<std::pair<Vertex, double>> cones;
    for (Vertex v : mesh.vertices()) {
        if (v.isBoundary()) continue;
        double k = allCurvatures[vIdx[v]];
        if (abs(k) > 1e-10) {
            cones.push_back(std::make_pair(v, k));
        }
    }
    return cones;
}

Vector<double> ConePlacer::DRminimize(const VertexData<double>& initialS) {
    Vector<double> u, v, s = initialS;

    double residual = 1;
    size_t iter     = 0;
    while (iter++ < 50) {
        u        = stepU(s);
        v        = stepV(s, u);
        residual = (u - v).norm();
        if (residual < epsilon) {
            break;
        } else {
            s = stepV(s, u, v);
        }
    }

    return rA_inv * u;
}

Vector<double> ConePlacer::stepU(const Vector<double>& s) {}

Vector<double> ConePlacer::stepV(const Vector<double>& s,
                                 const Vector<double>& u) {}

Vector<double> ConePlacer::stepS(const Vector<double>& s,
                                 const Vector<double>& u,
                                 const Vector<double>& v) {}
