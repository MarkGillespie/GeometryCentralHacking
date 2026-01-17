#include "geometrycentral/surface/homology_generators.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/elementary_geometry.h"

#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geom;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
EdgeData<bool> psEdgeOrientations;

void vizGenerators(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                   const HomologyGenerators& generators) {
    //== primal
    for (size_t iP = 0; iP < generators.primalGenerators.size(); iP++) {
        std::vector<Vector3> pPrimal;
        std::vector<std::array<size_t, 2>> sPrimal;
        std::vector<size_t> iPrimal;

        size_t nG0 = pPrimal.size();
        for (Halfedge ij : generators.primalGenerators[iP]) {
            sPrimal.push_back({pPrimal.size(), pPrimal.size() + 1});
            pPrimal.push_back(geom.vertexPositions[ij.tailVertex()]);
            iPrimal.push_back(iP);
        }
        // push final vertex
        pPrimal.push_back(
            geom.vertexPositions
                [generators.primalGenerators[iP].back().tipVertex()]);
        iPrimal.push_back(iP);

        polyscope::registerCurveNetwork(
            "primal generator " + std::to_string(iP), pPrimal, sPrimal)
            ->addNodeScalarQuantity("index (primal)", iPrimal)
            ->setMapRange({0, generators.primalGenerators.size()})
            ->setEnabled(true);
    }

    //== dual
    auto fPos = [&](Halfedge ij) -> Vector3 {
        if (ij.face().isBoundaryLoop()) {
            return SurfacePoint(ij, .5).interpolate(geom.vertexPositions);
        } else {
            return SurfacePoint(ij.face(), Vector3{1, 1, 1} / 3.)
                .interpolate(geom.vertexPositions);
        }
    };
    for (size_t iP = 0; iP < generators.dualGenerators.size(); iP++) {
        std::vector<Vector3> pDual;
        std::vector<std::array<size_t, 2>> sDual;
        std::vector<size_t> iDual;

        size_t nG0 = pDual.size();
        for (Halfedge ij : generators.dualGenerators[iP]) {
            sDual.push_back({pDual.size() + 0, pDual.size() + 1});
            pDual.push_back(fPos(ij));
            iDual.push_back(iP);
        }
        // push final vertex
        pDual.push_back(fPos(generators.dualGenerators[iP].back().twin()));
        iDual.push_back(iP);

        polyscope::registerCurveNetwork("dual generator " + std::to_string(iP),
                                        pDual, sDual)
            ->addNodeScalarQuantity("index (dual)", iDual)
            ->setMapRange({0, generators.dualGenerators.size()})
            ->setEnabled(true);
    }
}

EdgeData<double>
circumcentricDualEdgeLengths(ManifoldSurfaceMesh& mesh,
                             IntrinsicGeometryInterface& geom) {
    geom.requireEdgeLengths();
    geom.requireEdgeCotanWeights();
    EdgeData<double> dualLengths(mesh);
    for (Edge e : mesh.edges())
        dualLengths[e] = geom.edgeCotanWeights[e] * geom.edgeLengths[e];
    geom.requireEdgeCotanWeights();
    geom.requireEdgeLengths();
    return dualLengths;
}

std::array<Vector2, 4> layoutDiamond(ManifoldSurfaceMesh& mesh,
                                     IntrinsicGeometryInterface& geom,
                                     Halfedge iHe) {

    // Conventions:
    //  - iHe points from vertex 2 to vertex 0, other vertices are numbered ccw
    //  - iHe is incident on face A, other is face B
    //  - halfedges within face are numbered CCW as A0, A1, A2 (etc),
    //    starting with iHe and twin(iHe)
    //  - When we lay out the triangle, p3 is at the origin and
    //    edge 3-0 is along the X-axis
    //  - flips is always ccw, so iHe points from vertex 3 --> 1 after

    // Gather index values
    Halfedge iHeA0 = iHe;
    Halfedge iHeA1 = iHeA0.next();
    Halfedge iHeA2 = iHeA1.next();
    Halfedge iHeB0 = iHe.twin();
    Halfedge iHeB1 = iHeB0.next();
    Halfedge iHeB2 = iHeB1.next();

    // Gather length values
    geom.requireEdgeLengths();
    const EdgeData<double>& l = geom.edgeLengths;
    double l01                = l[iHeA1.edge()];
    double l12                = l[iHeA2.edge()];
    double l23                = l[iHeB1.edge()];
    double l30                = l[iHeB2.edge()];
    double l02                = l[iHeA0.edge()];
    geom.unrequireEdgeLengths();

    // Lay out the vertices of the diamond
    Vector2 p3{0., 0.};
    Vector2 p0{l30, 0.};
    Vector2 p2 = layoutTriangleVertex(
        p3, p0, l02,
        l23); // involves more arithmetic than strictly necessary
    Vector2 p1 = layoutTriangleVertex(p2, p0, l01, l12);

    return {p0, p1, p2, p3};
}

EdgeData<double> barycentricDualEdgeLengths(ManifoldSurfaceMesh& mesh,
                                            IntrinsicGeometryInterface& geom) {
    EdgeData<double> dualLengths(mesh);

    for (Edge e : mesh.edges()) {
        // TODO: fix?
        if (e.isBoundary()) dualLengths[e] = 1;

        std::array<Vector2, 4> p = layoutDiamond(mesh, geom, e.halfedge());
        Vector2 pi = p[2], pj = p[0], pk = p[1], pl = p[3];
        Vector2 pijk   = (pi + pj + pk) / 3.;
        Vector2 pjil   = (pj + pi + pl) / 3.;
        dualLengths[e] = (pijk - pjil).norm();
    }

    return dualLengths;
}

FaceData<Vector3>
primalOneFormToFaceVectors(ManifoldSurfaceMesh& mesh,
                           VertexPositionGeometry& geom,
                           const EdgeData<double>& primalOneForm) {
    FaceData<Vector3> faceVectors(mesh);

    geom.requireFaceNormals();
    geom.requireFaceAreas();
    for (Face f : mesh.faces()) {
        Vector3 n      = geom.faceNormals[f];
        faceVectors[f] = Vector3::zero();
        for (Halfedge ij : f.adjacentHalfedges()) {
            Vector3 pi = geom.vertexPositions[ij.tailVertex()],
                    pj = geom.vertexPositions[ij.tipVertex()],
                    pk = geom.vertexPositions[ij.next().tipVertex()];

            double sign = (ij.orientation() ? 1.0 : -1.0);

            faceVectors[f] += sign * primalOneForm[ij.edge()] *
                              (cross(pk - pj, n) - cross(pi - pk, n));
        }

        faceVectors[f] = faceVectors[f] / (6. * geom.faceAreas[f]);
    }
    geom.unrequireFaceAreas();
    geom.unrequireFaceNormals();

    return faceVectors;
}

FaceData<Vector3>
dualOneFormToFaceVectors(ManifoldSurfaceMesh& mesh,
                         VertexPositionGeometry& geom,
                         const EdgeData<double>& dualOneForm) {
    FaceData<Vector3> faceVectors(mesh);
    EdgeData<double> primalOneForm(mesh,
                                   geom.hodge1Inverse * dualOneForm.raw());

    geom.requireFaceNormals();
    geom.requireFaceAreas();
    for (Face f : mesh.faces()) {
        Vector3 n      = geom.faceNormals[f];
        faceVectors[f] = Vector3::zero();
        for (Halfedge ij : f.adjacentHalfedges()) {
            Vector3 pi = geom.vertexPositions[ij.tailVertex()],
                    pj = geom.vertexPositions[ij.tipVertex()],
                    pk = geom.vertexPositions[ij.next().tipVertex()];

            double sign = (ij.orientation() ? 1.0 : -1.0);

            faceVectors[f] += sign * primalOneForm[ij.edge()] *
                              (cross(pk - pj, n) - cross(pi - pk, n));
        }

        faceVectors[f] = -cross(n, faceVectors[f]) / (6. * geom.faceAreas[f]);
        if (faceVectors[f].norm() > 1.) { // FIXME: remove this
            faceVectors[f] = unit(faceVectors[f]);
        }
    }
    geom.unrequireFaceAreas();
    geom.unrequireFaceNormals();

    return faceVectors;
}

// Returns harmonic 1-forms dual to the given generators
HarmonicGenerators
local_computeHarmonicGenerators(ManifoldSurfaceMesh& mesh,
                                IntrinsicGeometryInterface& geom,
                                const HomologyGenerators& generators) {
    auto sign = [&](Halfedge ij) -> double {
        return ij.orientation() ? 1. : -1.;
    };
    HarmonicGenerators result;
    geom.requireCotanLaplacian();
    geom.requireDECOperators();
    SparseMatrix<double> L0       = geom.cotanLaplacian;
    const SparseMatrix<double>&d0 = geom.d0, &d1 = geom.d1;
    const SparseMatrix<double>&hodge1 = geom.hodge1,
          hodge1Inv                   = geom.hodge1Inverse;
    SparseMatrix<double> L2           = d1 * hodge1Inv * d1.transpose();

    std::vector<Eigen::Triplet<double>> d1DirichletTriplets, d1NeumannTriplets;
    geom.requireFaceIndices();
    geom.requireEdgeIndices();
    const FaceData<size_t>& fIdx = geom.faceIndices;
    const EdgeData<size_t>& eIdx = geom.edgeIndices;
    for (Face f : mesh.faces()) {
        size_t iF = fIdx[f];
        for (Halfedge h : f.adjacentHalfedges()) {
            size_t iE = eIdx[h.edge()];
            if (h.edge().isBoundary()) {
                d1DirichletTriplets.emplace_back(iF, iE, 1);
                // d1NeumannTriplets.emplace_back(iF, iE, 0); // don't need to
                // set zero coefficient
            } else {
                d1DirichletTriplets.emplace_back(iF, iE, sign(h));
                d1NeumannTriplets.emplace_back(iF, iE, sign(h));
            }
        }
    }
    geom.unrequireEdgeIndices();
    geom.unrequireFaceIndices();

    size_t nE = mesh.nEdges(), nF = mesh.nFaces();
    SparseMatrix<double> d1Dirichlet(nF, nE), d1Neumann(nF, nE);
    d1Dirichlet.setFromTriplets(d1DirichletTriplets.begin(),
                                d1DirichletTriplets.end());
    d1Neumann.setFromTriplets(d1NeumannTriplets.begin(),
                              d1NeumannTriplets.end());

    SparseMatrix<double> L2Dirichlet =
        d1Dirichlet * hodge1Inv * d1Dirichlet.transpose();
    SparseMatrix<double> L2Neumann =
        d1Neumann * hodge1Inv * d1Neumann.transpose();

    VertexData<bool> isInteriorVertex(mesh, true);
    for (BoundaryLoop b : mesh.boundaryLoops()) {
        for (Vertex i : b.adjacentVertices()) isInteriorVertex[i] = false;
    }

    BlockDecompositionResult<double> decomp0 =
        blockDecomposeSquare(L0, isInteriorVertex.raw(), false);
    SparseMatrix<double> L0ii        = decomp0.AA;
    const SparseMatrix<double>& L0ib = decomp0.AB;

    //===== get primal harmonic generators by solving for jump across dual
    // generators
    result.primalGenerators.reserve(generators.dualGenerators.size());
    for (const std::vector<Halfedge>& dualGenerator :
         generators.dualGenerators) {
        EdgeData<double> jump(mesh, 0);
        for (Halfedge ij : dualGenerator) jump[ij.edge()] += sign(ij);

        // Solve for a jump-harmonic function w/ given jump and appropriate
        // boundary conditions
        Vector<double> alpha;
        switch (generators.primalType) {
        case HomologyType::Absolute: { // impose zero-Neumann boundary condition
                                       // on potential
            Vector<double> rhs = d0.transpose() * hodge1 * jump.raw();
            alpha              = solvePositiveDefinite(L0, rhs);
            break;
        }
        case HomologyType::Relative: { // impose zero-Dirichlet boundary
                                       // condition on potential
            Vector<double> rhs = d0.transpose() * hodge1 * jump.raw();
            alpha              = solvePositiveDefinite(L0, rhs);

            Vector<double> fullRHS = d0.transpose() * hodge1 * jump.raw(), iRHS,
                           bRHS;
            decomposeVector(decomp0, fullRHS, iRHS, bRHS);

            Vector<double> iPotential = solvePositiveDefinite(L0ii, iRHS);
            Vector<double> bPotential = Vector<double>::Zero(bRHS.size());
            alpha = reassembleVector(decomp0, iPotential, bPotential);
            break;
        }
        }

        EdgeData<double> gamma(mesh, d0 * alpha);
        for (Halfedge ij : dualGenerator) gamma[ij.edge()] -= sign(ij);
        result.primalGenerators.push_back(gamma);
    }

    //===== get dual harmonic generators by solving for jump across primal
    // generators
    result.dualGenerators.reserve(generators.primalGenerators.size());
    size_t iG = 0;
    for (const std::vector<Halfedge>& primalGenerator :
         generators.primalGenerators) {
        EdgeData<double> jump(mesh, 0);
        for (Halfedge ij : primalGenerator) jump[ij.edge()] += sign(ij);

        // Let jump = *^1 d0 α + d1^T β + γ,
        // where α has zero-Neumann boundary conditions, and β has
        // zero-Dirichlet. Then d0^T * d0 α = d0^T jump, and d1 *^{-1} d1^T β =
        // d1 *^{-1}  jump

        // Solve for a jump-harmonic function w/ given jump and appropriate
        // boundary conditions
        EdgeData<double> gamma;
        std::cout << "dual generator type: " << generators.dualType
                  << std::endl;
        switch (generators.dualType) {
        case HomologyType::Absolute: { // impose zero-Neumann boundary condition
                                       // on potential
            Vector<double> rhs  = d1Neumann * hodge1Inv * jump.raw();
            Vector<double> beta = solvePositiveDefinite(L2Neumann, rhs);
            polyscope::getSurfaceMesh("mesh")->addFaceScalarQuantity(
                "beta " + std::to_string(iG), beta);
            gamma = EdgeData<double>(mesh, d1Neumann.transpose() * beta);
            break;
        }
        case HomologyType::Relative: { // impose zero-Dirichlet boundary
                                       // condition on potential
            HERE();
            Vector<double> rhs  = d1Dirichlet * hodge1Inv * jump.raw();
            Vector<double> beta = solvePositiveDefinite(L2Dirichlet, rhs);
            polyscope::getSurfaceMesh("mesh")->addFaceScalarQuantity(
                "beta " + std::to_string(iG), beta);
            gamma = EdgeData<double>(mesh, d1Dirichlet.transpose() * beta);
            break;
        }
        }

        // { // alpha; impose zero-Neumann boundary condition on potential
        //   Vector<double> rhs = d0.transpose() * jump.raw();
        //   alpha = solvePositiveDefinite(L0, rhs);
        // }

        for (Halfedge ij : primalGenerator) gamma[ij.edge()] -= sign(ij);
        result.dualGenerators.push_back(gamma);

        // EdgeData<double> gamma(mesh, jump.raw() - hodge1 * d0 * alpha -
        // d1Dirichlet.transpose() * beta);
        // result.dualGenerators.push_back(-gamma);
        iG++;
    }

    // for (size_t iG = 0; iG < generators.dualGenerators.size(); iG++) {
    //   const std::vector<Halfedge>& dualGenerator =
    //   generators.dualGenerators[iG]; const EdgeData<double>& gamma =
    //   result.primalGenerators[iG]; for (size_t iH = 0; iH <
    //   generators.primalGenerators.size(); iH++) {
    //     const std::vector<Halfedge>& primalGenerator =
    //     generators.primalGenerators[iH]; double gammaIntegral = 0; for
    //     (Halfedge ij : primalGenerator) gammaIntegral += sign(ij) *
    //     gamma[ij.edge()]; std::cout << "int gamma (" << iG << ", " << iH <<
    //     "): " << gammaIntegral << std::endl;

    //     double intersectionCount = 0;
    //     for (Halfedge ij : primalGenerator) {
    //       for (Halfedge ab : dualGenerator) {
    //         if (ij.edge() == ab.edge()) intersectionCount -= sign(ij) *
    //         sign(ab);
    //       }
    //     }
    //     std::cout << "intersect (" << iG << ", " << iH << "): " <<
    //     intersectionCount << std::endl;
    //   }
    // }

    //===== dual generators
    result.dualGenerators.reserve(generators.dualGenerators.size());
    geom.unrequireDECOperators();
    geom.unrequireCotanLaplacian();
    return result;
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    HomologyGeneratorOptions opt = defaultHomologyGeneratorOptions;
    if (ImGui::Button("absolute primal generators")) {
        opt.generatorType = HomologyGeneratorType::AbsolutePrimal;
        vizGenerators(*mesh, *geom, computeHomologyGenerators(*mesh, opt));
    }
    if (ImGui::Button("relative dual generators")) {
        opt.generatorType = HomologyGeneratorType::RelativeDual;
        vizGenerators(*mesh, *geom, computeHomologyGenerators(*mesh, opt));
    }
    if (ImGui::Button("absolute primal relative dual generators")) {
        opt.generatorType = HomologyGeneratorType::AbsolutePrimalRelativeDual;
        HomologyGenerators homologyGenerators =
            computeHomologyGenerators(*mesh, opt);
        vizGenerators(*mesh, *geom, homologyGenerators);

        HarmonicGenerators harmonicGenerators =
            computeHarmonicGenerators(*mesh, *geom, homologyGenerators);
        psMesh->addOneFormIntrinsicVectorQuantity(
            "harmonic 0 (primal)", harmonicGenerators.primalGenerators[0],
            psEdgeOrientations);
        psMesh->addFaceVectorQuantity(
            "harmonic 0 (primal; vec)",
            primalOneFormToFaceVectors(*mesh, *geom,
                                       harmonicGenerators.primalGenerators[0]));
        psMesh->addFaceVectorQuantity(
            "harmonic 0 (dual; vec)",
            dualOneFormToFaceVectors(*mesh, *geom,
                                     harmonicGenerators.dualGenerators[0]));
    }
    if (ImGui::Button("relative primal generators")) {
        opt.generatorType = HomologyGeneratorType::RelativePrimal;
        vizGenerators(*mesh, *geom, computeHomologyGenerators(*mesh, opt));
    }
    if (ImGui::Button("absolute dual generators")) {
        opt.generatorType = HomologyGeneratorType::AbsoluteDual;
        vizGenerators(*mesh, *geom, computeHomologyGenerators(*mesh, opt));
    }
    if (ImGui::Button("relative primal absolute dual generators")) {
        opt.generatorType = HomologyGeneratorType::RelativePrimalAbsoluteDual;
        HomologyGenerators homologyGenerators =
            computeHomologyGenerators(*mesh, opt);
        vizGenerators(*mesh, *geom, homologyGenerators);

        HarmonicGenerators harmonicGenerators =
            computeHarmonicGenerators(*mesh, *geom, homologyGenerators);
        psMesh->addOneFormIntrinsicVectorQuantity(
            "harmonic 0 (primal)", harmonicGenerators.primalGenerators[0],
            psEdgeOrientations);
        psMesh->addFaceVectorQuantity(
            "harmonic 0 (primal; vec)",
            primalOneFormToFaceVectors(*mesh, *geom,
                                       harmonicGenerators.primalGenerators[0]));
        psMesh->addFaceVectorQuantity(
            "harmonic 0 (dual; vec)",
            dualOneFormToFaceVectors(*mesh, *geom,
                                     harmonicGenerators.dualGenerators[0]));
    }
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("Geometry program");
    args::Positional<std::string> inputFilename(parser, "mesh",
                                                "Mesh to be processed.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (const args::Help&) {
        std::cout << parser;
        return 0;
    } catch (const args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string filename = "../../meshes/bunny_small.obj";
    // Make sure a mesh name was given
    if (inputFilename) {
        filename = args::get(inputFilename);
    }

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = myCallback;

    // Load mesh
    std::tie(mesh, geom) = readManifoldSurfaceMesh(filename);
    std::cout << "Genus: " << mesh->genus() << std::endl;

    psEdgeOrientations = EdgeData<bool>(*mesh);
    for (Edge e : mesh->edges()) {
        psEdgeOrientations[e] =
            (e.firstVertex().getIndex() < e.secondVertex().getIndex());
    }

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh("mesh", geom->vertexPositions,
                                            mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
