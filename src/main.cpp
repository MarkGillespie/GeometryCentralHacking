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
            "harmonic 0", harmonicGenerators.primalGenerators[0],
            psEdgeOrientations);
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
        vizGenerators(*mesh, *geom, computeHomologyGenerators(*mesh, opt));
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
    psMesh = polyscope::registerSurfaceMesh(
        polyscope::guessNiceNameFromPath(filename), geom->vertexPositions,
        mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
