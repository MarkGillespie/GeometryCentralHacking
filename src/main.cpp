#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/remeshing.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "homology_basis.h"
#include "trivial_connections.h"
#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geom;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    // if (ImGui::Button("Remesh")) {
    // }
}

Vector3 center(Vertex v) { return geom->vertexPositions[v]; }
Vector3 center(Edge e) {
    SurfacePoint p(e, 0.5);
    return p.interpolate(geom->vertexPositions);
}
Vector3 center(Face f) {
    SurfacePoint p(f, Vector3{1, 1, 1} / 3.);
    return p.interpolate(geom->vertexPositions);
}
Vector3 faceCenter(Halfedge he) {
    // if he.face() is an ordinary face, return its center. If he.face() is a
    // boundary loop, then just use he's center
    return he.isInterior() ? center(he.face()) : center(he.edge());
}
Vector3 dir(Edge e) {
    return (center(e.halfedge().vertex()) -
            center(e.halfedge().twin().vertex()))
        .normalize();
}
Vector3 dualDir(Edge e) {
    return (center(e.halfedge().face()) - center(e.halfedge().twin().face()))
        .normalize();
}

Vector3 inR3(const SurfacePoint& p) {
    // SurfacePoint pIn = opt->intTri->equivalentPointOnInput(p);
    // return opt->intTri->equivalentPointOnInput(p).interpolate(
    //     geom->vertexPositions);
    return p.interpolate(geom->vertexPositions);
}

void plotPrimalTree(const VertexData<Halfedge>& primalTree) {
    std::vector<Vector3> nodePositions;
    std::vector<std::array<size_t, 2>> edgeIndices;
    std::vector<Vector3> edgeOrientations;

    for (Vertex v : mesh->vertices()) {
        if (primalTree[v] == Halfedge()) continue;
        Halfedge he = primalTree[v];
        edgeIndices.push_back({nodePositions.size(), nodePositions.size() + 1});
        nodePositions.push_back(center(he.tailVertex()));
        nodePositions.push_back(center(he.tipVertex()));
        double sign = he.orientation() ? 1 : -1;
        edgeOrientations.push_back(sign * dir(he.edge()));
    }

    polyscope::registerCurveNetwork("primal tree", nodePositions, edgeIndices)
        ->addEdgeVectorQuantity("orientation", edgeOrientations);
}

void plotDualTree(const FaceAndBoundaryData<Halfedge>& dualTree) {
    std::vector<Vector3> nodePositions;
    std::vector<std::array<size_t, 2>> edgeIndices;
    std::vector<Vector3> edgeOrientations;

    auto drawHalfedge = [&](Halfedge he) {
        if (he == Halfedge()) return;
        edgeIndices.push_back({nodePositions.size(), nodePositions.size() + 1});
        nodePositions.push_back(faceCenter(he));
        nodePositions.push_back(faceCenter(he.twin()));
        double sign = he.orientation() ? 1 : -1;
        edgeOrientations.push_back(sign * dualDir(he.edge()));
    };

    for (Face f : mesh->faces()) drawHalfedge(dualTree[f]);
    for (BoundaryLoop b : mesh->boundaryLoops()) {
        drawHalfedge(dualTree[b]);
        if (dualTree[b] != Halfedge()) {
            WATCH2(b, dualTree[b]);
        }
    }

    polyscope::registerCurveNetwork("dual tree", nodePositions, edgeIndices)
        ->addEdgeVectorQuantity("orientation", edgeOrientations);
}

void plotCut(const std::vector<Edge>& cut) {
    std::vector<Vector3> nodePositions;
    std::vector<std::array<size_t, 2>> edgeIndices;
    std::vector<Vector3> edgeOrientations;

    for (Edge e : cut) {
        edgeIndices.push_back({nodePositions.size(), nodePositions.size() + 1});
        nodePositions.push_back(center(e.halfedge().tailVertex()));
        nodePositions.push_back(center(e.halfedge().tipVertex()));
    }

    polyscope::registerCurveNetwork("cut", nodePositions, edgeIndices);
}

void plotPrimalCurves(const std::vector<std::vector<Halfedge>>& curves) {
    std::vector<Vector3> nodePositions;
    std::vector<std::array<size_t, 2>> edgeIndices;
    std::vector<Vector3> edgeOrientations;

    for (const std::vector<Halfedge>& curve : curves) {
        for (Halfedge he : curve) {
            edgeIndices.push_back(
                {nodePositions.size(), nodePositions.size() + 1});
            nodePositions.push_back(center(he.tailVertex()));
            nodePositions.push_back(center(he.tipVertex()));
            double sign = he.orientation() ? 1 : -1;
            edgeOrientations.push_back(sign * dir(he.edge()));
        }
    }

    polyscope::registerCurveNetwork("primal loops", nodePositions, edgeIndices)
        ->addEdgeVectorQuantity("orientation", edgeOrientations);
}

void plotDualCurves(const std::vector<std::vector<Halfedge>>& curves) {
    std::vector<Vector3> nodePositions;
    std::vector<std::array<size_t, 2>> edgeIndices;
    std::vector<Vector3> edgeOrientations;

    for (const std::vector<Halfedge>& curve : curves) {
        for (Halfedge he : curve) {
            edgeIndices.push_back(
                {nodePositions.size(), nodePositions.size() + 1});
            nodePositions.push_back(center(he.face()));
            nodePositions.push_back(center(he.twin().face()));
            double sign = he.orientation() ? 1 : -1;
            edgeOrientations.push_back(sign * dualDir(he.edge()));
        }
    }

    polyscope::registerCurveNetwork("dual loops", nodePositions, edgeIndices)
        ->addEdgeVectorQuantity("orientation", edgeOrientations);
}

void findParallelField(bool face = false) {
    if (face) {
        TrivialConnection tc(*mesh, *geom);
        VertexData<double> singularity(*mesh, 0);
        singularity[0]   = 1;
        singularity[100] = 1;
        singularity *= 0.5 * mesh->eulerCharacteristic();
        tc.computeFaceConnection(singularity);
        FaceData<Vector2> parallelField =
            tc.parallelTransport(mesh->face(0), Vector2{1, 0});
        psMesh->addFaceIntrinsicVectorQuantity("face field", parallelField);
    } else {
        TrivialConnection tc(*mesh, *geom);
        FaceData<double> singularity(*mesh, 0);
        singularity[0]    = 1;
        singularity[1546] = 1;
        singularity *= 0.5 * mesh->eulerCharacteristic();
        tc.computeVertexConnection(singularity);
        VertexData<Vector2> parallelField =
            tc.parallelTransport(mesh->vertex(0), Vector2{1, 0});
        psMesh->addVertexIntrinsicVectorQuantity("vertex field", parallelField);
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
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string filename = "../../meshes/spot.obj";
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

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh("mesh", geom->vertexPositions,
                                            mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

    auto primalTree = buildPrimalSpanningTree(*mesh);
    plotPrimalTree(primalTree);
    auto dualTree = buildDualSpanningCotree(*mesh, primalTree);
    plotDualTree(dualTree);

    // auto primalCurves = computePrimalHomologyBasis(*mesh);
    // plotPrimalCurves(primalCurves);

    auto bothCurves = computePrimalAndDualHomologyBases(*mesh);
    plotPrimalCurves(std::get<0>(bothCurves));
    plotDualCurves(std::get<1>(bothCurves));

    // auto dualTree = buildDualSpanningCotree(*mesh);
    // plotDualTree(dualTree);
    // auto primalTree = buildPrimalSpanningTree(*mesh, dualTree);
    // plotPrimalTree(primalTree);
    // std::vector<std::vector<Halfedge>>
    // dualLoops =
    //     computeDualHomologyBasis(*mesh);
    // plotDualCurves(dualLoops);

    // geom->requireEdgeLengths();
    // auto primalLoops = computePrimalHomologyBasis(*mesh, geom->edgeLengths);
    // plotPrimalCurves(primalLoops);
    // geom->unrequireEdgeLengths();


    // geom->requireEdgeLengths();
    // auto bothCurves =
    //     computePrimalAndDualHomologyBases(*mesh, geom->edgeLengths);
    // plotPrimalCurves(std::get<0>(bothCurves));
    // plotDualCurves(std::get<1>(bothCurves));
    // geom->unrequireEdgeLengths();

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
