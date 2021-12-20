#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

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
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {}

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
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(filename);
    std::cout << "Genus: " << mesh->genus() << std::endl;

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
        polyscope::guessNiceNameFromPath(filename), geometry->vertexPositions,
        mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    IntegerCoordinatesIntrinsicTriangulation icit(*mesh, *geometry);

    auto drawIntrinsicEdges = [&](std::string name) {
        EdgeData<std::vector<SurfacePoint>> tracedIntrinsicEdges =
            icit.traceAllIntrinsicEdgesAlongInput();
        std::vector<Vector3> positions;
        std::vector<std::array<size_t, 2>> edgeEdges;
        for (Edge e : icit.intrinsicMesh->edges()) {
            for (size_t iP = 0; iP + 1 < tracedIntrinsicEdges[e].size(); iP++)
                edgeEdges.push_back(
                    {positions.size() + iP, positions.size() + iP + 1});
            for (const SurfacePoint& q : tracedIntrinsicEdges[e])
                positions.push_back(q.interpolate(geometry->vertexPositions));
        }
        return polyscope::registerCurveNetwork(name, positions, edgeEdges);
    };

    icit.flipToDelaunay();
    drawIntrinsicEdges("idt");

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
