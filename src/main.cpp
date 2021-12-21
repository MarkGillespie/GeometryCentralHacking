#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "ConePlacer.h"
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
    psMesh = polyscope::registerSurfaceMesh("mesh", geometry->vertexPositions,
                                            mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

    ConePlacer fang(*mesh, *geometry);
    std::vector<std::pair<Vertex, double>> cones = fang.optimalCones();

    std::vector<Vector3> conePositions;
    for (const std::pair<Vertex, double>& cone : cones)
        conePositions.push_back(geometry->vertexPositions[cone.first]);

    polyscope::registerPointCloud("cones", conePositions);

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
