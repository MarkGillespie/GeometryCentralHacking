#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geom;

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
    std::tie(mesh, geom) = readSurfaceMesh(filename);

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
        polyscope::guessNiceNameFromPath(filename), geom->vertexPositions,
        mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    std::cout << "Orig mesh: nV = " << mesh->nVertices()
              << ", nE = " << mesh->nEdges() << ", nF = " << mesh->nFaces()
              << ", nH = " << mesh->nHalfedges() << std::endl;

    std::vector<std::vector<size_t>> polygons;
    std::vector<std::vector<std::tuple<size_t, size_t>>> twins;

    HalfedgeData<size_t> iFaceHe(*mesh);
    for (Face f : mesh->faces()) {
        Halfedge h     = f.halfedge();
        Halfedge hCurr = h;
        size_t iH      = 0;
        do {
            iFaceHe[hCurr] = iH;
            iH++;
            hCurr = hCurr.next();
        } while (hCurr != h);
    }
    geom->requireVertexIndices();
    geom->requireFaceIndices();
    for (Face f : mesh->faces()) {
        Halfedge h     = f.halfedge();
        Halfedge hCurr = h;
        size_t iH      = 0;
        polygons.push_back(std::vector<size_t>{});
        twins.push_back(std::vector<std::tuple<size_t, size_t>>{});
        do {
            polygons.back().push_back(geom->vertexIndices[hCurr.tailVertex()]);
            size_t iTwinFace = hCurr.twin() == hCurr
                                   ? INVALID_IND
                                   : geom->faceIndices[hCurr.twin().face()];
            size_t iTwinFaceHe =
                hCurr.twin() == hCurr ? INVALID_IND : iFaceHe[hCurr.twin()];
            twins.back().push_back(std::make_tuple(iTwinFace, iTwinFaceHe));
            hCurr = hCurr.next();
        } while (hCurr != h);
    }
    geom->unrequireFaceIndices();
    geom->unrequireVertexIndices();

    SurfaceMesh recon(polygons, twins);
    std::cout << "New mesh: nV = " << recon.nVertices()
              << ", nE = " << recon.nEdges() << ", nF = " << recon.nFaces()
              << ", nH = " << recon.nHalfedges() << std::endl;

    polyscope::registerSurfaceMesh("recon", geom->vertexPositions,
                                   recon.getFaceVertexList(),
                                   polyscopePermutations(recon));

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
