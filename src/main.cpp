#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/remeshing.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

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

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    if (ImGui::Button("Remesh")) {
        RemeshOptions options;
        options.boundaryCondition = RemeshBoundaryCondition::Free;
        options.targetEdgeLength  = -2;
        options.smoothStyle       = RemeshSmoothStyle::Laplacian;
        remesh(*mesh, *geom, options);

        psMesh = polyscope::registerSurfaceMesh(
            "reremeshed", geom->vertexPositions, mesh->getFaceVertexList(),
            polyscopePermutations(*mesh));
    }
}

Vector3 vertexNormal(VertexPositionGeometry& geom, Vertex v) {
    Vector3 totalNormal = Vector3::zero();
    for (Corner c : v.adjacentCorners()) {
        Vector3 cornerNormal = geom.cornerAngle(c) * geom.faceNormal(c.face());
        totalNormal += cornerNormal;
    }
    return normalize(totalNormal);
}

Vector3 boundaryVertexTangent(VertexPositionGeometry& geom, Vertex v) {
    if (v.isBoundary()) {
        auto edgeVec = [&](Edge e) -> Vector3 {
            return (geom.vertexPositions[e.halfedge().tipVertex()] -
                    geom.vertexPositions[e.halfedge().tailVertex()])
                .normalize();
        };

        Vector3 totalTangent = Vector3::zero();
        for (Edge e : v.adjacentEdges()) {
            if (e.isBoundary()) {
                totalTangent += edgeVec(e);
            }
        }
        return totalTangent.normalize();
    } else {
        return Vector3::zero();
    }
}

Vector3 projectToPlane(Vector3 v, Vector3 norm) {
    return v - norm * dot(norm, v);
}
Vector3 projectToLine(Vector3 v, Vector3 tangent) {
    return tangent * dot(tangent, v);
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
    psMesh = polyscope::registerSurfaceMesh("original", geom->vertexPositions,
                                            mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

    RemeshOptions options;
    options.boundaryCondition = RemeshBoundaryCondition::Fixed;
    options.targetEdgeLength  = -2;
    options.smoothStyle       = RemeshSmoothStyle::Laplacian;
    options.maxIterations     = 1000;
    remesh(*mesh, *geom, options);

    psMesh = polyscope::registerSurfaceMesh(
        "remeshed (fixed)", geom->vertexPositions, mesh->getFaceVertexList(),
        polyscopePermutations(*mesh));

    std::tie(mesh, geom)      = readManifoldSurfaceMesh(filename);
    options.boundaryCondition = RemeshBoundaryCondition::Tangential;
    remesh(*mesh, *geom, options);
    psMesh = polyscope::registerSurfaceMesh(
        "remeshed (tangential)", geom->vertexPositions,
        mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    std::tie(mesh, geom)      = readManifoldSurfaceMesh(filename);
    options.boundaryCondition = RemeshBoundaryCondition::Free;
    remesh(*mesh, *geom, options);
    psMesh = polyscope::registerSurfaceMesh(
        "remeshed (free)", geom->vertexPositions, mesh->getFaceVertexList(),
        polyscopePermutations(*mesh));

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
