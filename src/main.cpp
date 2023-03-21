#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/utilities.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
// std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geom;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {}

Vector3 inR3(SurfacePoint& p) { return p.interpolate(geom->vertexPositions); }

std::array<Vector3, 2> getEdgeTangentBasis(Edge e) {
    geom->requireFaceNormals();
    Vector3 N = (geom->faceNormals[e.halfedge().face()] +
                 geom->faceNormals[e.halfedge().twin().face()])
                    .normalize();
    geom->unrequireFaceNormals();

    Vector3 T = (geom->vertexPositions[e.halfedge().tipVertex()] -
                 geom->vertexPositions[e.halfedge().tailVertex()])
                    .normalize();

    Vector3 B = cross(N, T);
    return {T, B};
}

Vector3 tangentInR3(const SurfacePoint& p, Vector2 v) {
    Vector3 result;
    switch (p.type) {
    case SurfacePointType::Vertex: {
        geom->requireVertexTangentBasis();
        auto basis = geom->vertexTangentBasis[p.vertex];
        result     = v[0] * basis[0] + v[1] * basis[1];
        geom->unrequireVertexTangentBasis();
        break;
    }
    case SurfacePointType::Edge: {
        auto basis = getEdgeTangentBasis(p.edge);
        result     = v[0] * basis[0] + v[1] * basis[1];
        break;
    }
    case SurfacePointType::Face: {
        geom->requireFaceTangentBasis();
        auto basis = geom->faceTangentBasis[p.face];
        result     = v[0] * basis[0] + v[1] * basis[1];
        geom->unrequireFaceTangentBasis();
        break;
    }
    }

    return result;
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
    // std::tie(mesh, geom) = readManifoldSurfaceMesh(filename);
    std::tie(mesh, geom) = readSurfaceMesh(filename);

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
        polyscope::guessNiceNameFromPath(filename), geom->vertexPositions,
        mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    GeodesicAlgorithmExact mmp(*mesh, *geom);
    SurfacePoint p(mesh->face(0), Vector3{1., 1., 1.} / 3.);
    SurfacePoint q(mesh->edge(57), 0.3);
    mmp.propagate({p});
    psMesh->addVertexScalarQuantity("dist", mmp.getDistanceFunction());
    // VertexData<Vector2> logMap(*mesh);
    // for (Vertex v : mesh->vertices()) {
    //     logMap[v] = mmp.getLog(v);
    // }

    // psMesh->addVertexParameterizationQuantity("log", logMap);

    // double maxErr = 0;
    // for (Vertex v : mesh->vertices()) {
    //     double trueDist = mmp.getDistance(v);
    //     double pathDist;
    //     mmp.traceBack(v, pathDist);
    //     maxErr = fmax(maxErr, abs(trueDist - pathDist));
    // }
    // WATCH(maxErr);

    std::vector<SurfacePoint> samplePts;
    // for (Vertex v : mesh->vertices()) samplePts.push_back(SurfacePoint(v));
    // for (Edge e : mesh->edges()) {
    //     samplePts.push_back(SurfacePoint(e, 0.33));
    //     samplePts.push_back(SurfacePoint(e, 0.66));
    // }
    // for (Face f : mesh->faces()) {
    //     for (size_t iP = 0; iP < 18; iP++) {
    //         samplePts.push_back(
    //             SurfacePoint(f, normalizeBarycentric(Vector3{
    //                                 unitRand(), unitRand(), unitRand()})));
    //     }
    // }
    // for (BoundaryLoop b : mesh->boundaryLoops()) {
    //     for (Halfedge he : b.adjacentHalfedges()) {
    //         if (he.face().isBoundaryLoop()) {
    //             he = he.twin();
    //         } else if (he.face().getIndex() == 34144) {
    //             std::cout << "wait, isn't this a boundary loop?" << vendl;
    //         }
    //         if (he.face().isBoundaryLoop()) {
    //             std::cout << "big problem" << vendl;
    //             continue;
    //         }

    //         // samplePts.push_back(SurfacePoint(he.tailVertex()));

    //         Face f = he.face();

    //         for (size_t iP = 0; iP < 0; iP++) {
    //             samplePts.push_back(
    //                 SurfacePoint(f, normalizeBarycentric(Vector3{
    //                                     unitRand(), unitRand(),
    //                                     unitRand()})));
    //         }
    //     }
    // }

    // std::vector<Vector3> samplePositions, sampleGradients;
    // std::vector<Vector2> sampleLogs;
    // std::vector<double> sampleDists;
    // for (SurfacePoint pt : samplePts) {
    //     sampleDists.push_back(mmp.getDistance(pt));
    //     samplePositions.push_back(inR3(pt));
    //     sampleLogs.push_back(mmp.getLog(pt));
    //     sampleGradients.push_back(tangentInR3(pt,
    //     mmp.getDistanceGradient(pt)));
    // }

    // auto psPoints = polyscope::registerPointCloud("Samples",
    // samplePositions); psPoints->addScalarQuantity("dist", sampleDists);
    // psPoints->addParameterizationQuantity("log", sampleLogs);
    // psPoints->addVectorQuantity("grad", sampleGradients);


    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
