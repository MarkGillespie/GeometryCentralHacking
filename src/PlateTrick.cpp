#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "ConePlacer.h"
#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using SO3        = Eigen::Matrix3d;
using Quaternion = Eigen::Vector4d;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

SO3 fromAxisAngle(double angle, Vector3 axis) {
    axis = axis.normalize();

    Eigen::Matrix3d K;
    // clang-format off
  K << 0,   -axis.z, axis.y,
    axis.z, 0    , -axis.x,
    -axis.y, axis.x,  0;
    // clang-format on

    return SO3::Identity() + sin(angle) * K + (1 - cos(angle)) * K * K;
}

double arg(Quaternion q) { return 2 * std::acos(q(0)); }
Vector3 imag(Quaternion q) { return Vector3{q(1), q(2), q(3)}; }
Vector3 axis(Quaternion q) { return imag(q) / sqrt(1 - q(0) * q(0)); }
// Vector3 fromEigen(Eigen::Vector3d v) { return Vector3{v(0), v(1), v(2)}; }

SO3 fromQuaternion(const Quaternion& q) {
    if (abs(q(0) - 1) < 1e-8) {
        return SO3::Identity();
    } else {
        return fromAxisAngle(arg(q), axis(q));
    }
}

polyscope::CurveNetwork* vizPath(std::vector<Quaternion> path,
                                 std::string name = "curve") {
    std::vector<Vector3> positions;
    std::vector<Vector3> tangents;

    for (Quaternion q : path) {
        SO3 R       = fromQuaternion(q);
        Vector3 pos = fromEigen(R.col(0));
        Vector3 dir = fromEigen(R.col(1));
        positions.push_back(pos);
        tangents.push_back(dir);
    }

    auto curve = polyscope::registerCurveNetworkLine(name, positions);
    curve->addNodeVectorQuantity("dir", tangents)->setEnabled(true);
    return curve;
}

void doTrick() {
    int res   = 500;
    double dt = 1. / (res - 1.);

    auto getPath = [&](double a) -> std::vector<Quaternion> {
        std::vector<Quaternion> path;
        for (size_t iS = 0; iS < (size_t)res; iS++) {
            double t      = iS * dt * 2 * M_PI;
            double offset = sqrt(fabs(2 * a * (1 - a) * (1 - cos(t))));
            // path.push_back(
            //     Quaternion{(1 - a) * cos(t) + a, 0, (1 - a) * sin(t),
            //     offset});
            path.push_back(
                Quaternion{(1 - a) * cos(t) + a, 0, offset, (1 - a) * sin(t)});
        }
        return path;
    };

    vizPath(getPath(0), "path 0");
    vizPath(getPath(0.1), "path 0.1");
    vizPath(getPath(0.25), "path 0.25");
    vizPath(getPath(0.5), "path 0.5");
    vizPath(getPath(0.9), "path 0.9");
    vizPath(getPath(0.95), "path 0.95");
}

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

    doTrick();

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
