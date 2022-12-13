#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"

#include "polyscope/surface_mesh.h"

#include "utils.h"

#include <queue>

using namespace geometrycentral;
using namespace geometrycentral::surface;

template <typename T>
class FaceAndBoundaryData {
  public:
    FaceAndBoundaryData() {}
    FaceAndBoundaryData(ManifoldSurfaceMesh& mesh) {
        faceData     = FaceData<T>(mesh);
        boundaryData = BoundaryLoopData<T>(mesh);
    }
    FaceAndBoundaryData(ManifoldSurfaceMesh& mesh, T defaultValue) {
        faceData     = FaceData<T>(mesh, defaultValue);
        boundaryData = BoundaryLoopData<T>(mesh, defaultValue);
    }

    // Access with an element pointer
    T& operator[](BoundaryLoop b) { return boundaryData[b]; }
    const T& operator[](BoundaryLoop b) const { return boundaryData[b]; }
    T& operator[](Face f) {
        return f.isBoundaryLoop() ? boundaryData[f.asBoundaryLoop()]
                                  : faceData[f];
    }
    const T& operator[](Face f) const {
        return f.isBoundaryLoop() ? boundaryData[f.asBoundaryLoop()]
                                  : faceData[f];
    }

    size_t size() const { return faceData.size() + boundaryData.size(); }

  protected:
    FaceData<T> faceData;
    BoundaryLoopData<T> boundaryData;
};

//==== Tree-Cotree
std::vector<std::vector<Halfedge>>
computePrimalHomologyBasis(ManifoldSurfaceMesh& mesh,
                           bool connectBoundary = false);

std::vector<std::vector<Halfedge>>
computePrimalHomologyBasis(ManifoldSurfaceMesh& mesh,
                           const EdgeData<double>& edgeCost,
                           bool connectBoundary = false);

std::vector<std::vector<Halfedge>>
computeDualHomologyBasis(ManifoldSurfaceMesh& mesh,
                         bool connectBoundary = false);

std::vector<std::vector<Halfedge>>
computeDualHomologyBasis(ManifoldSurfaceMesh& mesh,
                         const EdgeData<double>& edgeCost,
                         bool connectBoundary = false);

void computeHomologyBases(
    ManifoldSurfaceMesh& mesh, const VertexData<Halfedge>& primalTree,
    const FaceAndBoundaryData<Halfedge>& dualTree,
    std::vector<std::vector<Halfedge>>* primalBasis = nullptr,
    std::vector<std::vector<Halfedge>>* dualBasis   = nullptr,
    bool ignoreBoundaryEdges                        = false);

// Encodes the primal spanning tree via a VertexData<Halfedge> primalParent
// where primalParent[v].vertex() is v's parent, and primalParent[v].edge()
// is the edge from the parent to v
// If dualTree is specified, the primal tree avoids those edges
// If connectBoundary is true, then the spanning tree is allowed to jump between
// vertices on the boundary of the mesh. This is useful for computing [relative
// homology generators](https://en.wikipedia.org/wiki/Relative_homology)
// (essentially homology generators which are allowed to start and end on the
// boundary of the mesh)
VertexData<Halfedge> buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                                             bool connectBoundary = false,
                                             Vertex root          = Vertex());
VertexData<Halfedge> buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                                             const EdgeData<double>& edgeWeight,
                                             bool connectBoundary = false,
                                             Vertex root          = Vertex());
VertexData<Halfedge>
buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                        const FaceAndBoundaryData<Halfedge>& dualTree,
                        bool connectBoundary = true, Vertex root = Vertex());
VertexData<Halfedge>
buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight,
                        const FaceAndBoundaryData<Halfedge>& dualTree,
                        bool connectBoundary = true, Vertex root = Vertex());
bool inPrimalSpanningTree(Halfedge he, const VertexData<Halfedge>& primalTree);

// Similarly, encodes the dual spanning tree via a FaceData<Halfedge> dualParent
// where dualParent[f].face() is f's parent and dualParent[f].edge() is the
// (primal) edge between f and its parent
// If primalTree is specified, the dual tree avoids those edges
// If connectBoundary is true, then the spanning tree is allowed to jump between
// faces on the boundary of the mesh. This is useful for computing [relative
// homology generators](https://en.wikipedia.org/wiki/Relative_homology)
// (essentially homology generators which are allowed to start and end on the
// boundary of the mesh)
FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh, bool connectBoundary = false,
                        Face root = Face());
FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight,
                        bool connectBoundary = false, Face root = Face());
FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh,
                        const VertexData<Halfedge>& primalTree,
                        bool connectBoundary = true, Face root = Face());
FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight,
                        const VertexData<Halfedge>& primalTree,
                        bool connectBoundary = true, Face root = Face());
bool inDualSpanningCotree(Halfedge he,
                          const FaceAndBoundaryData<Halfedge>& dualTree);

std::vector<Halfedge> primalTreeLoop(Halfedge he,
                                     const VertexData<Halfedge>& primalTree,
                                     bool trimLoop = true);
std::vector<Halfedge>
dualTreeLoop(Halfedge he, const FaceAndBoundaryData<Halfedge>& dualTree,
             bool trimLoop = true);

// Walk up dual tree from he1 and he2 and return path between them
std::vector<Halfedge>
dualTreePath(Halfedge heSrc, Halfedge heDst,
             const FaceAndBoundaryData<Halfedge>& dualTree,
             bool trimLoop = true);
