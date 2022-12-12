#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"

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
computePrimalHomologyBasis(ManifoldSurfaceMesh& mesh);
std::vector<std::vector<Halfedge>>
computePrimalHomologyBasis(ManifoldSurfaceMesh& mesh,
                           const EdgeData<double>& edgeCost);
std::vector<std::vector<Halfedge>>
computeDualHomologyBasis(ManifoldSurfaceMesh& mesh);
std::tuple<std::vector<std::vector<Halfedge>>,
           std::vector<std::vector<Halfedge>>>
computePrimalAndDualHomologyBases(ManifoldSurfaceMesh& mesh,
                                  const EdgeData<double>& edgeCost);
std::tuple<std::vector<std::vector<Halfedge>>,
           std::vector<std::vector<Halfedge>>>
computePrimalAndDualHomologyBases(ManifoldSurfaceMesh& mesh);

// Encodes the primal spanning tree via a VertexData<Halfedge> primalParent
// where primalParent[v].vertex() is v's parent, and primalParent[v].edge()
// is the edge from the parent to v
// If dualTree is specified, the primal tree avoids those edges
VertexData<Halfedge> buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh);
VertexData<Halfedge>
buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight);
VertexData<Halfedge>
buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                        const FaceAndBoundaryData<Halfedge>& dualTree);
VertexData<Halfedge>
buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight,
                        const FaceAndBoundaryData<Halfedge>& dualTree);
bool inPrimalSpanningTree(Halfedge he, const VertexData<Halfedge>& primalTree);

// Similarly, encodes the dual spanning tree via a FaceData<Halfedge> dualParent
// where dualParent[f].face() is f's parent and dualParent[f].edge() is the
// (primal) edge between f and its parent
// If primalTree is specified, the dual tree avoids those edges
FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh);
FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight);
FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh,
                        const VertexData<Halfedge>& primalTree);
FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight,
                        const VertexData<Halfedge>& primalTree);
bool inDualSpanningCotree(Halfedge he,
                          const FaceAndBoundaryData<Halfedge>& dualTree);

std::vector<Halfedge> primalTreeLoop(Halfedge he,
                                     const VertexData<Halfedge>& primalTree);
std::vector<Halfedge>
dualTreeLoop(Halfedge he, const FaceAndBoundaryData<Halfedge>& dualTree);

// Walk up dual tree from he1 and he2 and return path between them
std::vector<Halfedge>
dualTreePath(Halfedge heSrc, Halfedge heDst,
             const FaceAndBoundaryData<Halfedge>& dualTree);
