
EdgeData<double> primalWeights, dualWeights; // if not set, build tree via BFS


// If opt contains initialized primalWeights, builds an MST. Otherwise, builds a
// spanning tree via BFS


VertexData<Halfedge> buildPrimalSpanningTree(
    ManifoldSurfaceMesh& mesh, const EdgeData<double>& primalWeights,
    const FaceData<Halfedge>* dualTree = nullptr,
    HomologyGeneratorOptions opt       = defaultHomologyGeneratorOptions);

FaceData<Halfedge> buildDualSpanningTree(
    ManifoldSurfaceMesh& mesh, const EdgeData<double>& dualWeights,
    const VertexData<Halfedge>* primalTree = nullptr,
    HomologyGeneratorOptions opt           = defaultHomologyGeneratorOptions);


#include <queue> // priority queue

if (opt.primalWeights.size() > 0)
    return buildPrimalSpanningTree(mesh, opt.primalWeights, dualTree, opt);
VertexData<Halfedge> buildPrimalSpanningTree(
    ManifoldSurfaceMesh& mesh, const EdgeData<double>& primalWeights,
    const FaceData<Halfedge>* dualTree, HomologyGeneratorOptions opt) {
    HomologyType homologyType = primalHomologyType(opt.generatorType);

    VertexData<Halfedge> primalTree(mesh, Halfedge());
    VertexData<bool> visited(mesh, false);

    auto inDualTreePtr = [&](Halfedge ij) -> bool {
        return dualTree && inDualTree(ij, *dualTree);
    };

    typedef std::pair<double, Vertex>
        DistanceVertex; // vertex w/ distance, gets sorted by distance
    std::priority_queue<DistanceVertex, std::vector<DistanceVertex>,
                        std::greater<DistanceVertex>>
        toVisit; // build min heap

    if (homologyType == HomologyType::Relative && mesh.hasBoundary()) {
        // if we're looking for relative generators, we should connect together
        // all boundary vertices. The easiest way to do that is just to start by
        // pushing all boundary vertices onto the queue
        for (BoundaryLoop b : mesh.boundaryLoops()) {
            for (Vertex v : b.adjacentVertices()) {
                toVisit.push(std::make_pair(0, v));
                visited[v] = true;
            }
        }
    } else {
        // if we want absolute generators, we should just pick an arbitrary root
        // vertex and start there
        Vertex root =
            (opt.primalRoot == Vertex()) ? mesh.vertex(0) : opt.primalRoot;
        toVisit.push(std::make_pair(0, root));
        visited[root] = true;
    }

    while (!toVisit.empty()) {
        Vertex i  = std::get<1>(toVisit.top());
        double wi = std::get<0>(toVisit.top());
        toVisit.pop();
        for (Halfedge ji : i.incomingHalfedges()) {
            Vertex j = ji.tailVertex();
            if (!inDualTreePtr(ji) && !visited[j]) {
                primalTree[j] = ji;
                toVisit.push(std::make_pair(wi + primalWeights[ji.edge()], j));
                visited[j] = true;
            }
        }
    }
    return primalTree;
}


if (opt.dualWeights.size() > 0)
    return buildDualSpanningTree(mesh, opt.dualWeights, primalTree, opt);
// we take the convention that tree[vertex].twin().vertex() is the parent
// vertex, and tree[face].twin().face() is the parent face
FaceData<Halfedge> buildDualSpanningTree(ManifoldSurfaceMesh& mesh,
                                         const EdgeData<double>& dualWeights,
                                         const VertexData<Halfedge>* primalTree,
                                         HomologyGeneratorOptions opt) {
    HomologyType homologyType = dualHomologyType(opt.generatorType);

    FaceData<Halfedge> dualTree(mesh, Halfedge());
    FaceData<bool> visited(mesh, false);

    auto inPrimalTreePtr = [&](Halfedge ij) -> bool {
        return primalTree && inPrimalTree(ij, *primalTree);
    };

    typedef std::pair<double, Face>
        DistanceFace; // face w/ distance, gets sorted by distance
    std::priority_queue<DistanceFace, std::vector<DistanceFace>,
                        std::greater<DistanceFace>>
        toVisit; // build min heap

    if (homologyType == HomologyType::Relative && mesh.hasBoundary()) {
        // if we're looking for relative generators, we should connect together
        // all boundary faces. The easiest way to do that is just to start by
        // pushing all boundary faces onto the queue
        for (BoundaryLoop b : mesh.boundaryLoops()) {
            for (Halfedge ij : b.adjacentHalfedges()) {
                if (inPrimalTreePtr(ij)) continue;
                toVisit.push(std::make_pair(0, ij.twin().face()));
                visited[ij.twin().face()] = true;
                dualTree[ij.twin().face()] =
                    ij.twin(); // mark boundary edges as used in dual tree
            }
        }
    } else {
        // if we want absolute generators, we should just pick an arbitrary root
        // vertex and start there
        Face root = (opt.dualRoot == Face()) ? mesh.face(0) : opt.dualRoot;
        toVisit.push(std::make_pair(0, root));
        visited[root] = true;
    }

    while (!toVisit.empty()) {
        Face i    = std::get<1>(toVisit.top());
        double wi = std::get<0>(toVisit.top());
        toVisit.pop();
        for (Halfedge ij : i.adjacentHalfedges()) {
            Face j = ij.twin().face();
            if (!j.isBoundaryLoop() && !inPrimalTreePtr(ij) && !visited[j]) {
                dualTree[j] = ij.twin();
                toVisit.push(std::make_pair(wi + dualWeights[ij.edge()], j));
                visited[j] = true;
            }
        }
    }

    return dualTree;
}
