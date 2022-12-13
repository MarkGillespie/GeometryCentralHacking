#include "homology_basis.h"

std::vector<std::vector<Halfedge>>
computePrimalHomologyBasis(ManifoldSurfaceMesh& mesh, bool connectBoundary) {
    EdgeData<double> cost(mesh, 1);
    return computePrimalHomologyBasis(mesh, cost, connectBoundary);
}

std::vector<std::vector<Halfedge>>
computePrimalHomologyBasis(ManifoldSurfaceMesh& mesh,
                           const EdgeData<double>& edgeCost,
                           bool connectBoundary) {

    VertexData<Halfedge> primalTree =
        buildPrimalSpanningTree(mesh, edgeCost, connectBoundary);
    FaceAndBoundaryData<Halfedge> dualTree =
        buildDualSpanningCotree(mesh, edgeCost, primalTree, !connectBoundary);

    std::vector<std::vector<Halfedge>> primalBasis;
    computeHomologyBases(mesh, primalTree, dualTree, &primalBasis, nullptr,
                         connectBoundary);

    return primalBasis;
}

std::vector<std::vector<Halfedge>>
computeDualHomologyBasis(ManifoldSurfaceMesh& mesh, bool connectBoundary) {
    EdgeData<double> cost(mesh, 1);
    return computeDualHomologyBasis(mesh, cost, connectBoundary);
}

std::vector<std::vector<Halfedge>>
computeDualHomologyBasis(ManifoldSurfaceMesh& mesh,
                         const EdgeData<double>& edgeCost,
                         bool connectBoundary) {
    FaceAndBoundaryData<Halfedge> dualTree =
        buildDualSpanningCotree(mesh, edgeCost, connectBoundary);
    VertexData<Halfedge> primalTree =
        buildPrimalSpanningTree(mesh, edgeCost, dualTree, !connectBoundary);

    std::vector<std::vector<Halfedge>> dualBasis;
    computeHomologyBases(mesh, primalTree, dualTree, nullptr, &dualBasis,
                         !connectBoundary);

    return dualBasis;
}

void computeHomologyBases(ManifoldSurfaceMesh& mesh,
                          const VertexData<Halfedge>& primalTree,
                          const FaceAndBoundaryData<Halfedge>& dualTree,
                          std::vector<std::vector<Halfedge>>* primalBasis,
                          std::vector<std::vector<Halfedge>>* dualBasis,
                          bool ignoreBoundaryEdges) {

    for (Edge e : mesh.edges()) {
        if (!inPrimalSpanningTree(e.halfedge(), primalTree) &&
            !inDualSpanningCotree(e.halfedge(), dualTree)) {
            if (ignoreBoundaryEdges && e.isBoundary()) continue;

            if (primalBasis) {
                primalBasis->push_back(
                    primalTreeLoop(e.halfedge(), primalTree));
            }

            if (dualBasis) {
                dualBasis->push_back(dualTreeLoop(e.halfedge(), dualTree));
            }
        }
    }
}

VertexData<Halfedge> buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                                             bool connectBoundary,
                                             Vertex root) {
    FaceAndBoundaryData<Halfedge> emptyDualTree;
    return buildPrimalSpanningTree(mesh, emptyDualTree, connectBoundary, root);
}

VertexData<Halfedge> buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                                             const EdgeData<double>& edgeWeight,
                                             bool connectBoundary,
                                             Vertex root) {

    VertexData<Halfedge> primalTree(mesh, Halfedge());

    // Distance estimates to build a shortest-paths tree
    VertexData<double> dist(mesh, std::numeric_limits<double>::infinity());

    VertexData<bool> visited(mesh, false);

    // make sure we start on a boundary vertex if connectBoundary is true
    // this is necessary for the loop tracing code later
    if (connectBoundary && (root == Vertex() || !root.isBoundary()) &&
        mesh.nBoundaryLoops() > 0) {
        root = mesh.boundaryLoop(0).halfedge().vertex();
    }
    if (root == Vertex()) root = mesh.vertex(0);
    std::priority_queue<std::tuple<double, Vertex>> toVisit;
    toVisit.push({0, root});
    dist[root] = 0;

    while (!toVisit.empty()) {
        Vertex v = std::get<1>(toVisit.top());
        toVisit.pop();

        // if we've already processed this vertex, skip it
        if (visited[v]) continue;

        visited[v] = true;

        if (connectBoundary && v.isBoundary()) {
            // If connectBoundary is true, propagate to the whole boundary once
            // we've hit one vertex
            for (BoundaryLoop b : mesh.boundaryLoops()) {
                for (Vertex vb : b.adjacentVertices()) {
                    visited[vb] = true;
                    dist[vb]    = dist[v];
                    for (Halfedge he : vb.outgoingHalfedges()) {
                        Vertex w = he.twin().vertex();
                        if (w.isBoundary() || visited[w]) continue;

                        double vDist = dist[v] + edgeWeight[he.edge()];
                        if (vDist < dist[w]) {
                            dist[w]       = vDist;
                            primalTree[w] = he;

                            // std::priority_queue is a max-heap, so we sort by
                            // negative cost to pop the closest vertex next
                            toVisit.push({-dist[w], w});
                        }
                    }
                }
            }
        } else {
            for (Halfedge he : v.outgoingHalfedges()) {
                Vertex w = he.twin().vertex();
                if (visited[w]) continue;

                double vDist = dist[v] + edgeWeight[he.edge()];
                if (vDist < dist[w]) {
                    dist[w]       = vDist;
                    primalTree[w] = he;

                    // std::priority_queue is a max-heap, so we sort by
                    // negative cost to pop the closest vertex next
                    toVisit.push({-dist[w], w});
                }
            }
        }
    }

    for (Vertex v : mesh.vertices()) {
        if (!visited[v]) {
            std::cout << "Error: primal spanning tree did not reach all faces. "
                         "Is the input mesh connected?"
                      << vendl;
        }
    }

    return primalTree;
}

VertexData<Halfedge>
buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                        const FaceAndBoundaryData<Halfedge>& dualTree,
                        bool connectBoundary, Vertex root) {
    EdgeData<double> edgeWeight(mesh, 1);
    return buildPrimalSpanningTree(mesh, edgeWeight, dualTree, connectBoundary,
                                   root);
}

VertexData<Halfedge>
buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight,
                        const FaceAndBoundaryData<Halfedge>& dualTree,
                        bool connectBoundary, Vertex root) {

    // The cost is the length of the loop that this edge makes when added to the
    // dual tree, or infinity for edges in the dual tree
    EdgeData<double> edgeCost(mesh);
    for (Edge e : mesh.edges()) {
        if (inDualSpanningCotree(e.halfedge(), dualTree)) {
            edgeCost[e] = std::numeric_limits<double>::infinity();
        } else {
            edgeCost[e] = 0;
            for (Halfedge he : dualTreeLoop(e.halfedge(), dualTree, false)) {
                edgeCost[e] += edgeWeight[he.edge()];
            }
        }
    }

    VertexData<Halfedge> primalTree(mesh, Halfedge());

    std::priority_queue<std::tuple<double, Vertex, Halfedge>> toVisit;
    if (connectBoundary) {
        // If we want to connect up the boundary vertices, we might as well
        // start there
        // This is necessary for the loop tracing code later
        for (BoundaryLoop b : mesh.boundaryLoops()) {
            for (Vertex v : b.adjacentVertices()) {
                toVisit.push(
                    {std::numeric_limits<double>::infinity(), v, Halfedge()});
            }
        }
    }
    if (toVisit.empty()) {
        if (root == Vertex()) root = mesh.vertex(0);
        toVisit.push({0, root, Halfedge()});
    } else {
        root = Vertex(); // ignore root if we have a boundary
    }

    while (!toVisit.empty()) {
        Vertex v          = std::get<1>(toVisit.top());
        Halfedge parentHe = std::get<2>(toVisit.top());
        toVisit.pop();

        // if we've already processed this vertex, skip it
        if (primalTree[v] != Halfedge()) continue;

        primalTree[v] = parentHe;

        for (Halfedge he : v.outgoingHalfedges()) {
            if (!std::isfinite(edgeCost[he.edge()])) continue;
            Vertex w = he.tipVertex();

            if (w != root && !(connectBoundary && w.isBoundary()) &&
                primalTree[w] == Halfedge()) {
                toVisit.push({edgeCost[he.edge()], w, he});
            }
        }
    }

    for (Vertex v : mesh.vertices()) {
        if (!v.isBoundary() && v != root && primalTree[v] == Halfedge()) {
            std::cout
                << "Error: primal spanning tree did not reach all vertices. "
                   "Is the input mesh connected?"
                << vendl;
            std::cout << "\t" << v << vendl;
        }
    }

    return primalTree;
}

bool inPrimalSpanningTree(Halfedge he, const VertexData<Halfedge>& primalTree) {
    Vertex v = he.tailVertex();
    Vertex w = he.tipVertex();
    return (primalTree[v] != Halfedge() && primalTree[v].edge() == he.edge()) ||
           (primalTree[w] != Halfedge() && primalTree[w].edge() == he.edge());
}

FaceAndBoundaryData<Halfedge> buildDualSpanningCotree(ManifoldSurfaceMesh& mesh,
                                                      bool connectBoundary,
                                                      Face root) {
    VertexData<Halfedge> emptyPrimalTree;
    return buildDualSpanningCotree(mesh, emptyPrimalTree, connectBoundary,
                                   root);
}

FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight,
                        bool connectBoundary, Face root) {

    FaceAndBoundaryData<Halfedge> dualTree(mesh, Halfedge());

    // Distance estimates to build a shortest-paths tree
    FaceAndBoundaryData<double> dist(mesh,
                                     std::numeric_limits<double>::infinity());

    FaceAndBoundaryData<bool> visited(mesh, false);

    std::priority_queue<std::tuple<double, Face>> toVisit;

    // make sure we start on a boundary face if connectBoundary is true
    // this is necessary for the loop tracing code later
    if (connectBoundary && mesh.nBoundaryLoops() > 0) {
        root = mesh.boundaryLoop(0).asFace();
    }
    if (root == Face()) root = mesh.face(0);
    toVisit.push({0, root});
    dist[root]     = 0;
    dualTree[root] = Halfedge();

    while (!toVisit.empty()) {
        Face f = std::get<1>(toVisit.top());
        toVisit.pop();

        // if we've already processed this face, skip it
        if (visited[f]) continue;
        visited[f] = true;

        if (connectBoundary && f.isBoundaryLoop()) {

            for (BoundaryLoop b : mesh.boundaryLoops()) {
                visited[b] = true;
                dist[b]    = dist[f];
                for (Halfedge he : b.adjacentHalfedges()) {
                    Face g = he.twin().face();
                    if (visited[g] || g.isBoundaryLoop()) continue;

                    double fDist = dist[f] + edgeWeight[he.edge()];
                    if (fDist < dist[g]) {
                        dist[g]     = fDist;
                        dualTree[g] = he;

                        // std::priority_queue is a max-heap, so we sort by
                        // negative cost to pop the closest face next
                        toVisit.push({-dist[g], g});
                    }
                }
            }
        } else {

            for (Halfedge he : f.adjacentHalfedges()) {
                Face g = he.twin().face();
                if (visited[g] || (!connectBoundary && g.isBoundaryLoop()))
                    continue;

                double fDist = dist[f] + edgeWeight[he.edge()];
                if (fDist < dist[g]) {
                    dist[g]     = fDist;
                    dualTree[g] = he;

                    // std::priority_queue is a max-heap, so we sort by
                    // negative cost to pop the closest face next
                    toVisit.push({-dist[g], g});
                }
            }
        }
    }

    for (Face f : mesh.faces()) {
        if (!visited[f]) {
            std::cout << "Error: dual spanning tree did not reach all faces. "
                         "Is the input mesh connected?"
                      << vendl;
        }
    }

    return dualTree;
}

FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh,
                        const VertexData<Halfedge>& primalTree,
                        bool connectBoundary, Face root) {
    EdgeData<double> edgeWeight(mesh, 1);
    return buildDualSpanningCotree(mesh, edgeWeight, primalTree,
                                   connectBoundary, root);
}

FaceAndBoundaryData<Halfedge> buildDualSpanningCotree(
    ManifoldSurfaceMesh& mesh, const EdgeData<double>& edgeWeight,
    const VertexData<Halfedge>& primalTree, bool connectBoundary, Face root) {

    // The cost is the length of the loop that this edge makes when added to the
    // primal tree, or infinity for edges in the primal tree
    EdgeData<double> edgeCost(mesh);
    for (Edge e : mesh.edges()) {
        if (inPrimalSpanningTree(e.halfedge(), primalTree)) {
            edgeCost[e] = std::numeric_limits<double>::infinity();
        } else {
            edgeCost[e] = 0;
            for (Halfedge he :
                 primalTreeLoop(e.halfedge(), primalTree, false)) {
                edgeCost[e] += edgeWeight[he.edge()];
            }
        }
    }

    FaceAndBoundaryData<Halfedge> dualTree(mesh, Halfedge());
    bool foundBoundary = false;

    // make sure we start on a boundary face if connectBoundary is true
    // this is necessary for the loop tracing code later
    if (connectBoundary && mesh.nBoundaryLoops() > 0) {
        root = mesh.boundaryLoop(0).asFace();
    }
    if (root == Face()) root = mesh.face(0);
    std::priority_queue<std::tuple<double, Face, Halfedge>> toVisit;
    toVisit.push({0, root, Halfedge()});

    while (!toVisit.empty()) {
        Face f            = std::get<1>(toVisit.top());
        Halfedge parentHe = std::get<2>(toVisit.top());
        toVisit.pop();

        if (connectBoundary && f.isBoundaryLoop()) {
            // if we've already processed the boundary faces, skip this one
            if (foundBoundary) continue;

            // Special case for boundary : loop over all boundary
            // halfedges, not just those on this loop
            foundBoundary = true;

            for (BoundaryLoop b : mesh.boundaryLoops()) {
                dualTree[b] = parentHe;
                for (Halfedge he : b.adjacentHalfedges()) {
                    if (!std::isfinite(edgeCost[he.edge()])) continue;
                    Face g = he.twin().face();

                    if (!g.isBoundaryLoop() && g != root &&
                        dualTree[g] == Halfedge()) {
                        toVisit.push({edgeCost[he.edge()], g, he});
                    }
                }
            }
        } else if (!f.isBoundaryLoop()) {
            // if we've already processed this face, skip it
            if (dualTree[f] != Halfedge()) continue;

            dualTree[f] = parentHe;

            for (Halfedge he : f.adjacentHalfedges()) {
                if (!std::isfinite(edgeCost[he.edge()])) continue;
                Face g = he.twin().face();

                if (!g.isBoundaryLoop() && g != root &&
                    dualTree[g] == Halfedge()) {
                    toVisit.push({edgeCost[he.edge()], g, he});
                } else if (g.isBoundaryLoop() && !foundBoundary) {
                    toVisit.push({edgeCost[he.edge()], g, he});
                }
            }
        }
    }

    for (Face f : mesh.faces()) {
        if (f != root && dualTree[f] == Halfedge()) {
            std::cout << "Error: dual spanning tree did not reach all faces. "
                         "Is the input mesh connected?"
                      << vendl;
        }
    }

    return dualTree;
}

bool inDualSpanningCotree(Halfedge he,
                          const FaceAndBoundaryData<Halfedge>& dualTree) {
    Face f = he.face();
    Face g = he.twin().face();
    return (dualTree[f] != Halfedge() && dualTree[f].edge() == he.edge()) ||
           (dualTree[g] != Halfedge() && dualTree[g].edge() == he.edge());
}

std::vector<Halfedge> primalTreeLoop(Halfedge he,
                                     const VertexData<Halfedge>& primalTree,
                                     bool trimLoop) {

    // Walk up along primal tree to extract generator
    Halfedge curr = primalTree[he.tipVertex()];
    std::vector<Halfedge> forwardPath{he};
    while (curr != Halfedge()) {
        forwardPath.push_back(curr.twin());
        curr = primalTree[curr.tailVertex()];
    }

    curr = he;
    std::vector<Halfedge> backwardPath;
    while (primalTree[curr.tailVertex()] != Halfedge()) {
        curr = primalTree[curr.tailVertex()];
        backwardPath.push_back(curr);
    }

    if (trimLoop && !forwardPath.empty() && !backwardPath.empty()) {
        // Identify and remove common edges at the ends of both paths
        size_t nShared = 0;
        size_t nF      = forwardPath.size() - 1;
        size_t nB      = backwardPath.size() - 1;
        while (nShared <= nF && nShared <= nB &&
               forwardPath[nF - nShared] == backwardPath[nB - nShared].twin()) {
            nShared++;
        }

        // Remove last nShared elements from paths
        // https://stackoverflow.com/questions/34452139/how-to-remove-several-elements-from-the-end-of-stdvector
        forwardPath.resize(forwardPath.size() - nShared);
        backwardPath.resize(backwardPath.size() - nShared);
    }

    std::vector<Halfedge> loop;
    // First go backwards along backwardPath
    loop.insert(loop.end(), backwardPath.rbegin(), backwardPath.rend());
    // Then go forwards along forwardPath
    loop.insert(loop.end(), forwardPath.begin(), forwardPath.end());
    return loop;
}

std::vector<Halfedge>
dualTreeLoop(Halfedge he, const FaceAndBoundaryData<Halfedge>& dualTree,
             bool trimLoop) {
    std::vector<Halfedge> path =
        dualTreePath(he.twin(), he, dualTree, trimLoop);
    path.push_back(he);
    return path;
}

std::vector<Halfedge>
dualTreePath(Halfedge heSrc, Halfedge heDst,
             const FaceAndBoundaryData<Halfedge>& dualTree, bool trimLoop) {
    // Walk up along primal tree to extract generator
    Halfedge curr = dualTree[heSrc.face()];
    std::vector<Halfedge> forwardPath{};
    if (curr != Halfedge()) {
        forwardPath.push_back(curr.twin());
        while (dualTree[curr.face()] != Halfedge()) {
            curr = dualTree[curr.face()];
            forwardPath.push_back(curr.twin());
        }
    }

    curr = dualTree[heDst.face()];
    std::vector<Halfedge> backwardPath;
    while (curr != Halfedge()) {
        backwardPath.push_back(curr);
        curr = dualTree[curr.face()];
    }

    if (trimLoop && !forwardPath.empty() && !backwardPath.empty()) {
        // Identify and remove common edges at the ends of both paths
        size_t nShared = 0;
        size_t nF      = forwardPath.size() - 1;
        size_t nB      = backwardPath.size() - 1;
        while (nShared <= nF && nShared <= nB &&
               forwardPath[nF - nShared] == backwardPath[nB - nShared].twin()) {
            nShared++;
        }

        // Remove last nShared elements from paths
        // https://stackoverflow.com/questions/34452139/how-to-remove-several-elements-from-the-end-of-stdvector
        forwardPath.resize(forwardPath.size() - nShared);
        backwardPath.resize(backwardPath.size() - nShared);
    }

    std::vector<Halfedge> loop;
    // Then go forwards along forwardPath
    loop.insert(loop.end(), forwardPath.begin(), forwardPath.end());
    // First go backwards along backwardPath
    loop.insert(loop.end(), backwardPath.rbegin(), backwardPath.rend());
    return loop;
}
