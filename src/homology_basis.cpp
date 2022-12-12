#include "homology_basis.h"

std::vector<std::vector<Halfedge>>
computePrimalHomologyBasis(ManifoldSurfaceMesh& mesh) {
    bool hasBoundary = mesh.nBoundaryLoops() > 0;

    // In meshes with boundary, we should construct the dual tree first
    // In meshes without boundary, it doesn't matter
    // FaceAndBoundaryData<Halfedge> dualTree     =
    // buildDualSpanningCotree(mesh); VertexData<Halfedge> primalTree =
    // buildPrimalSpanningTree(mesh, dualTree);
    VertexData<Halfedge> primalTree = buildPrimalSpanningTree(mesh);
    FaceAndBoundaryData<Halfedge> dualTree =
        buildDualSpanningCotree(mesh, primalTree);

    std::vector<std::vector<Halfedge>> primalGenerators;
    // bool first = true; // in meshes with boundary, we ignore one boundary
    // edge
    for (Edge e : mesh.edges()) {
        if (!inPrimalSpanningTree(e.halfedge(), primalTree) &&
            !inDualSpanningCotree(e.halfedge(), dualTree)) {
            // if (hasBoundary && first) {
            //     first = false;
            //     continue;
            // }
            // first = false;
            primalGenerators.push_back(
                primalTreeLoop(e.halfedge(), primalTree));
        }
    }

    return primalGenerators;
}

// TODO: computing dual tree first messes with greedy strategy
std::vector<std::vector<Halfedge>>
computePrimalHomologyBasis(ManifoldSurfaceMesh& mesh,
                           const EdgeData<double>& edgeCost) {
    bool hasBoundary = mesh.nBoundaryLoops() > 0;

    // In meshes with boundary, we should construct the dual tree first
    // In meshes without boundary, it doesn't matter
    VertexData<Halfedge> primalTree = buildPrimalSpanningTree(mesh, edgeCost);
    FaceAndBoundaryData<Halfedge> dualTree =
        buildDualSpanningCotree(mesh, edgeCost, primalTree);

    std::vector<std::vector<Halfedge>> primalGenerators;
    // bool first = true; // in meshes with boundary, we ignore one boundary
    // edge
    for (Edge e : mesh.edges()) {
        if (!inPrimalSpanningTree(e.halfedge(), primalTree) &&
            !inDualSpanningCotree(e.halfedge(), dualTree)) {
            // if (hasBoundary && first) {
            //     first = false;
            //     continue;
            // }
            // first = false;
            primalGenerators.push_back(
                primalTreeLoop(e.halfedge(), primalTree));
        }
    }

    return primalGenerators;
}

// TODO: dual tree cotree on mesh with boundary
std::vector<std::vector<Halfedge>>
computeDualHomologyBasis(ManifoldSurfaceMesh& mesh) {
    // order doesn't matter in a mesh without boundary
    VertexData<Halfedge> primalTree = buildPrimalSpanningTree(mesh);
    FaceAndBoundaryData<Halfedge> dualTree =
        buildDualSpanningCotree(mesh, primalTree);

    std::vector<std::vector<Halfedge>> dualGenerators;
    for (Edge e : mesh.edges()) {
        if (!e.isBoundary() &&
            !inPrimalSpanningTree(e.halfedge(), primalTree) &&
            !inDualSpanningCotree(e.halfedge(), dualTree)) {
            dualGenerators.push_back(dualTreeLoop(e.halfedge(), dualTree));
        }
    }
    return dualGenerators;
}

std::tuple<std::vector<std::vector<Halfedge>>,
           std::vector<std::vector<Halfedge>>>
computePrimalAndDualHomologyBases(ManifoldSurfaceMesh& mesh) {
    bool hasBoundary = mesh.nBoundaryLoops() > 0;

    // In meshes with boundary, we should construct the dual tree first
    // In meshes without boundary, it doesn't matter
    VertexData<Halfedge> primalTree = buildPrimalSpanningTree(mesh);
    FaceAndBoundaryData<Halfedge> dualTree =
        buildDualSpanningCotree(mesh, primalTree);

    std::vector<std::vector<Halfedge>> primalGenerators, dualGenerators;
    // in meshes with boundary, we ignore one boundary edge
    Halfedge bdyHalfedge;
    for (Edge e : mesh.edges()) {
        if (!inPrimalSpanningTree(e.halfedge(), primalTree) &&
            !inDualSpanningCotree(e.halfedge(), dualTree)) {

            // if (hasBoundary && bdyHalfedge == Halfedge()) {
            //     bdyHalfedge = e.halfedge();
            //     continue;
            // }

            primalGenerators.push_back(
                primalTreeLoop(e.halfedge(), primalTree));

            // if (e.isBoundary()) {
            //     dualGenerators.push_back(
            //         dualTreePath(e.halfedge(), bdyHalfedge, dualTree));
            // } else {
            dualGenerators.push_back(dualTreeLoop(e.halfedge(), dualTree));
            // }
        }
    }
    return std::make_tuple(primalGenerators, dualGenerators);
}

std::tuple<std::vector<std::vector<Halfedge>>,
           std::vector<std::vector<Halfedge>>>
computePrimalAndDualHomologyBases(ManifoldSurfaceMesh& mesh,
                                  const EdgeData<double>& edgeCost) {
    bool hasBoundary = mesh.nBoundaryLoops() > 0;

    // In meshes with boundary, we should construct the dual tree first
    // In meshes without boundary, it doesn't matter
    FaceAndBoundaryData<Halfedge> dualTree =
        buildDualSpanningCotree(mesh, edgeCost);
    VertexData<Halfedge> primalTree =
        buildPrimalSpanningTree(mesh, edgeCost, dualTree);

    std::vector<std::vector<Halfedge>> primalGenerators, dualGenerators;
    // in meshes with boundary, we ignore one boundary edge
    Halfedge bdyHalfedge;
    for (Edge e : mesh.edges()) {
        if (!inPrimalSpanningTree(e.halfedge(), primalTree) &&
            !inDualSpanningCotree(e.halfedge(), dualTree)) {

            if (hasBoundary && bdyHalfedge == Halfedge()) {
                bdyHalfedge = e.halfedge();
                continue;
            }

            primalGenerators.push_back(
                primalTreeLoop(e.halfedge(), primalTree));

            if (e.isBoundary()) {
                dualGenerators.push_back(
                    dualTreePath(e.halfedge(), bdyHalfedge, dualTree));
            } else {
                dualGenerators.push_back(dualTreeLoop(e.halfedge(), dualTree));
            }
        }
    }
    return std::make_tuple(primalGenerators, dualGenerators);
}

VertexData<Halfedge> buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh) {
    FaceAndBoundaryData<Halfedge> emptyDualTree;
    return buildPrimalSpanningTree(mesh, emptyDualTree);
}

VertexData<Halfedge>
buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight) {

    VertexData<Halfedge> primalTree(mesh, Halfedge());

    // Distance estimates to build a shortest-paths tree
    VertexData<double> dist(mesh, std::numeric_limits<double>::infinity());

    // Most vertices not in the tree can be identified because
    // primalTree[v] == Halfedge(). But this is not true for the root, so we
    // keep track of this data in a separate array
    VertexData<bool> inTree(mesh, false);

    // Loop over all vertices to build a spanning forest on disconnected meshes
    for (Vertex v : mesh.vertices()) {
        // if (!inTree[v] && !v.isBoundary()) {
        if (!inTree[v]) {
            inTree[v] = true;

            Vertex root = v;
            std::priority_queue<std::pair<double, Vertex>> toVisit;
            toVisit.push({0, root});
            dist[root] = 0;

            while (!toVisit.empty()) {
                Vertex v = std::get<1>(toVisit.top());
                toVisit.pop();

                for (Halfedge he : v.outgoingHalfedges()) {
                    Vertex w = he.twin().vertex();
                    dist[w]  = fmin(dist[w], dist[v] + edgeWeight[he.edge()]);
                    if (w != root && primalTree[w] == Halfedge()) {
                        primalTree[w] = he;
                        // std::priority_queue is a max-heap, so we sort by
                        // negative cost to pop the closest vertex next
                        toVisit.push({-dist[w], w});
                        inTree[w] = true;
                    }
                }
            }
        }
    }
    return primalTree;
}

VertexData<Halfedge>
buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                        const FaceAndBoundaryData<Halfedge>& dualTree) {
    auto mayCross = [&](Halfedge he) -> bool {
        return dualTree.size() == 0 || !inDualSpanningCotree(he, dualTree);
    };

    VertexData<Halfedge> primalTree(mesh, Halfedge());

    // Most vertices not in the tree can be identified because
    // primalTree[v] == Halfedge(). But this is not true for the root, so we
    // keep track of this data in a separate array
    VertexData<bool> inTree(mesh, false);

    // Loop over all vertices to build a spanning forest on disconnected meshes
    for (Vertex v : mesh.vertices()) {
        if (!inTree[v]) {
            inTree[v] = true;

            Vertex root = v;
            std::queue<Vertex> toVisit;
            toVisit.push(root);

            while (!toVisit.empty()) {
                Vertex v = toVisit.front();
                toVisit.pop();

                for (Halfedge he : v.outgoingHalfedges()) {
                    if (!mayCross(he)) continue;
                    Vertex w = he.twin().vertex();
                    // if (w.isBoundary()) continue;
                    if (w != root && primalTree[w] == Halfedge()) {
                        primalTree[w] = he;
                        toVisit.push(w);
                        inTree[w] = true;
                    }
                }
            }
        }
    }
    return primalTree;
}

VertexData<Halfedge>
buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight,
                        const FaceAndBoundaryData<Halfedge>& dualTree) {

    auto mayCross = [&](Halfedge he) -> bool {
        return dualTree.size() == 0 || !inDualSpanningCotree(he, dualTree);
    };

    EdgeData<double> edgeCost(mesh);
    //  the cost of an edge is the length of the loop that this edge makes when
    //  added to the dual tree
    for (Edge e : mesh.edges()) {
        if (!inDualSpanningCotree(e.halfedge(), dualTree)) {
            for (Halfedge he : dualTreeLoop(e.halfedge(), dualTree)) {
                edgeCost[e] += edgeWeight[he.edge()];
            }
        }
    }

    VertexData<Halfedge> primalTree(mesh, Halfedge());

    // Most vertices not in the tree can be identified because
    // primalTree[v] == Halfedge(). But this is not true for the root, so we
    // keep track of this data in a separate array
    VertexData<bool> inTree(mesh, false);

    // Loop over all vertices to build a spanning forest on disconnected meshes
    for (Vertex v : mesh.vertices()) {
        // if (!inTree[v] && !v.isBoundary()) {
        if (!inTree[v]) {
            inTree[v] = true;

            Vertex root = v;
            std::priority_queue<std::pair<double, Vertex>> toVisit;
            toVisit.push({0, root});

            while (!toVisit.empty()) {
                Vertex v = std::get<1>(toVisit.top());
                toVisit.pop();

                for (Halfedge he : v.outgoingHalfedges()) {
                    if (!mayCross(he)) continue;
                    Vertex w = he.twin().vertex();
                    if (w != root && primalTree[w] == Halfedge()) {
                        primalTree[w] = he;
                        // std::priority_queue is a max-heap, so we sort by
                        // negative cost to pop the closest vertex next
                        toVisit.push({edgeCost[he.edge()], w});
                        inTree[w] = true;
                    }
                }
            }
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

FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh) {
    VertexData<Halfedge> emptyPrimalTree;
    return buildDualSpanningCotree(mesh, emptyPrimalTree);
}

FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight) {
    FaceAndBoundaryData<Halfedge> dualTree(mesh, Halfedge());

    // Most faces not in the tree can be identified because
    // dualTree[f] == Halfedge(). But this is not true for the root, so we
    // keep track of this data in a separate array
    FaceAndBoundaryData<bool> inTree(mesh, false);

    // Distance estimates to build a shortest-paths tree
    FaceAndBoundaryData<double> dist(mesh,
                                     std::numeric_limits<double>::infinity());

    // Loop over all faces to build a spanning forest on disconnected meshes
    for (Face f : mesh.faces()) {
        if (!inTree[f]) {
            inTree[f] = true;

            Face root = f;
            std::priority_queue<std::pair<double, Face>> toVisit;
            toVisit.push({0, root});
            dist[root] = 0;

            while (!toVisit.empty()) {
                Face f = std::get<1>(toVisit.top());
                toVisit.pop();

                for (Halfedge he : f.adjacentHalfedges()) {
                    Face g  = he.twin().face();
                    dist[g] = fmin(dist[g], dist[f] + edgeWeight[he.edge()]);
                    if (g != root && dualTree[g] == Halfedge()) {
                        dualTree[g] = he;
                        // std::priority_queue is a max-heap, so we sort by
                        // negative cost to pop the closest vertex next
                        toVisit.push({-dist[g], g});
                        inTree[g] = true;
                    }
                }
            }
        }
    }
    return dualTree;
}

FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh,
                        const VertexData<Halfedge>& primalTree) {

    auto mayCross = [&](Halfedge he) -> bool {
        return primalTree.size() == 0 || !inPrimalSpanningTree(he, primalTree);
    };

    FaceAndBoundaryData<Halfedge> dualTree(mesh, Halfedge());

    // Most faces not in the tree can be identified because
    // dualTree[f] == Halfedge(). But this is not true for the root, so we
    // keep track of this data in a separate array
    FaceData<bool> inTree(mesh, false);

    bool foundBoundary = false;

    // TODO: connect boundary edges
    // Loop over all faces to build a spanning forest on disconnected meshes
    for (Face f : mesh.faces()) {
        if (!inTree[f]) {
            inTree[f] = true;

            Face root = f;
            std::queue<Face> toVisit;
            toVisit.push(root);

            while (!toVisit.empty()) {
                Face f = toVisit.front();
                toVisit.pop();

                if (f.isBoundaryLoop()) {
                    // Special case for boundary : loop over all boundary
                    // halfedges, not just those on this loop

                    for (BoundaryLoop b : mesh.boundaryLoops()) {
                        for (Halfedge he : b.adjacentHalfedges()) {
                            if (!mayCross(he)) continue;

                            Face g = he.twin().face();
                            if (g.isBoundaryLoop()) {
                                continue;
                            } else if (g != root && dualTree[g] == Halfedge()) {
                                dualTree[g] = he;
                                toVisit.push(g);
                                inTree[g] = true;
                            }
                        }
                    }

                } else {
                    for (Halfedge he : f.adjacentHalfedges()) {
                        if (!mayCross(he)) continue;

                        Face g = he.twin().face();
                        if (g.isBoundaryLoop() && !foundBoundary) {
                            foundBoundary = true;
                            dualTree[g]   = he;
                            toVisit.push(g);
                        } else if (g != root && dualTree[g] == Halfedge()) {
                            dualTree[g] = he;
                            toVisit.push(g);
                            inTree[g] = true;
                        }
                    }
                }
            }
        }

        // TODO: delete this
        break;
    }
    return dualTree;
}

FaceAndBoundaryData<Halfedge>
buildDualSpanningCotree(ManifoldSurfaceMesh& mesh,
                        const EdgeData<double>& edgeWeight,
                        const VertexData<Halfedge>& primalTree) {

    auto mayCross = [&](Halfedge he) -> bool {
        return primalTree.size() == 0 || !inPrimalSpanningTree(he, primalTree);
    };

    EdgeData<double> edgeCost(mesh);
    // The cost is the length of the loop that this edge makes when added to the
    // primal tree
    for (Edge e : mesh.edges()) {
        if (!inPrimalSpanningTree(e.halfedge(), primalTree)) {
            for (Halfedge he : primalTreeLoop(e.halfedge(), primalTree)) {
                edgeCost[e] += edgeWeight[he.edge()];
            }
        }
    }

    FaceAndBoundaryData<Halfedge> dualTree(mesh, Halfedge());

    // Most faces not in the tree can be identified because
    // dualTree[f] == Halfedge(). But this is not true for the root, so we
    // keep track of this data in a separate array
    FaceAndBoundaryData<bool> inTree(mesh, false);

    // TODO: connect boundary edges
    // Loop over all faces to build a spanning forest on disconnected meshes
    for (Face f : mesh.faces()) {
        if (!inTree[f]) {
            inTree[f] = true;

            Face root = f;
            std::priority_queue<std::pair<double, Face>> toVisit;
            toVisit.push({0, root});

            while (!toVisit.empty()) {
                Face f = std::get<1>(toVisit.top());
                toVisit.pop();

                for (Halfedge he : f.adjacentHalfedges()) {
                    if (!mayCross(he)) continue;

                    Face g = he.twin().face();
                    if (g != root && dualTree[g] == Halfedge()) {
                        dualTree[g] = he;
                        toVisit.push({edgeCost[he.edge()], g});
                        inTree[g] = true;
                    }
                }
            }
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
                                     const VertexData<Halfedge>& primalTree) {

    // Walk up along primal tree to extract generator
    Halfedge curr = primalTree[he.tipVertex()];
    std::vector<Halfedge> forwardPath{he};
    if (curr != Halfedge()) {
        forwardPath.push_back(curr.twin());
        while (primalTree[curr.tailVertex()] != Halfedge()) {
            curr = primalTree[curr.tailVertex()];
            forwardPath.push_back(curr.twin());
        }
    }

    curr = he;
    std::vector<Halfedge> backwardPath;
    while (primalTree[curr.tailVertex()] != Halfedge()) {
        curr = primalTree[curr.tailVertex()];
        backwardPath.push_back(curr);
    }

    // Identify and remove common edges at the ends of both paths
    size_t nShared = 0;
    size_t nF      = forwardPath.size() - 1;
    size_t nB      = backwardPath.size() - 1;
    while (nShared < nF && nShared < nB &&
           forwardPath[nF - nShared] == backwardPath[nB - nShared].twin()) {
        nShared++;
    }

    // Remove last nShared elements from paths
    // https://stackoverflow.com/questions/34452139/how-to-remove-several-elements-from-the-end-of-stdvector
    forwardPath.resize(forwardPath.size() - nShared);
    backwardPath.resize(backwardPath.size() - nShared);


    std::vector<Halfedge> loop;
    // First go backwards along backwardPath
    loop.insert(loop.end(), backwardPath.rbegin(), backwardPath.rend());
    // Then go forwards along forwardPath
    loop.insert(loop.end(), forwardPath.begin(), forwardPath.end());
    return loop;
}

std::vector<Halfedge>
dualTreeLoop(Halfedge he, const FaceAndBoundaryData<Halfedge>& dualTree) {
    std::vector<Halfedge> path = dualTreePath(he.twin(), he, dualTree);
    path.push_back(he);
    return path;
}

std::vector<Halfedge>
dualTreePath(Halfedge heSrc, Halfedge heDst,
             const FaceAndBoundaryData<Halfedge>& dualTree) {
    // Walk up along primal tree to extract generator
    Halfedge curr = dualTree[heSrc.face()];
    std::vector<Halfedge> forwardPath{};
    if (curr != Halfedge()) {
        forwardPath.push_back(curr.twin());
        do {
            curr = dualTree[curr.face()];
            forwardPath.push_back(curr.twin());
        } while (dualTree[curr.face()] != Halfedge());
    }

    curr = dualTree[heDst.face()];
    std::vector<Halfedge> backwardPath{curr};
    do {
        curr = dualTree[curr.face()];
        backwardPath.push_back(curr);
    } while (dualTree[curr.face()] != Halfedge());

    // Identify and remove common edges at the ends of both paths
    size_t nShared = 0;
    size_t nF      = forwardPath.size() - 1;
    size_t nB      = backwardPath.size() - 1;
    while (nShared < nF && nShared < nB &&
           forwardPath[nF - nShared] == backwardPath[nB - nShared].twin()) {
        nShared++;
    }

    // Remove last nShared elements from paths
    // https://stackoverflow.com/questions/34452139/how-to-remove-several-elements-from-the-end-of-stdvector
    forwardPath.resize(forwardPath.size() - nShared);
    backwardPath.resize(backwardPath.size() - nShared);


    std::vector<Halfedge> loop;
    // Then go forwards along forwardPath
    loop.insert(loop.end(), forwardPath.begin(), forwardPath.end());
    // First go backwards along backwardPath
    loop.insert(loop.end(), backwardPath.rbegin(), backwardPath.rend());
    return loop;
}
