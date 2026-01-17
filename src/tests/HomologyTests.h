#include "test_utils.h"

#include "geometrycentral/surface/homology_generators.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

class HomologyTest : public ::testing::Test {
  public:
    static std::vector<std::unique_ptr<ManifoldSurfaceMesh>> meshPtrs;
    static std::vector<std::unique_ptr<VertexPositionGeometry>> geomPtrs;

  protected:
    static void SetUpTestSuite() {
        {
            auto doubleTorus =
                readManifoldSurfaceMesh("../data/double-torus-small.obj");
            meshPtrs.emplace_back(std::move(std::get<0>(doubleTorus)));
            geomPtrs.emplace_back(std::move(std::get<1>(doubleTorus)));
        }
        {
            auto annulus = readManifoldSurfaceMesh("../data/annulus.obj");
            meshPtrs.emplace_back(std::move(std::get<0>(annulus)));
            geomPtrs.emplace_back(std::move(std::get<1>(annulus)));
        }
    }
};

std::vector<std::unique_ptr<ManifoldSurfaceMesh>> HomologyTest::meshPtrs;
std::vector<std::unique_ptr<VertexPositionGeometry>> HomologyTest::geomPtrs;

TEST_F(HomologyTest, LoopsClosed) {
    for (size_t iM = 0; iM < meshPtrs.size(); iM++) {
        ManifoldSurfaceMesh& mesh = *meshPtrs[iM];
        auto primalLoopClosed     = [&](const std::vector<Halfedge>& primalLoop,
                                    bool ignoreBoundary = false) -> bool {
            VertexData<int> loopBdy(mesh, 0);
            for (Halfedge ij : primalLoop) {
                loopBdy[ij.tailVertex()] -= 1;
                loopBdy[ij.tipVertex()] += 1;
            }
            for (Vertex i : mesh.vertices())
                if (loopBdy[i] != 0 && !(ignoreBoundary && i.isBoundary()))
                    return false;
            return true;
        };
        auto allPrimalLoopsClosed =
            [&](const std::vector<std::vector<Halfedge>>& primalLoops,
                bool ignoreBoundary = false) -> bool {
            for (const std::vector<Halfedge>& loop : primalLoops) {
                if (!primalLoopClosed(loop, ignoreBoundary)) return false;
            }
            return true;
        };
        auto dualLoopClosed = [&](const std::vector<Halfedge>& dualLoop,
                                  bool ignoreBoundary = false) -> bool {
            FaceData<int> loopBdy(mesh, 0);
            for (Halfedge ij : dualLoop) {
                if (!ij.face().isBoundaryLoop()) loopBdy[ij.face()] -= 1;
                if (!ij.twin().face().isBoundaryLoop())
                    loopBdy[ij.twin().face()] += 1;
            }
            for (Face f : mesh.faces()) {
                if (loopBdy[f] != 0) {
                    if (ignoreBoundary) { // check if f is boundary
                        bool isBdy = false;
                        for (Edge e : f.adjacentEdges()) {
                            if (e.isBoundary()) {
                                isBdy = true;
                                break;
                            }
                        }
                        if (!isBdy) return false;
                    } else {
                        return false;
                    }
                }
            }
            return true;
        };
        auto allDualLoopsClosed =
            [&](const std::vector<std::vector<Halfedge>>& dualLoops,
                bool ignoreBoundary = false) -> bool {
            for (const std::vector<Halfedge>& loop : dualLoops) {
                if (!dualLoopClosed(loop, ignoreBoundary)) return false;
            }
            return true;
        };

        HomologyGeneratorOptions opt = defaultHomologyGeneratorOptions;
        opt.generatorType            = HomologyGeneratorType::AbsolutePrimal;
        ASSERT_TRUE(allPrimalLoopsClosed(
            computeHomologyGenerators(mesh, opt).primalGenerators, false));

        opt.generatorType = HomologyGeneratorType::RelativePrimal;
        ASSERT_TRUE(allPrimalLoopsClosed(
            computeHomologyGenerators(mesh, opt).primalGenerators, true));

        opt.generatorType = HomologyGeneratorType::AbsolutePrimalRelativeDual;
        ASSERT_TRUE(allPrimalLoopsClosed(
            computeHomologyGenerators(mesh, opt).primalGenerators, false));
        ASSERT_TRUE(allDualLoopsClosed(
            computeHomologyGenerators(mesh, opt).dualGenerators, true));

        opt.generatorType = HomologyGeneratorType::RelativePrimalAbsoluteDual;
        ASSERT_TRUE(allPrimalLoopsClosed(
            computeHomologyGenerators(mesh, opt).primalGenerators, true));
        ASSERT_TRUE(allDualLoopsClosed(
            computeHomologyGenerators(mesh, opt).dualGenerators, false));

        opt.generatorType = HomologyGeneratorType::AbsoluteDual;
        ASSERT_TRUE(allDualLoopsClosed(
            computeHomologyGenerators(mesh, opt).dualGenerators, false));

        opt.generatorType = HomologyGeneratorType::RelativeDual;
        ASSERT_TRUE(allDualLoopsClosed(
            computeHomologyGenerators(mesh, opt).dualGenerators, true));
    }
}

TEST_F(HomologyTest, HomologyBasisDuality) {
    for (size_t iM = 0; iM < meshPtrs.size(); iM++) {
        ManifoldSurfaceMesh& mesh = *(meshPtrs[iM]);

        // https://math.stackexchange.com/q/361121
        // |V| - |E| + |F| = 2 - 2g - b
        // we get 2g + (b-1) harmonic forms, unless b = 0, in which case there
        // are 2g
        int dim = 2 + mesh.nEdges() - mesh.nVertices() - mesh.nFaces() -
                  mesh.nBoundaryLoops();
        if (mesh.nBoundaryLoops() > 0) dim += mesh.nBoundaryLoops() - 1;

        auto intersectionProduct =
            [&](const std::vector<Halfedge>& primalLoop,
                const std::vector<Halfedge>& dualLoop) -> int {
            int result = 0;
            for (Halfedge hP : primalLoop) {
                for (Halfedge hD : dualLoop) {
                    if (hP == hD) {
                        result -= 1;
                    } else if (hP == hD.twin()) {
                        result += 1;
                    }
                }
            }
            return result;
        };

        auto basesAreDual =
            [&](const std::vector<std::vector<Halfedge>>& primalLoops,
                const std::vector<std::vector<Halfedge>>& dualLoops) -> bool {
            for (size_t iP = 0; iP < primalLoops.size(); iP++) {
                for (size_t iD = 0; iD < dualLoops.size(); iD++) {
                    int product =
                        intersectionProduct(primalLoops[iP], dualLoops[iD]);
                    if ((iP == iD) && (product != 1)) return false;
                    if ((iP != iD) && (product != 0)) return false;
                }
            }
            return true;
        };

        HomologyGeneratorOptions opt = defaultHomologyGeneratorOptions;
        opt.generatorType = HomologyGeneratorType::AbsolutePrimalRelativeDual;
        HomologyGenerators gen = computeHomologyGenerators(mesh, opt);
        ASSERT_TRUE(basesAreDual(gen.primalGenerators, gen.dualGenerators));
        ASSERT_EQ(gen.primalGenerators.size(), dim);
        ASSERT_EQ(gen.dualGenerators.size(), dim);

        opt.generatorType = HomologyGeneratorType::RelativePrimalAbsoluteDual;
        gen               = computeHomologyGenerators(mesh, opt);
        ASSERT_TRUE(basesAreDual(gen.primalGenerators, gen.dualGenerators));
        ASSERT_EQ(gen.primalGenerators.size(), dim);
        ASSERT_EQ(gen.dualGenerators.size(), dim);
    }
}

TEST_F(HomologyTest, DualHomologyBasis) {
    for (size_t iM = 0; iM < meshPtrs.size(); iM++) {
        ManifoldSurfaceMesh& mesh = *(meshPtrs[iM]);

        auto intersectionProduct =
            [&](const std::vector<Halfedge>& primalLoop,
                const std::vector<Halfedge>& dualLoop) -> int {
            int result = 0;
            for (Halfedge hP : primalLoop) {
                for (Halfedge hD : dualLoop) {
                    if (hP == hD) {
                        result -= 1;
                    } else if (hP == hD.twin()) {
                        result += 1;
                    }
                }
            }
            return result;
        };

        auto basesAreDual =
            [&](const std::vector<std::vector<Halfedge>>& primalLoops,
                const std::vector<std::vector<Halfedge>>& dualLoops) -> bool {
            for (size_t iP = 0; iP < primalLoops.size(); iP++) {
                for (size_t iD = 0; iD < dualLoops.size(); iD++) {
                    int product =
                        intersectionProduct(primalLoops[iP], dualLoops[iD]);
                    if ((iP == iD) && (product != 1)) return false;
                    if ((iP != iD) && (product != 0)) return false;
                }
            }
            return true;
        };

        HomologyGeneratorOptions opt = defaultHomologyGeneratorOptions;
        opt.generatorType = HomologyGeneratorType::AbsolutePrimalRelativeDual;
        HomologyGenerators gen = computeHomologyGenerators(mesh, opt);
        ASSERT_TRUE(basesAreDual(gen.primalGenerators, gen.dualGenerators));

        opt.generatorType = HomologyGeneratorType::RelativePrimalAbsoluteDual;
        gen               = computeHomologyGenerators(mesh, opt);
        ASSERT_TRUE(basesAreDual(gen.primalGenerators, gen.dualGenerators));
    }
}

TEST_F(HomologyTest, HarmonicBasis) {
    for (size_t iM = 0; iM < meshPtrs.size(); iM++) {
        ManifoldSurfaceMesh& mesh    = *(meshPtrs[iM]);
        VertexPositionGeometry& geom = *(geomPtrs[iM]);

        // https://math.stackexchange.com/q/361121
        // |V| - |E| + |F| = 2 - 2g - b
        // we get 2g + (b-1) harmonic forms, unless b = 0, in which case there
        // are 2g
        int dim = 2 + mesh.nEdges() - mesh.nVertices() - mesh.nFaces() -
                  mesh.nBoundaryLoops();
        if (mesh.nBoundaryLoops() > 0) dim += mesh.nBoundaryLoops() - 1;

        geom.requireDECOperators();
        const SparseMatrix<double>& d0        = geom.d0;
        const SparseMatrix<double>& d1        = geom.d1;
        const SparseMatrix<double>& hodge1    = geom.hodge1;
        const SparseMatrix<double>& hodge1Inv = geom.hodge1Inverse;

        auto sign = [](Halfedge ij) -> double {
            return ij.orientation() ? 1 : -1;
        };

        auto primalFormsAreClosed =
            [&](const std::vector<EdgeData<double>>& primalForms) -> bool {
            for (const EdgeData<double>& form : primalForms) {
                if ((d1 * form.raw()).norm() > 1e-8) return false;
            }
            return true;
        };
        auto primalFormsAreCoclosed =
            [&](const std::vector<EdgeData<double>>& primalForms,
                bool ignoreBoundary = false) -> bool {
            for (const EdgeData<double>& form : primalForms) {
                Vector<double> delForm = d0.transpose() * hodge1 * form.raw();
                if ((delForm).norm() > 1e-8) {
                    if (ignoreBoundary) {
                        for (Vertex v : mesh.vertices()) {
                            if (!v.isBoundary() &&
                                abs(delForm(v.getIndex())) > 1e-5)
                                return false;
                        }
                    } else {
                        return false;
                    }
                }
            }
            return true;
        };

        auto dualFormsAreClosed =
            [&](const std::vector<EdgeData<double>>& dualForms,
                bool ignoreBoundary = false) -> bool {
            for (const EdgeData<double>& form : dualForms) {
                Vector<double> dForm = d0.transpose() * form.raw();
                if ((dForm).norm() > 1e-8) {
                    if (ignoreBoundary) {
                        for (Vertex v : mesh.vertices()) {
                            if (!v.isBoundary() &&
                                abs(dForm(v.getIndex())) > 1e-5)
                                return false;
                        }
                    } else {
                        return false;
                    }
                }
            }
            return true;
        };
        auto dualFormsAreCoclosed =
            [&](const std::vector<EdgeData<double>>& dualForms) -> bool {
            for (const EdgeData<double>& form : dualForms) {
                if ((d1 * hodge1Inv * form.raw()).norm() > 1e-8) return false;
            }
            return true;
        };

        auto integrateForm = [&](const EdgeData<double>& form,
                                 const std::vector<Halfedge>& path) -> double {
            double result = 0;
            for (Halfedge ij : path) result += sign(ij) * form[ij.edge()];
            return result;
        };

        auto basesAreDual =
            [&](const std::vector<EdgeData<double>>& forms,
                const std::vector<std::vector<Halfedge>>& loops) -> bool {
            for (size_t iF = 0; iF < forms.size(); iF++) {
                for (size_t iL = 0; iL < loops.size(); iL++) {
                    double integral = integrateForm(forms[iF], loops[iL]);
                    if ((iF == iL) && (fabs(integral - 1.) > 1e-5))
                        return false;
                    if ((iF != iL) && (fabs(integral - 0.) > 1e-5))
                        return false;
                }
            }
            return true;
        };

        HomologyGeneratorOptions opt = defaultHomologyGeneratorOptions;
        opt.generatorType = HomologyGeneratorType::AbsolutePrimalRelativeDual;
        HomologyGenerators loops = computeHomologyGenerators(mesh, opt);
        HarmonicGenerators forms = computeHarmonicGenerators(mesh, geom, loops);

        ASSERT_TRUE(
            basesAreDual(forms.primalGenerators, loops.primalGenerators));
        {
            for (size_t iF = 0; iF < forms.dualGenerators.size(); iF++) {
                for (size_t iL = 0; iL < loops.dualGenerators.size(); iL++) {
                    double integral = integrateForm(forms.dualGenerators[iF],
                                                    loops.dualGenerators[iL]);
                    if (iF == iL) EXPECT_NEAR(integral, 1, 1e-8);
                    if (iF != iL) EXPECT_NEAR(integral, 0, 1e-8);
                }
            }
        }
        EXPECT_TRUE(basesAreDual(forms.dualGenerators, loops.dualGenerators));
        EXPECT_TRUE(primalFormsAreClosed(forms.primalGenerators));
        EXPECT_TRUE(primalFormsAreCoclosed(forms.primalGenerators, false));
        EXPECT_TRUE(dualFormsAreClosed(forms.dualGenerators, false));
        EXPECT_TRUE(dualFormsAreCoclosed(forms.dualGenerators));
        EXPECT_EQ(forms.primalGenerators.size(), dim);
        EXPECT_EQ(loops.primalGenerators.size(), dim);

        opt.generatorType = HomologyGeneratorType::RelativePrimalAbsoluteDual;
        loops             = computeHomologyGenerators(mesh, opt);
        forms             = computeHarmonicGenerators(mesh, geom, loops);
        EXPECT_TRUE(
            basesAreDual(forms.primalGenerators, loops.primalGenerators));
        EXPECT_TRUE(basesAreDual(forms.dualGenerators, loops.dualGenerators));
        EXPECT_TRUE(primalFormsAreClosed(forms.primalGenerators));
        EXPECT_TRUE(primalFormsAreCoclosed(forms.primalGenerators, true));
        EXPECT_TRUE(dualFormsAreClosed(forms.dualGenerators, true));
        EXPECT_TRUE(dualFormsAreCoclosed(forms.dualGenerators));
        EXPECT_EQ(forms.dualGenerators.size(), dim);
        EXPECT_EQ(loops.dualGenerators.size(), dim);

        geom.unrequireDECOperators();
    }
}
