"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                let index = 0;
                for (let v of mesh.vertices) {
                        v.index = index++;
                }
                index = 0;
                for (let e of mesh.edges) {
                        e.index = index++;
                }
                index = 0;
                for (let f of mesh.faces) {
                        f.index = index++;
                }
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                let vn = mesh.vertices.length;
                let en = mesh.edges.length;

                let M = new Triplet(en, vn);
                for (let e of mesh.edges) {
                        let v1 = e.halfedge.vertex;
                        let v2 = e.halfedge.twin.vertex;
                        M.addEntry(1, e.index, v1.index);
                        M.addEntry(1, e.index, v2.index);
                }
                // for (let v of mesh.vertices) {
                //         for (let e of v.adjacentEdges()) {
                //                 M.addEntry(1, e.index, v.index);
                //         }
                // }

                return SparseMatrix.fromTriplet(M);
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                let en = mesh.edges.length;
                let fn = mesh.faces.length;

                let M = new Triplet(fn, en);
                for (let f of mesh.faces) {
                        for (let e of f.adjacentEdges()) {
                                M.addEntry(1, f.index, e.index);
                        }
                }

                return SparseMatrix.fromTriplet(M);
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                let vn = this.mesh.vertices.length;

                let vec = DenseMatrix.zeros(vn);
                for (let v of subset.vertices) {
                        vec.set(1, v);
                }

                return vec;
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                let en = this.mesh.edges.length;

                let vec = DenseMatrix.zeros(en);
                for (let e of subset.edges) {
                        vec.set(1, e);
                }

                return vec;
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                let fn = this.mesh.faces.length;

                let vec = DenseMatrix.zeros(fn);
                for (let f of subset.faces) {
                        vec.set(1, f);
                }

                return vec;
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                let s = MeshSubset.deepCopy(subset);

                let edgeVec = this.A0.timesDense(this.buildVertexVector(s));
                for (let i = 0; i < edgeVec.nRows(); i++) {
                        if (edgeVec.get(i)) {
                                s.addEdge(i);
                        }
                }

                let faceVec = this.A1.timesDense(this.buildEdgeVector(s));
                for (let i = 0; i < faceVec.nRows(); i++) {
                        if (faceVec.get(i)) {
                                s.addFace(i);
                        }
                }

                return s;
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                let s = MeshSubset.deepCopy(subset);

                let edgeVec = this.A1.transpose().timesDense(this.buildFaceVector(s));
                for (let i = 0; i < edgeVec.nRows(); i++) {
                        if (edgeVec.get(i)) {
                                s.addEdge(i);
                        }
                }

                let vertexVec = this.A0.transpose().timesDense(this.buildEdgeVector(s));
                for (let i = 0; i < vertexVec.nRows(); i++) {
                        if (vertexVec.get(i)) {
                                s.addVertex(i);
                        }
                }

                return s;
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                let s1 = this.closure(this.star(subset));
                let s2 = this.star(this.closure(subset));
                s1.deleteSubset(s2);

                return s1;
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                return subset.equals(this.closure(subset));
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                if (!this.isComplex(subset)) {
                        return -1;
                }

                let dim = (subset.faces.size > 0) ? 2 : (subset.edges.size > 0) ? 1 : 0;
                if (dim > 0) {
                        let vertexVec = this.buildVertexVector(subset);
                        let edgeVertexVec = this.A0.transpose().timesDense(this.buildEdgeVector(subset));
                        for (let i = 0; i < this.mesh.vertices.length; i++) {
                                if (vertexVec.get(i) && !edgeVertexVec.get(i)) {
                                        return -1;
                                }
                        }
                }
                if (dim > 1) {
                        let edgeVec = this.buildEdgeVector(subset);
                        let faceEdgeVec = this.A1.transpose().timesDense(this.buildFaceVector(subset));
                        for (let i = 0; i < this.mesh.edges.length; i++) {
                                if (edgeVec.get(i) && !faceEdgeVec.get(i)) {
                                        return -1;
                                }
                        }
                }

                return dim;
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                let dim = (subset.faces.size > 0) ? 2 : (subset.edges.size > 0) ? 1 : 0;

                let s = new MeshSubset();
                if (dim == 1) {
                        let edgeVertexVec = this.A0.transpose().timesDense(this.buildEdgeVector(subset));
                        for (let i = 0; i < this.mesh.vertices.length; i++) {
                                if (edgeVertexVec.get(i) == 1) {
                                        s.addVertex(i);
                                }
                        }
                } else if (dim == 2) {
                        let faceEdgeVec = this.A1.transpose().timesDense(this.buildFaceVector(subset));
                        for (let i = 0; i < this.mesh.edges.length; i++) {
                                if (faceEdgeVec.get(i) == 1) {
                                        s.addEdge(i);
                                }
                        }
                }
                s = this.closure(s);

                return s;
        }
}
