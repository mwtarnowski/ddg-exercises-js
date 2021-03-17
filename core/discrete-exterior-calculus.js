"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {
		let vn = geometry.mesh.vertices.length;

		let M = new Triplet(vn, vn);
		for (let v of geometry.mesh.vertices) {
			let i = vertexIndex[v];
			let dual = geometry.barycentricDualArea(v);
			M.addEntry(dual, i, i);
		}

		return SparseMatrix.fromTriplet(M);
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
		let en = geometry.mesh.edges.length;

		let M = new Triplet(en, en);
		for (let e of geometry.mesh.edges) {
			let i = edgeIndex[e];
			let dual = (geometry.cotan(e.halfedge) + geometry.cotan(e.halfedge.twin)) / 2;
			M.addEntry(dual, i, i);
		}

		return SparseMatrix.fromTriplet(M);
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
		let fn = geometry.mesh.faces.length;

		let M = new Triplet(fn, fn);
		for (let f of geometry.mesh.faces) {
			let i = faceIndex[f];
			let dual = 1 / geometry.area(f);
			M.addEntry(dual, i, i);
		}

		return SparseMatrix.fromTriplet(M);
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		let vn = geometry.mesh.vertices.length;
		let en = geometry.mesh.edges.length;

		let M = new Triplet(en, vn);
		for (let e of geometry.mesh.edges) {
			let i = edgeIndex[e];
			let j1 = vertexIndex[e.halfedge.vertex];  // target vertex (+1)
			M.addEntry( 1, i, j1);
			let j2 = vertexIndex[e.halfedge.twin.vertex];  // source vertex (-1)
			M.addEntry(-1, i, j2);
		}

		return SparseMatrix.fromTriplet(M);
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
		let en = geometry.mesh.edges.length;
		let fn = geometry.mesh.faces.length;

		let M = new Triplet(fn, en);
		for (let f of geometry.mesh.faces) {
			let i = faceIndex[f];
			for (let h of f.adjacentHalfedges()) {
				let j = edgeIndex[h.edge];
				if (h.edge.halfedge === h) {  // same orientation (+1)
					M.addEntry( 1, i, j);
				} else {  // opposite orientation (-1)
					M.addEntry(-1, i, j);
				}
			}
		}

		return SparseMatrix.fromTriplet(M);
	}
}
