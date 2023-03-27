#include "DelaunayMesh.hpp"
#include "../Types/Matrix.hpp"
#include <iostream>

using namespace gfd;

const double MINLENGTHSQ = 1e-8; // limit for minimal edge length
const double ORTOGONALSQ = 1e-13; // limit for linear dependency
const double RADIUSSCALE = 1.0000000001;

uint DelaunayMesh::searchNode(const Vector4 &p, uint curr) const {
	if(m_nsize == 0) return NONE;

	// check starter node
	if(curr >= m_nsize) curr = 0;

	// travel and find the nearest node
	double currsq = getRadiusSq(p, curr);
	while(true) {
		const uint prev = curr;
		const Buffer<uint> &e = getNodeEdges(prev);
		for(uint i=0; i<e.size(); i++) {
			const uint next = getEdgeOtherNode(e[i], prev);
			const double nextsq = getRadiusSq(p, next);
			if(nextsq < currsq) {
				curr = next;
				currsq = nextsq;
			}
		}
		if(prev == curr) break;
	}
	return curr;
}

uint DelaunayMesh::insertNode(const Vector4 &p, const double w, uint near, const bool forced) {
	// inserts a node into a convex delaunay mesh with constant and positive definite metric
	// p = new node position
	// w = new node weight
	// near = previous node where to start the search (can be used for optimization)
	// forced = (true -> create always new node) (false -> allow merging with previous nodes)
	// assumes that all flags are zero
	uint i, j, k, l;
	if(m_nsize == 0) { // current mesh is empty
		const uint node = addNode(p); // this is the first node
		setNodeWeight(node, w);
		return node;
	}

	// Check if the node already exists (or is too close to any other)
	near = searchNode(p, near);
	if(!forced) {
		if(w == 0.0 && m_w.empty()) { // weights are zero -> simple test
			if(getLengthSq(p, near) < MINLENGTHSQ) return near;
		}
		else { // a non-zero weight exist -> more complicated test
			if(getRadiusSq(p, near) < w + MINLENGTHSQ) return near;
			uint ns = 1;
			Buffer<uint> n(ns, near);
			for(i=0; i<ns; i++) {
				const uint ni = n[i];
				const Vector4 pi = getNodePosition(ni);
				const Vector4 di = getTransformed(p - pi, 0);
				if(di.dot(p - pi) + w < getNodeWeight(ni) + MINLENGTHSQ) return ni;
				const Buffer<uint> &ei = getNodeEdges(ni);
				for(j=0; j<ei.size(); j++) {
					const uint nj = getEdgeOtherNode(ei[j], ni);
					if(di.dot(getNodePosition(nj) - pi) > 0.0) n.gatherOnce(nj, ns);
				}
			}
		}
	}

	// Create new node
	const uint node = addNode(p);
	setNodeWeight(node, w);

	if(m_esize == 0) { // 0-dimensional convex mesh (only one node exists -> add the second one)
		addEdge(near, node);
		return node;
	}

	if(m_fsize == 0) { // 1-dimensional convex mesh (mesh on a straight line)
		uint curr = getNodeAnyEdge(near);
		if(getEdgeDeviation(curr, p).lensq() > ORTOGONALSQ) { // the new node increase the mesh dimension to 2
			increaseDimension();
			return node;
		}
		if(includesCell(p, curr)) { // insert node inside the mesh. replace curr with two edges
			CellSet detached;
			const Buffer<uint> en = getEdgeNodes(curr);
			detachEdge(curr, detached);
			for(j=0; j<en.size(); j++) attachEdge(en[j], node, detached);
		}
		else { // insert node outside of the mesh
			addEdge(curr, node);
		}
		return node;
	}

	if(m_bsize == 0) { // 2-dimensional convex mesh (mesh on a plane)
		uint curr = getNodeAnyFace(near);
		if(getFaceDeviation(curr, p).lensq() > ORTOGONALSQ) { // the new node increase the mesh dimension to 3
			increaseDimension();
			return node;
		}
		uint fes = 0;
		Buffer<uint> fe;
		CellSet detached;
		if(includesCell(p, curr)) { // insert node inside the mesh.
			fe = getFaceEdges(curr);
			fes = fe.size();
			detachFace(curr, detached);
		}
		else { // insert node outside of the mesh
			fe.gather(curr, fes);
			Buffer<uint> fn = getEdgeNodes(curr);
			uint fns = fn.size();
			while(fns > 0) {
				const Buffer<uint> &ne = getNodeEdges(fn[0]);
				for(i=0; i<ne.size(); i++) {
					const Buffer<uint> &ef = getEdgeFaces(ne[i]);
					if(ef.size() != 1) continue; // ne[i] is not on boundary
					if(fe.includes(ne[i], fes)) continue; // ne[i] is already on the list fe

					// the unique next boundary edge is found
					const Vector4 v = getEdgeDeviation(ne[i], p);
					const Buffer<uint> &en = getEdgeNodes(ne[i]);
					if(v.lensq() > ORTOGONALSQ && v.dot(getFaceDirection(ef[0], ne[i])) < 0.0) {
						fe.gather(ne[i], fes);
						for(j=0; j<en.size(); j++) fn.gatherOrUngather(en[j], fns);
						break;
					}
					for(j=0; j<en.size(); j++) fn.ungather(en[j], fns);
					break;
				}
			}
		}

		// remove recursively
		for(i=0; i<fes; i++) {
			const Buffer<uint> &ef = getEdgeFaces(fe[i]);
			if(getEdgeDeviation(fe[i], p).lensq() > ORTOGONALSQ) { // p is linearly independent of edge
				if(ef.empty()) continue; // fe[i] is a boundary edge
				const Vector4 fp = getFacePosition(ef[0]);
				const double fsq = getRadiusSq(fp, getEdgeNodes(fe[i]));
				if(!isInsideSphere(fp, fsq, node)) continue; // node is outside the face radius
			}

			// remove edge ->
			const Buffer<uint> efe = (ef.empty() ? Buffer<uint>(1, fe[i]) : getFaceEdges(ef[0]));
			for(j=0; j<efe.size(); j++) {
				if(fe.gatherOrUngather(efe[j], fes)) continue;
				if(i + 1 != 0) i--;
				detachEdgeRecursive(efe[j], detached);
			}
		}

		// insert faces
		for(i=0; i<fes; i++) {
			const Buffer<uint> en = getEdgeNodes(fe[i]);
			Buffer<uint> ee(en.size() + 1);
			for(j=0; j<en.size(); j++) ee[j] = attachEdge(en[j], node, detached);
			ee[j] = fe[i];
			fe[i] = attachFace(ee, detached);
		}
		for(i=0; i<fes; i++) mergeFace(fe[i], detached);

		// remove unnecessary edges
		removeDetached(detached);

		return node;
	}

	if(m_qsize == 0) { // 3-dimensional convex mesh
		uint curr = getNodeAnyBody(near);
		if(getBodyDeviation(curr, p).lensq() > ORTOGONALSQ) { // the new node increase the mesh dimension to 4
			increaseDimension();
			return node;
		}
		uint bfs = 0;
		Buffer<uint> bf;
		CellSet detached;
		if(includesCell(p, curr)) { // insert node inside the mesh.
			bf = getBodyFaces(curr);
			bfs = bf.size();
			detachBody(curr, detached);
		}
		else { // insert node outside of the mesh
			bf.gather(curr, bfs);
			Buffer<uint> be = getFaceEdges(curr);
			uint bes = be.size();
			while(bes > 0) {
				const Buffer<uint> &ef = getEdgeFaces(be[0]);
				for(i=0; i<ef.size(); i++) {
					const Buffer<uint> &fb = getFaceBodies(ef[i]);
					if(fb.size() != 1) continue; // ef[i] is not on boundary
					if(bf.includes(ef[i], bfs)) continue; // ef[i] is already on the list bf

					// the unique next boundary face is found
					const Vector4 v = getFaceDeviation(ef[i], p);
					const Buffer<uint> &fe = getFaceEdges(ef[i]);
					if(v.lensq() > ORTOGONALSQ && v.dot(getBodyDirection(fb[0], ef[i])) < 0.0) {
						bf.gather(ef[i], bfs);
						for(j=0; j<fe.size(); j++) be.gatherOrUngather(fe[j], bes);
						break;
					}
					for(j=0; j<fe.size(); j++) be.ungather(fe[j], bes);
					break;
				}
			}
		}

		// remove recursively
		for(i=0; i<bfs; i++) {
			const Buffer<uint> &fb = getFaceBodies(bf[i]);
			if(getFaceDeviation(bf[i], p).lensq() > ORTOGONALSQ) { // p is linearly independent of face
				if(fb.empty()) continue; // bf[i] is a boundary face
				const Vector4 bp = getBodyPosition(fb[0]);
				const double bsq = getRadiusSq(bp, getFaceNodes(bf[i]));
				if(!isInsideSphere(bp, bsq, node)) continue; // node is outside the body radius
			}

			// remove face ->
			const Buffer<uint> fbf = (fb.empty() ? Buffer<uint>(1, bf[i]) : getBodyFaces(fb[0]));
			for(j=0; j<fbf.size(); j++) {
				if(bf.gatherOrUngather(fbf[j], bfs)) continue;
				if(i + 1 != 0) i--;
				detachFaceRecursive(fbf[j], detached);
			}
		}

		// insert bodies
		for(i=0; i<bfs; i++) {
			const Buffer<uint> fe = getFaceEdges(bf[i]);
			Buffer<uint> ff(fe.size() + 1);
			for(j=0; j<fe.size(); j++) {
				const Buffer<uint> en = getEdgeNodes(fe[j]);
				Buffer<uint> ee(en.size() + 1);
				for(k=0; k<en.size(); k++) ee[k] = attachEdge(en[k], node, detached);
				ee[k] = fe[j];
				ff[j] = attachFace(ee, detached);
			}
			ff[j] = bf[i];
			bf[i] = attachBody(ff, detached);
		}
		for(i=0; i<bfs; i++) mergeBody(bf[i], detached);

		// remove unnecessary faces and edges
		removeDetached(detached);
		return node;
	}

	// 4-dimensional mesh
	uint curr = getNodeAnyQuad(near);
	uint qbs = 0;
	Buffer<uint> qb;
	CellSet detached;
	if(includesCell(p, curr)) { // insert node inside the mesh.
		qb = getQuadBodies(curr);
		qbs = qb.size();
		removeQuad(curr);
	}
	else { // insert node outside of the mesh
		qb.gather(curr, qbs);
		Buffer<uint> qf = getBodyFaces(curr);
		uint qfs = qf.size();
		while(qfs > 0) {
			const Buffer<uint> &fb = getFaceBodies(qf[0]);
			for(i=0; i<fb.size(); i++) {
				const Buffer<uint> &bq = getBodyQuads(fb[i]);
				if(bq.size() != 1) continue; // fb[i] is not on boundary
				if(qb.includes(fb[i], qbs)) continue; // fb[i] is already on the list qb

				// the unique next boundary body is found
				const Vector4 v = getBodyDeviation(fb[i], p);
				const Buffer<uint> &bf = getBodyFaces(fb[i]);
				if(v.lensq() > ORTOGONALSQ && v.dot(getQuadDirection(bq[0], fb[i])) < 0.0) {
					qb.gather(fb[i], qbs);
					for(j=0; j<bf.size(); j++) qf.gatherOrUngather(bf[j], qfs);
					break;
				}
				for(j=0; j<bf.size(); j++) qf.ungather(bf[j], qfs);
				break;
			}
		}
	}

	// remove recursively
	for(i=0; i<qbs; i++) {
		const Buffer<uint> &bq = getBodyQuads(qb[i]);
		if(getBodyDeviation(qb[i], p).lensq() > ORTOGONALSQ) { // p is linearly independent of body
			if(bq.empty()) continue; // qb[i] is a boundary body
			const Vector4 qp = getQuadPosition(bq[0]);
			const double qsq = getRadiusSq(qp, getBodyNodes(qb[i]));
			if(!isInsideSphere(qp, qsq, node)) continue; // node is outside the quad radius
		}

		// remove body ->
		const Buffer<uint> bqb = (bq.empty() ? Buffer<uint>(1, qb[i]) : getQuadBodies(bq[0]));
		for(j=0; j<bqb.size(); j++) {
			if(qb.gatherOrUngather(bqb[j], qbs)) continue;
			if(i + 1 != 0) i--;
			detachBodyRecursive(bqb[j], detached);
		}
	}

	// insert quads
	for(i=0; i<qbs; i++) {
		const Buffer<uint> bf = getBodyFaces(qb[i]);
		Buffer<uint> bb(bf.size() + 1);
		for(j=0; j<bf.size(); j++) {
			const Buffer<uint> fe = getFaceEdges(bf[j]);
			Buffer<uint> ff(fe.size() + 1);
			for(k=0; k<fe.size(); k++) {
				const Buffer<uint> en = getEdgeNodes(fe[k]);
				Buffer<uint> ee(en.size() + 1);
				for(l=0; l<en.size(); l++) ee[l] = attachEdge(en[l], node, detached);
				ee[l] = fe[k];
				ff[k] = attachFace(ee, detached);
			}
			ff[k] = bf[j];
			bb[j] = attachBody(ff, detached);
		}
		bb[j] = qb[i];
		qb[i] = attachQuad(bb, detached);
	}
	for(i=0; i<qbs; i++) mergeQuad(qb[i], detached);

	removeDetached(detached);
	return node;
}

bool DelaunayMesh::eraseNode(const uint n) {
	uint i, j, k;
	if(n >= m_nsize) return false;

	if(m_esize == 0) { // 0-dimensional convex mesh (only one node exists)
		removeNode(n);
		return true;
	}

	if(m_fsize == 0) { // 1-dimensional convex mesh (mesh on a straight line)
		const Buffer<uint> &ne = getNodeEdges(n);
		if(ne.size() >= 2) addEdge(getEdgeOtherNode(ne[0], n), getEdgeOtherNode(ne[1], n));
		removeNode(n);
		return true;
	}

	CellSet detached;
	const Vector4 np = getNodePosition(n);
	if(m_bsize == 0) {
		// remove node and store hole boundary
		uint es = 0;
		Buffer<uint> e;
		const Buffer<uint> nf = getNodeFaces(n);
		for(i=0; i<nf.size(); i++) {
			const Buffer<uint> &fe = getFaceEdges(nf[i]);
			for(j=0; j<fe.size(); j++) e.gatherOnce(fe[j], es);
		}
		const Buffer<uint> &ne = getNodeEdges(n);
		for(i=0; i<ne.size(); i++) e.ungather(ne[i], es);
		detachNode(n, detached);

		// store boundary nodes
		uint nns = 0;
		Buffer<uint> nn;
		for(i=0; i<es; i++) {
			const Buffer<uint> &en = getEdgeNodes(e[i]);
			for(j=0; j<en.size(); j++) nn.gatherOnce(en[j], nns);
		}

		// create new mesh
		uint node = 0;
		DelaunayMesh mesh(m_dim);
		mesh.setMetric(getMetric());
		for(i=0; i<nns; i++) node = mesh.insertNode(getNodePosition(nn[i]), getNodeWeight(nn[i]), node, true);

		// remove elements outside
		uint fs = 0;
		Buffer<uint> f;
		Buffer<bool> estay(mesh.getEdgeSize(), false);
		Buffer<bool> fstay(mesh.getFaceSize(), false);
		for(i=0; i<es; i++) {
			// find identical edge from the mesh
			Buffer<uint> en = getEdgeNodes(e[i]);
			for(j=0; j<en.size(); j++) {
				for(k=0; nn[k] != en[j]; k++);
				en[j] = k;
			}
			const Buffer<uint> &ne = mesh.getNodeEdges(en[0]);
			for(j=0; !en.isAnagram(mesh.getEdgeNodes(ne[j])); j++);
			estay[ne[j]] = true;

			// gather inside faces into f
			const Vector4 v = mesh.getEdgeDeviation(ne[j], np);
			const Buffer<uint> &ef = mesh.getEdgeFaces(ne[j]);
			for(k=0; k<ef.size(); k++) {
				if(v.dot(mesh.getFaceDirection(ef[k], ne[j])) > 0.0) f.gatherOnce(ef[k], fs);
			}
		}
		for(i=0; i<fs; i++) { // try to spread out removable area
			fstay[f[i]] = true;
			const Buffer<uint> &fe = mesh.getFaceEdges(f[i]);
			for(j=0; j<fe.size(); j++) {
				if(estay[fe[j]]) continue;
				estay[fe[j]] = true;

				const Buffer<uint> &ef = mesh.getEdgeFaces(fe[j]);
				for(k=0; k<ef.size(); k++) f.gatherOnce(ef[k], fs);
			}
		}
		for(i=fstay.size(); i-->0; ) if(!fstay[i]) mesh.removeFace(i);
		for(i=estay.size(); i-->0; ) if(!estay[i]) mesh.removeEdge(i);

		// copy elements from mesh
		Buffer<uint> ee(mesh.getEdgeSize());
		for(i=0; i<ee.size(); i++) {
			const Buffer<uint> &en = mesh.getEdgeNodes(i);
			ee[i] = attachEdge(nn[en[0]], nn[en[1]], detached);
		}
		for(i=0; i<mesh.getFaceSize(); i++) {
			const Buffer<uint> fe = mesh.getFaceEdges(i);
			for(j=0; j<fe.size(); j++) fe[j] = ee[fe[j]];
			attachFace(fe, detached);
		}
		removeDetached(detached);
		return true;
	}

	if(m_qsize == 0) {
		// remove node and store hole boundary
		uint fs = 0;
		Buffer<uint> f;
		const Buffer<uint> nb = getNodeBodies(n);
		for(i=0; i<nb.size(); i++) {
			const Buffer<uint> &bf = getBodyFaces(nb[i]);
			for(j=0; j<bf.size(); j++) f.gatherOnce(bf[j], fs);
		}
		const Buffer<uint> nf = getNodeFaces(n);
		for(i=0; i<nf.size(); i++) f.ungather(nf[i], fs);
		detachNode(n, detached);

		// store boundary nodes
		uint nns = 0;
		Buffer<uint> nn;
		for(i=0; i<fs; i++) {
			const Buffer<uint> fn = getFaceNodes(f[i]);
			for(j=0; j<fn.size(); j++) nn.gatherOnce(fn[j], nns);
		}

		// create new mesh
		uint node = 0;
		DelaunayMesh mesh(m_dim);
		mesh.setMetric(getMetric());
		for(i=0; i<nns; i++) node = mesh.insertNode(getNodePosition(nn[i]), getNodeWeight(nn[i]), node, true);

		// remove elements outside
		uint bs = 0;
		Buffer<uint> b;
		Buffer<bool> estay(mesh.getEdgeSize(), false);
		Buffer<bool> fstay(mesh.getFaceSize(), false);
		Buffer<bool> bstay(mesh.getBodySize(), false);
		for(i=0; i<fs; i++) {
			// find identical face from the mesh
			Buffer<uint> fn = getFaceNodes(f[i]);
			for(j=0; j<fn.size(); j++) {
				for(k=0; nn[k] != fn[j]; k++);
				fn[j] = k;
			}
			const Buffer<uint> nf = mesh.getNodeFaces(fn[0]);
			for(j=0; !fn.isAnagram(mesh.getFaceNodes(nf[j])); j++);
			fstay[nf[j]] = true;

			// gather inside bodies into b
			const Vector4 v = mesh.getFaceDeviation(nf[j], np);
			const Buffer<uint> &fb = mesh.getFaceBodies(nf[j]);
			for(k=0; k<fb.size(); k++) {
				if(v.dot(mesh.getBodyDirection(fb[k], nf[j])) > 0.0) b.gatherOnce(fb[k], bs);
			}
		}
		for(i=0; i<bs; i++) { // try to spread out removable area
			bstay[b[i]] = true;
			const Buffer<uint> &bf = mesh.getBodyFaces(b[i]);
			for(j=0; j<bf.size(); j++) {
				const Buffer<uint> &fe = mesh.getFaceEdges(bf[j]);
				for(k=0; k<fe.size(); k++) estay[fe[k]] = true;

				if(fstay[bf[j]]) continue;
				fstay[bf[j]] = true;

				const Buffer<uint> &fb = mesh.getFaceBodies(bf[j]);
				for(k=0; k<fb.size(); k++) b.gatherOnce(fb[k], bs);
			}
		}
		for(i=bstay.size(); i-->0; ) if(!bstay[i]) mesh.removeBody(i);
		for(i=fstay.size(); i-->0; ) if(!fstay[i]) mesh.removeFace(i);
		for(i=estay.size(); i-->0; ) if(!estay[i]) mesh.removeEdge(i);

		// copy elements from the mesh
		Buffer<uint> ee(mesh.getEdgeSize());
		for(i=0; i<ee.size(); i++) {
			const Buffer<uint> &en = mesh.getEdgeNodes(i);
			ee[i] = attachEdge(nn[en[0]], nn[en[1]], detached);
		}
		Buffer<uint> ff(mesh.getFaceSize());
		for(i=0; i<ff.size(); i++) {
			const Buffer<uint> fe = mesh.getFaceEdges(i);
			for(j=0; j<fe.size(); j++) fe[j] = ee[fe[j]];
			ff[i] = attachFace(fe, detached);
		}
		for(i=0; i<mesh.getBodySize(); i++) {
			const Buffer<uint> bf = mesh.getBodyFaces(i);
			for(j=0; j<bf.size(); j++) bf[j] = ff[bf[j]];
			attachBody(bf, detached);
		}
		removeDetached(detached);
		return true;
	}


	// remove node and store hole boundary
	uint bs = 0;
	Buffer<uint> b;
	const Buffer<uint> nq = getNodeQuads(n);
	for(i=0; i<nq.size(); i++) {
		const Buffer<uint> &qb = getQuadBodies(nq[i]);
		for(j=0; j<qb.size(); j++) b.gatherOnce(qb[j], bs);
	}
	const Buffer<uint> nb = getNodeBodies(n);
	for(i=0; i<nb.size(); i++) b.ungather(nb[i], bs);
	detachNode(n, detached);

	// store boundary nodes
	uint nns = 0;
	Buffer<uint> nn;
	for(i=0; i<bs; i++) {
		const Buffer<uint> bn = getBodyNodes(b[i]);
		for(j=0; j<bn.size(); j++) nn.gatherOnce(bn[j], nns);
	}

	// create new mesh
	uint node = 0;
	DelaunayMesh mesh(m_dim);
	mesh.setMetric(getMetric());
	for(i=0; i<nns; i++) node = mesh.insertNode(getNodePosition(nn[i]), getNodeWeight(nn[i]), node, true);

	// remove elements outside
	uint qs = 0;
	Buffer<uint> q;
	Buffer<bool> estay(mesh.getEdgeSize(), false);
	Buffer<bool> fstay(mesh.getFaceSize(), false);
	Buffer<bool> bstay(mesh.getBodySize(), false);
	Buffer<bool> qstay(mesh.getQuadSize(), false);
	for(i=0; i<bs; i++) {
		// find identical body from the mesh
		Buffer<uint> bn = getBodyNodes(b[i]);
		for(j=0; j<bn.size(); j++) {
			for(k=0; nn[k] != bn[j]; k++);
			bn[j] = k;
		}
		const Buffer<uint> nb = mesh.getNodeBodies(bn[0]);
		for(j=0; !bn.isAnagram(mesh.getBodyNodes(nb[j])); j++);
		bstay[nb[j]] = true;

		// gather inside quads into q
		const Vector4 v = mesh.getBodyDeviation(nb[j], np);
		const Buffer<uint> &bq = mesh.getBodyQuads(nb[j]);
		for(k=0; k<bq.size(); k++) {
			if(v.dot(mesh.getQuadDirection(bq[k], nb[j])) > 0.0) q.gatherOnce(bq[k], qs);
		}
	}
	for(i=0; i<qs; i++) { // try to spread out removable area
		qstay[q[i]] = true;
		const Buffer<uint> &qb = mesh.getQuadBodies(q[i]);
		for(j=0; j<qb.size(); j++) {
			const Buffer<uint> &bf = mesh.getBodyFaces(qb[j]);
			for(k=0; k<bf.size(); k++) fstay[bf[k]] = true;
			const Buffer<uint> be = mesh.getBodyEdges(qb[j]);
			for(k=0; k<be.size(); k++) estay[be[k]] = true;

			if(bstay[qb[j]]) continue;
			bstay[qb[j]] = true;

			const Buffer<uint> &bq = mesh.getBodyQuads(qb[j]);
			for(k=0; k<bq.size(); k++) q.gatherOnce(bq[k], qs);
		}
	}
	for(i=qstay.size(); i-->0; ) if(!qstay[i]) mesh.removeQuad(i);
	for(i=bstay.size(); i-->0; ) if(!bstay[i]) mesh.removeBody(i);
	for(i=fstay.size(); i-->0; ) if(!fstay[i]) mesh.removeFace(i);
	for(i=estay.size(); i-->0; ) if(!estay[i]) mesh.removeEdge(i);

	// copy elements from mesh
	Buffer<uint> ee(mesh.getEdgeSize());
	for(i=0; i<ee.size(); i++) {
		const Buffer<uint> &en = mesh.getEdgeNodes(i);
		ee[i] = attachEdge(nn[en[0]], nn[en[1]], detached);
	}
	Buffer<uint> ff(mesh.getFaceSize());
	for(i=0; i<ff.size(); i++) {
		const Buffer<uint> fe = mesh.getFaceEdges(i);
		for(j=0; j<fe.size(); j++) fe[j] = ee[fe[j]];
		ff[i] = attachFace(fe, detached);
	}
	Buffer<uint> bb(mesh.getBodySize());
	for(i=0; i<bb.size(); i++) {
		const Buffer<uint> bf = mesh.getBodyFaces(i);
		for(j=0; j<bf.size(); j++) bf[j] = ff[bf[j]];
		bb[i] = attachBody(bf, detached);
	}
	for(i=0; i<mesh.getQuadSize(); i++) {
		const Buffer<uint> qb = mesh.getQuadBodies(i);
		for(j=0; j<qb.size(); j++) qb[j] = bb[qb[j]];
		attachQuad(qb, detached);
	}
	removeDetached(detached);
	return true;
}

bool DelaunayMesh::insertMesh(const DelaunayMesh &mesh) {
	if(getDimension() != mesh.getDimension()) return false; // dimensions must match

	// erase nodes that prevent insertion
	Buffer<uint> nn(mesh.getNodeSize(), NONE);
	Buffer<uint> ee(mesh.getEdgeSize(), NONE);
	Buffer<uint> ff(mesh.getFaceSize(), NONE);
	Buffer<uint> bb(mesh.getBodySize(), NONE);
	uint i, j, k, near;
	Vector4 p;
	if(m_dim == 4)
	{
		near = 0;
		for(i=0; i<mesh.getBodySize(); i++)
		{
			const Buffer<uint> &c = mesh.getBodyQuads(i);
			if(c.size() >= 2) continue; // ignore inner cells to consider only boundaries

			// mark boundary nodes
			const Buffer<uint> n = mesh.getBodyNodes(i);
			for(j=0; j<n.size(); j++) nn[n[j]] = 0;
			const Buffer<uint> e = mesh.getBodyEdges(i);
			for(j=0; j<e.size(); j++) ee[e[j]] = 0;
			const Buffer<uint> &f = mesh.getBodyFaces(i);
			for(j=0; j<f.size(); j++) ff[f[j]] = 0;
			bb[i] = 0;

			// erase invalid nodes
			if(c.empty()) p = mesh.getBodyPosition(i);
			else p = mesh.getQuadPosition(c[0]);
			near = eraseNodesInside(p, mesh.getRadiusSq(p, n), near);
		}
	}
	if(m_dim >= 3)
	{
		near = 0;
		for(i=0; i<mesh.getFaceSize(); i++)
		{
			const Buffer<uint> &c = mesh.getFaceBodies(i);
			if(c.size() >= 2) continue; // ignore inner cells to consider only boundaries
			if(m_dim > 3 && !c.empty()) continue; // this condition is applied only for optimization

			// mark boundary nodes
			const Buffer<uint> n = mesh.getFaceNodes(i);
			for(j=0; j<n.size(); j++) nn[n[j]] = 0;
			const Buffer<uint> &e = mesh.getFaceEdges(i);
			for(j=0; j<e.size(); j++) ee[e[j]] = 0;
			ff[i] = 0;

			// erase invalid nodes
			if(c.empty()) p = mesh.getFacePosition(i);
			else p = mesh.getBodyPosition(c[0]);
			near = eraseNodesInside(p, mesh.getRadiusSq(p, n), near);
		}
	}
	if(m_dim >= 2)
	{
		near = 0;
		for(i=0; i<mesh.getEdgeSize(); i++)
		{
			const Buffer<uint> &c = mesh.getEdgeFaces(i);
			if(c.size() >= 2) continue; // ignore inner cells to consider only boundaries
			if(m_dim > 2 && !c.empty()) continue; // this condition is applied only for optimization

			// mark boundary nodes
			const Buffer<uint> &n = mesh.getEdgeNodes(i);
			for(j=0; j<n.size(); j++) nn[n[j]] = 0;
			ee[i] = 0;

			// erase invalid nodes
			if(c.empty()) p = mesh.getEdgePosition(i);
			else p = mesh.getFacePosition(c[0]);
			near = eraseNodesInside(p, mesh.getRadiusSq(p, n), near);
		}
	}
	near = 0;
	for(i=0; i<mesh.getNodeSize(); i++)
	{
		const Buffer<uint> &c = mesh.getNodeEdges(i);
		if(c.size() >= 2) continue; // ignore inner cells to consider only boundaries
		if(m_dim > 1 && !c.empty()) continue; // this condition is applied only for optimization

		// mark boundary nodes
		nn[i] = 0;

		// erase invalid nodes
		if(c.empty()) p = mesh.getNodePosition(i);
		else p = mesh.getEdgePosition(c[0]);
		near = eraseNodesInside(p, mesh.getRadiusSq(p, i), near);
	}

	// insert boundary cells
	for(i=0; i<nn.size(); i++)
	{
		if(nn[i] == NONE) continue;
		nn[i] = insertNode(mesh.getNodePosition(i), getNodeWeight(i), near, true);
		setNodeFlag(nn[i], mesh.getNodeFlag(i));
		near = nn[i];
	}
	for(i=0; i<ee.size(); i++)
	{
		if(ee[i] == NONE) continue;
		const Buffer<uint> &n = mesh.getEdgeNodes(i);
		ee[i] = addEdge(nn[n[0]], nn[n[1]]);
		setEdgeFlag(ee[i], mesh.getEdgeFlag(i));
	}
	for(i=0; i<ff.size(); i++)
	{
		if(ff[i] == NONE) continue;
		Buffer<uint> e = mesh.getFaceEdges(i);
		for(j=0; j<e.size(); j++) e[j] = ee[e[j]];
		ff[i] = addFace(e);
		setFaceFlag(ff[i], mesh.getFaceFlag(i));
	}
	for(i=0; i<bb.size(); i++)
	{
		if(bb[i] == NONE) continue;
		Buffer<uint> f = mesh.getBodyFaces(i);
		for(j=0; j<f.size(); j++) f[j] = ff[f[j]];
		bb[i] = addBody(f);
		setBodyFlag(bb[i], mesh.getBodyFlag(i));
	}

	// remove overlapping cells
	CellSet detached;
	uint rems = 0;
	Buffer<uint> rem;
	switch(m_dim) {
	case 1: {
		Buffer<char> boun(getNodeSize(), 0);
		for(i=0; i<nn.size(); i++)
		{
			if(nn[i] == NONE) continue;
			const Buffer<uint> &c = mesh.getNodeEdges(i);
			if(c.size() != 1) continue; // consider boundary cells with one parent cell

			const Vector4 v = mesh.getEdgeAverage(c[0]) - mesh.getNodePosition(i);
			const Buffer<uint> &next = getNodeEdges(nn[i]);
			for(j=0; j<next.size(); j++)
			{
				if(v.dot(getEdgeAverage(next[j]) - getNodePosition(nn[i])) < 0.0) continue;
				rem.gatherOnce(next[j], rems);
				boun[nn[i]] = 1;
				break;
			}
		}
		// spread removable area
		for(i=0; i<rems; i++)
		{
			const Buffer<uint> &c = getEdgeNodes(rem[i]);
			for(j=0; j<c.size(); j++)
			{
				if(boun[c[j]] != 0) continue;
				boun[c[j]] = -1;
				const Buffer<uint> &cc = getNodeEdges(c[j]);
				for(k=0; k<cc.size(); k++) rem.gatherOnce(cc[k], rems);
			}
			detachEdge(rem[i], detached);
		}
		for(i=0; i<boun.size(); i++) { if(boun[i] < 0) detachNode(i, detached); } // remove recursively interior
		break;
	}
	case 2: {
		Buffer<char> boun(getEdgeSize(), 0);
		for(i=0; i<ee.size(); i++)
		{
			if(ee[i] == NONE) continue;
			const Buffer<uint> &c = mesh.getEdgeFaces(i);
			if(c.size() != 1) continue; // consider boundary cells with one parent cell

			const Vector4 v = mesh.getEdgeDeviation(i, mesh.getNodeAverage(mesh.getFaceNodes(c[0])));
			const Buffer<uint> &next = getEdgeFaces(ee[i]);
			for(j=0; j<next.size(); j++)
			{
				if(v.dot(getNodeAverage(getFaceNodes(next[j])) - getNodePosition(getEdgeAnyNode(ee[i]))) < 0.0) continue;
				rem.gatherOnce(next[j], rems);
				boun[ee[i]] = 1;
				break;
			}
		}
		// spread removable area
		for(i=0; i<rems; i++)
		{
			const Buffer<uint> &c = getFaceEdges(rem[i]);
			for(j=0; j<c.size(); j++)
			{
				if(boun[c[j]] != 0) continue;
				boun[c[j]] = -1;
				const Buffer<uint> &cc = getEdgeFaces(c[j]);
				for(k=0; k<cc.size(); k++) rem.gatherOnce(cc[k], rems);
			}
			detachFace(rem[i], detached);
		}
		for(i=0; i<boun.size(); i++) { if(boun[i] < 0) detachEdgeRecursive(i, detached); } // remove recursively interior
		break;
	}
	case 3: {
		Buffer<char> boun(getFaceSize(), 0);
		for(i=0; i<ff.size(); i++)
		{
			if(ff[i] == NONE) continue;
			const Buffer<uint> &c = mesh.getFaceBodies(i);
			if(c.size() != 1) continue; // consider boundary cells with one parent cell

			const Vector4 v = mesh.getFaceDeviation(i, mesh.getNodeAverage(mesh.getBodyNodes(c[0])));
			const Buffer<uint> &next = getFaceBodies(ff[i]);
			for(j=0; j<next.size(); j++)
			{
				if(v.dot(getNodeAverage(getBodyNodes(next[j])) - getNodePosition(getFaceAnyNode(ff[i]))) < 0.0) continue;
				rem.gatherOnce(next[j], rems);
				boun[ff[i]] = 1;
				break;
			}
		}
		for(i=0; i<rems; i++)
		{
			const Buffer<uint> &c = getBodyFaces(rem[i]);
			for(j=0; j<c.size(); j++)
			{
				if(boun[c[j]] != 0) continue;
				boun[c[j]] = -1;
				const Buffer<uint> &cc = getFaceBodies(c[j]);
				for(k=0; k<cc.size(); k++) rem.gatherOnce(cc[k], rems);
			}
			detachBody(rem[i], detached);
		}
		for(i=0; i<boun.size(); i++) { if(boun[i] < 0) detachFaceRecursive(i, detached); } // remove recursively interior
		break;
	}
	default: {
		Buffer<char> boun(getBodySize(), 0);
		for(i=0; i<bb.size(); i++)
		{
			if(bb[i] == NONE) continue;
			const Buffer<uint> &c = mesh.getBodyQuads(i);
			if(c.size() != 1) continue; // consider boundary cells with one parent cell

			const Vector4 v = mesh.getBodyDeviation(i, mesh.getNodeAverage(mesh.getQuadNodes(c[0])));
			const Buffer<uint> &next = getBodyQuads(bb[i]);
			for(j=0; j<next.size(); j++)
			{
				if(v.dot(getNodeAverage(getQuadNodes(next[j])) - getNodePosition(getBodyAnyNode(bb[i]))) < 0.0) continue;
				rem.gatherOnce(next[j], rems);
				boun[bb[i]] = 1;
				break;
			}
		}
		// spread removable area
		for(i=0; i<rems; i++)
		{
			const Buffer<uint> &c = getQuadBodies(rem[i]);
			for(j=0; j<c.size(); j++)
			{
				if(boun[c[j]] != 0) continue;
				boun[c[j]] = -1;
				const Buffer<uint> &cc = getBodyQuads(c[j]);
				for(k=0; k<cc.size(); k++) rem.gatherOnce(cc[k], rems);
			}
			detachQuad(rem[i], detached);
		}
		for(i=0; i<boun.size(); i++) { if(boun[i] < 0) detachBodyRecursive(i, detached); } // remove recursively interior
		break;
	}
	}

	// add missing cells
	for(i=0; i<nn.size(); i++)
	{
		if(nn[i] != NONE) continue;
		nn[i] = attachNode(mesh.getNodePosition(i), detached);
		setNodeWeight(nn[i], mesh.getNodeWeight(i));
		setNodeFlag(nn[i], mesh.getNodeFlag(i));
	}
	for(i=0; i<ee.size(); i++)
	{
		if(ee[i] != NONE) continue;
		const Buffer<uint> &n = mesh.getEdgeNodes(i);
		ee[i] = attachEdge(nn[n[0]], nn[n[1]], detached);
		setEdgeFlag(ee[i], mesh.getEdgeFlag(i));
	}
	for(i=0; i<ff.size(); i++)
	{
		if(ff[i] != NONE) continue;
		Buffer<uint> e = mesh.getFaceEdges(i);
		for(j=0; j<e.size(); j++) e[j] = ee[e[j]];
		ff[i] = attachFace(e, detached);
		setFaceFlag(ff[i], mesh.getFaceFlag(i));
	}
	for(i=0; i<bb.size(); i++)
	{
		if(bb[i] != NONE) continue;
		Buffer<uint> f = mesh.getBodyFaces(i);
		for(j=0; j<f.size(); j++) f[j] = ff[f[j]];
		bb[i] = attachBody(f, detached);
		setBodyFlag(bb[i], mesh.getBodyFlag(i));
	}
	for(i=0; i<mesh.getQuadSize(); i++)
	{
		Buffer<uint> b = mesh.getQuadBodies(i);
		for(j=0; j<b.size(); j++) b[j] = bb[b[j]];
		const uint qq = attachQuad(b, detached);
		setQuadFlag(qq, mesh.getQuadFlag(i));
	}

	removeDetached(detached);
	return true;
}

uint DelaunayMesh::attachNode(const Vector4 &p, CellSet &detached) {
	// use detached slot if possible
	uint res;
	if(detached.ns > 0) res = detached.n[--detached.ns];
	else {
		// check if m_n is full -> resize the table
		res = m_nsize++;
		if(m_nsize > m_n.size()) resizeNodeBuffer(2 * m_nsize);
	}

	// create new node
	setNodePosition(res, p);
	return res;
}

uint DelaunayMesh::attachEdge(const uint n0, const uint n1, CellSet &detached) {
	// check if edge already exists
	uint res = findEdge(n0, n1);
	if(res != NONE) return res;

	// use detached slot if possible
	if(detached.es > 0) res = detached.e[--detached.es];
	else {
		// check if m_e is full -> resize the table
		res = m_esize++;
		if(m_esize > m_e.size()) resizeEdgeBuffer(2 * m_esize);
	}

	// create new edge
	m_e[res].n.resize(2);
	m_e[res].n[0] = n0;
	m_e[res].n[1] = n1;
	m_n[n0].e.push_back(res);
	m_n[n1].e.push_back(res);
	return res;
}

uint DelaunayMesh::attachFace(const Buffer<uint> &e, CellSet &detached) {
	// check if face already exists
	uint res = findFace(e);
	if(res != NONE) return res;

	// use detached slot if possible
	if(detached.fs > 0) res = detached.f[--detached.fs];
	else {
		// check if m_f is full -> resize the table
		res = m_fsize++;
		if(m_fsize > m_f.size()) resizeFaceBuffer(2 * m_fsize);
	}

	// create new face
	m_f[res].e = e;
	for(uint i=0; i<e.size(); i++) m_e[e[i]].f.push_back(res);
	orderFaceEdges(res);
	return res;
}

uint DelaunayMesh::attachBody(const Buffer<uint> &f, CellSet &detached) {
	// check if body already exists
	uint res = findBody(f);
	if(res != NONE) return res;

	// use detached slot if possible
	if(detached.bs > 0) res = detached.b[--detached.bs];
	else {
		// check if m_b is full -> resize the table
		res = m_bsize++;
		if(m_bsize > m_b.size()) resizeBodyBuffer(2 * m_bsize);
	}

	// create new body
	m_b[res].f = f;
	for(uint i=0; i<f.size(); i++) m_f[f[i]].b.push_back(res);
	orderBodyFaces(res);
	return res;
}

uint DelaunayMesh::attachQuad(const Buffer<uint> &b, CellSet &detached) {
	// check if quad already exists
	uint res = findQuad(b);
	if(res != NONE) return res;

	// use detached slot if possible
	if(detached.qs > 0) res = detached.q[--detached.qs];
	else {
		// check if m_q is full -> resize the table
		res = m_qsize++;
		if(m_qsize > m_q.size()) resizeQuadBuffer(2 * m_qsize);
	}

	// create new quad
	m_q[res].b = b;
	for(uint i=0; i<b.size(); i++) m_b[b[i]].q.push_back(res);
	orderQuadBodies(res);
	return res;
}

void DelaunayMesh::detachNode(const uint n, CellSet &detached) {
	detached.n.gather(n, detached.ns);
	setNodeFlag(n, 0);
	setNodeWeight(n, 0.0);

	// detach linked edges
	Buffer<uint> &e = m_n[n].e;
	for(uint i=e.size(); i-->0; ) detachEdge(e[i], detached);
}

void DelaunayMesh::detachEdge(const uint e, CellSet &detached) {
	uint i;
	detached.e.gather(e, detached.es);
	setEdgeFlag(e, 0);

	// detach linked faces
	Buffer<uint> &f = m_e[e].f;
	for(i=f.size(); i-->0; ) detachFace(f[i], detached);

	// remove links from nodes
	Buffer<uint> &n = m_e[e].n;
	for(i=n.size(); i-->0; ) m_n[n[i]].e.eraseFirst(e);
	n.clear();
}

void DelaunayMesh::detachFace(const uint f, CellSet &detached) {
	uint i;
	detached.f.gather(f, detached.fs);
	setFaceFlag(f, 0);

	// detach linked bodies
	Buffer<uint> &b = m_f[f].b;
	for(i=b.size(); i-->0; ) detachBody(b[i], detached);

	// remove links from edges
	Buffer<uint> &e = m_f[f].e;
	for(i=e.size(); i-->0; ) m_e[e[i]].f.eraseFirst(f);
	e.clear();
}

void DelaunayMesh::detachBody(const uint b, CellSet &detached) {
	uint i;
	detached.b.gather(b, detached.bs);
	setBodyFlag(b, 0);

	// detach linked quads
	Buffer<uint> &q = m_b[b].q;
	for(i=q.size(); i-->0; ) detachQuad(q[i], detached);

	// remove links from faces
	Buffer<uint> &f = m_b[b].f;
	for(i=f.size(); i-->0; ) m_f[f[i]].b.eraseFirst(b);
	f.clear();
}

void DelaunayMesh::detachQuad(const uint q, CellSet &detached) {
	detached.q.gather(q, detached.qs);
	setQuadFlag(q, 0);

	// remove links from bodies
	Buffer<uint> &b = m_q[q].b;
	for(uint i=b.size(); i-->0; ) m_b[b[i]].q.eraseFirst(q);
	b.clear();
}

void DelaunayMesh::detachEdgeRecursive(const uint e, CellSet &detached) {
	const Buffer<uint> n = getEdgeNodes(e);
	detachEdge(e, detached);
	for(uint i=0; i<n.size(); i++) {
		if(getNodeEdges(n[i]).empty()) detachNode(n[i], detached);
	}
}

void DelaunayMesh::detachFaceRecursive(const uint f, CellSet &detached) {
	const Buffer<uint> e = getFaceEdges(f);
	detachFace(f, detached);
	for(uint i=0; i<e.size(); i++) {
		if(getEdgeFaces(e[i]).empty()) detachEdgeRecursive(e[i], detached);
	}
}

void DelaunayMesh::detachBodyRecursive(const uint b, CellSet &detached) {
	const Buffer<uint> f = getBodyFaces(b);
	detachBody(b, detached);
	for(uint i=0; i<f.size(); i++) {
		if(getFaceBodies(f[i]).empty()) detachFaceRecursive(f[i], detached);
	}
}

void DelaunayMesh::detachQuadRecursive(const uint q, CellSet &detached) {
	const Buffer<uint> b = getQuadBodies(q);
	detachQuad(q, detached);
	for(uint i=0; i<b.size(); i++) {
		if(getBodyQuads(b[i]).empty()) detachBodyRecursive(b[i], detached);
	}
}

void DelaunayMesh::removeDetached(CellSet &detached) {
	if(detached.ns > 0) {
		const uint s = m_nsize - detached.ns; // size arfer remove
		for(uint i=0; i<detached.ns; ) { // re-organize detached cells
			const uint j = detached.n[i]; // id of detached cell
			if(j >= s) { detached.n[i] = detached.n[j-s]; detached.n[j-s] = j; } // conditionally swap i:th and (j-s)th = j
			if(j <= i + s) ++i; // go forward
		}
		while(detached.ns > 0) removeNode(detached.n[--detached.ns]);
	}
	if(detached.es > 0) {
		const uint s = m_esize - detached.es; // size arfer remove
		for(uint i=0; i<detached.es; ) { // re-organize detached cells
			const uint j = detached.e[i]; // id of detached cell
			if(j >= s) { detached.e[i] = detached.e[j-s]; detached.e[j-s] = j; } // conditionally swap i:th and (j-s)th = j
			if(j <= i + s) ++i; // go forward
		}
		while(detached.es > 0) removeEdge(detached.e[--detached.es]);
	}
	if(detached.fs > 0) {
		const uint s = m_fsize - detached.fs; // size arfer remove
		for(uint i=0; i<detached.fs; ) { // re-organize detached cells
			const uint j = detached.f[i]; // id of detached cell
			if(j >= s) { detached.f[i] = detached.f[j-s]; detached.f[j-s] = j; } // conditionally swap i:th and (j-s)th = j
			if(j <= i + s) ++i; // go forward
		}
		while(detached.fs > 0) removeFace(detached.f[--detached.fs]);
	}
	if(detached.bs > 0) {
		const uint s = m_bsize - detached.bs; // size arfer remove
		for(uint i=0; i<detached.bs; ) { // re-organize detached cells
			const uint j = detached.b[i]; // id of detached cell
			if(j >= s) { detached.b[i] = detached.b[j-s]; detached.b[j-s] = j; } // conditionally swap i:th and (j-s)th = j
			if(j <= i + s) ++i; // go forward
		}
		while(detached.bs > 0) removeBody(detached.b[--detached.bs]);
	}
	if(detached.qs > 0) {
		const uint s = m_qsize - detached.qs; // size arfer remove
		for(uint i=0; i<detached.qs; ) { // re-organize detached cells
			const uint j = detached.q[i]; // id of detached cell
			if(j >= s) { detached.q[i] = detached.q[j-s]; detached.q[j-s] = j; } // conditionally swap i:th and (j-s)th = j
			if(j <= i + s) ++i; // go forward
		}
		while(detached.qs > 0) removeQuad(detached.q[--detached.qs]);
	}
}

bool DelaunayMesh::includesCell(const Vector4 &p, uint &curr) const {
	if(m_nsize == 0) return false;
	if(m_esize == 0) { // 0-dimensional mesh
		if(curr >= m_nsize) curr = 0;
		return (getNodePosition(curr) - p).lensq() < ORTOGONALSQ;
	}
	uint i;
	if(m_fsize == 0) { // 1-dimensional mesh
		if(curr >= m_esize) curr = 0;
		uint prev = NONE;
		while(true) {
			const Buffer<uint> &n = getEdgeNodes(curr);
			for(i=0; i<n.size(); i++) {
				if(n[i] == prev) continue;
				const Vector4 v = p - getNodePosition(n[i]);
				if(v.lensq() < ORTOGONALSQ || v.dot(getEdgeDirection(curr, n[i])) > 0.0) continue;

				const Buffer<uint> &e = getNodeEdges(n[i]);
				if(e.size() < 2) {
					curr = n[i];
					return false;
				}
				if(e[0] == curr) curr = e[1];
				else curr = e[0];
				prev = n[i];
				break;
			}
			if(i == n.size()) return true;
		}
	}
	if(m_bsize == 0) { // 2-dimensional mesh
		if(curr >= m_fsize) curr = 0;
		uint prev = NONE;
		while(true) {
			const Buffer<uint> &e = getFaceEdges(curr);
			for(i=0; i<e.size(); i++) {
				if(e[i] == prev) continue;
				const Vector4 v = getEdgeDeviation(e[i], p);
				if(v.lensq() < ORTOGONALSQ || v.dot(getFaceDirection(curr, e[i])) > 0.0) continue;

				const Buffer<uint> &f = getEdgeFaces(e[i]);
				if(f.size() < 2) {
					curr = e[i];
					return false;
				}
				if(f[0] == curr) curr = f[1];
				else curr = f[0];
				prev = e[i];
				break;
			}
			if(i == e.size()) return true;
		}
	}
	if(m_qsize == 0) { // 3-dimensional mesh
		if(curr >= m_bsize) curr = 0;
		uint prev = NONE;
		while(true) {
			const Buffer<uint> &f = getBodyFaces(curr);
			for(i=0; i<f.size(); i++) {
				if(f[i] == prev) continue;
				const Vector4 v = getFaceDeviation(f[i], p);
				if(v.lensq() < ORTOGONALSQ || v.dot(getBodyDirection(curr, f[i])) > 0.0) continue;

				const Buffer<uint> &b = getFaceBodies(f[i]);
				if(b.size() < 2) {
					curr = f[i];
					return false;
				}
				if(b[0] == curr) curr = b[1];
				else curr = b[0];
				prev = f[i];
				break;
			}
			if(i == f.size()) return true;
		}
	}
	// 4-dimensional mesh
	if(curr >= m_qsize) curr = 0;
	uint prev = NONE;
	while(true) {
		const Buffer<uint> &b = getQuadBodies(curr);
		for(i=0; i<b.size(); i++) {
			if(b[i] == prev) continue;
			const Vector4 v = getBodyDeviation(b[i], p);
			if(v.lensq() < ORTOGONALSQ || v.dot(getQuadDirection(curr, b[i])) > 0.0) continue;

			const Buffer<uint> &q = getBodyQuads(b[i]);
			if(q.size() < 2) {
				curr = b[i];
				return false;
			}
			if(q[0] == curr) curr = q[1];
			else curr = q[0];
			prev = b[i];
			break;
		}
		if(i == b.size()) return true;
	}
}

//Finds the index of the element containing p in 2D or 3D mesh. Added by Jonni Lohi, University of Jyväskylä, 2023.
bool DelaunayMesh::findElement(const Vector4& p, uint& elementIndex) const {
	uint prev = elementIndex;
	if (m_bsize == 0) { //assume 2D mesh
		if (!includesCell(p, elementIndex)) {
			uint node = searchNode(p, getFaceNodes(prev)[0]);
			const Buffer<uint> f = getNodeFaces(node);
			for (uint i = 0; i < f.size(); i++) {
				uint face = f[i];
				if (includesCell(p, face)) {
					elementIndex = face;
					return true;
				}
			}
			elementIndex = prev;
			return false;
		}
		return true;
	}
	else { //assume 3D mesh
		if (!includesCell(p, elementIndex)) {
			uint node = searchNode(p, getBodyNodes(prev)[0]);
			const Buffer<uint> b = getNodeBodies(node);
			for (uint i = 0; i < b.size(); i++) {
				uint body = b[i];
				if (includesCell(p, body)) {
					elementIndex = body;
					return true;
				}
			}
			elementIndex = prev;
			return false;
		}
		return true;
	}
}

void DelaunayMesh::increaseDimension() {
	uint i, j;
	const uint nsize = m_nsize - 1;
	const uint esize = m_esize;
	const uint fsize = m_fsize;
	const uint bsize = m_bsize;

	for(i=0; i<nsize; i++) {
		addEdge(i, nsize);
	}
	for(i=0; i<esize; i++) {
		const Buffer<uint> &en = getEdgeNodes(i);
		Buffer<uint> fe(en.size() + 1);
		for(j=0; j<en.size(); j++) fe[j] = esize + en[j];
		fe[j] = i;
		addFace(fe);
	}
	for(i=0; i<fsize; i++) {
		const Buffer<uint> &fe = getFaceEdges(i);
		Buffer<uint> bf(fe.size() + 1);
		for(j=0; j<fe.size(); j++) bf[j] = fsize + fe[j];
		bf[j] = i;
		addBody(bf);
	}
	for(i=0; i<bsize; i++) {
		const Buffer<uint> &bf = getBodyFaces(i);
		Buffer<uint> qb(bf.size() + 1);
		for(j=0; j<bf.size(); j++) qb[j] = bsize + bf[j];
		qb[j] = i;
		addQuad(qb);
	}
}

uint DelaunayMesh::mergeFace(const uint f, CellSet &detached) {
	uint i, j;
	Buffer<uint> &e = m_f[f].e;
	uint es = e.size();
	if(es < 3) return f; // f is detached already
	const Buffer<uint> &fb = getFaceBodies(f);
	Vector4 fp(0,0,0,0);
	double fsq = 0.0;
	if(fb.empty()) {
		fp = getFacePosition(f);
		fsq = getRadiusSq(fp, getFaceNodes(f));
	}
	for(i=0; i<es; i++) {
		const Buffer<uint> &ef = getEdgeFaces(e[i]);
		if(ef.size() != 2) continue; // no unique face to merge
		const uint otherf = (ef[0] == f ? ef[1] : ef[0]);

		if(fb.empty()) {
			if(!isInsideSphere(fp, fsq, getFaceNodes(otherf))) continue;
		}
		else {
			Buffer<uint> &otherfb = m_f[otherf].b;
			if(!fb.isAnagram(otherfb)) continue;
			if(fb.size() == 1 && getFaceDeviation(f, getNodeAverage(getFaceNodes(otherf))).lensq() > ORTOGONALSQ) continue;
			for(j=0; j<fb.size(); j++) m_b[fb[j]].f.eraseFirst(otherf); // disconnect otherf from bodies
			otherfb.clear(); // disconnect bodies from otherf
		}

		// merge face ->
		Buffer<uint> &otherfe = m_f[otherf].e;
		for(j=otherfe.size(); j-->0; ) {
			if(e.gatherOrUngather(otherfe[j], es)) { // insert otherfe[j] to the boundary of f
				m_e[otherfe[j]].f.replaceFirst(otherf, f);
				otherfe.erase(j);
				continue;
			}
			// remove otherfe[j] from the boundary of f
			m_e[otherfe[j]].f.eraseFirst(f);
			if(i + 1 != 0) i--;
		}
		detachFaceRecursive(otherf, detached);
		e.resize(es);
		orderFaceEdges(f);
		for(j=0; j<fb.size(); j++) orderBodyFaces(fb[j]);
	}
	return f;
}

uint DelaunayMesh::mergeBody(const uint b, CellSet &detached) {
	uint i, j;

	Buffer<uint> &f = m_b[b].f;
	uint fs = f.size();
	if(fs < 4) return b; // b is detached already
	bool merged = false;
	const Buffer<uint> &bq = getBodyQuads(b);
	Vector4 bp(0,0,0,0);
	double bsq = 0.0;
	if(bq.empty()) {
		bp = getBodyPosition(b);
		bsq = getRadiusSq(bp, getBodyNodes(b));
	}
	for(i=0; i<fs; i++) {
		const Buffer<uint> &fb = getFaceBodies(f[i]);
		if(fb.size() != 2) continue; // no unique body to merge
		const uint otherb = (fb[0] == b ? fb[1] : fb[0]);

		if(bq.empty()) {
			if(!isInsideSphere(bp, bsq, getBodyNodes(otherb))) continue;
		}
		else {
			Buffer<uint> &otherbq = m_b[otherb].q;
			if(!bq.isAnagram(otherbq)) continue;
			if(bq.size() == 1 && getBodyDeviation(b, getNodeAverage(getBodyNodes(otherb))).lensq() > ORTOGONALSQ) continue;
			for(j=0; j<bq.size(); j++) m_q[bq[j]].b.eraseFirst(otherb); // disconnect otherb from quads
			otherbq.clear(); // disconnect quads from otherb
		}

		// merge body ->
		Buffer<uint> otherbf = m_b[otherb].f;
		for(j=otherbf.size(); j-->0; ) {
			if(f.gatherOrUngather(otherbf[j], fs)) { // insert otherbf[j] to the boundary of b
				m_f[otherbf[j]].b.replaceFirst(otherb, b);
				otherbf.erase(j);
				continue;
			}
			// remove otherbf[j] from the boundary of b
			m_f[otherbf[j]].b.eraseFirst(b);
			if(i + 1 != 0) i--;
		}
		detachBodyRecursive(otherb, detached);
		f.resize(fs);
		orderBodyFaces(b);
		for(j=0; j<bq.size(); j++) orderQuadBodies(bq[j]);
		merged = true;
	}
	if(merged) {
		const Buffer<uint> ff = f;
		for(i=0; i<ff.size(); i++) mergeFace(ff[i], detached);
	}
	return b;
}

uint DelaunayMesh::mergeQuad(const uint q, CellSet &detached) {
	uint i, j;
	Buffer<uint> &b = m_q[q].b;
	uint bs = b.size();
	if(bs < 5) return q; // q is detached already
	bool merged = false;
	const Vector4 qp = getQuadPosition(q);
	const double qsq = getRadiusSq(qp, getQuadNodes(q));
	for(i=0; i<bs; i++) {
		const Buffer<uint> &bq = getBodyQuads(b[i]);
		if(bq.size() != 2) continue; // no unique quad to merge
		const uint otherq = (bq[0] == q ? bq[1] : bq[0]);

		if(!isInsideSphere(qp, qsq, getQuadNodes(otherq))) continue;

		// merge quad ->
		Buffer<uint> otherqb = m_q[otherq].b;
		for(j=otherqb.size(); j-->0; ) {
			if(b.gatherOrUngather(otherqb[j], bs)) { // insert otherqb[j] to the boundary of q
				m_b[otherqb[j]].q.replaceFirst(otherq, q);
				otherqb.erase(j);
				continue;
			}
			// remove otherqb[j] from the boundary of q
			m_b[otherqb[j]].q.eraseFirst(q);
			if(i + 1 != 0) i--;
		}
		detachQuadRecursive(otherq, detached);
		b.resize(bs);
		orderQuadBodies(q);
		merged = true;
	}
	if(merged) {
		const Buffer<uint> bb = b;
		for(i=0; i<bb.size(); i++) mergeBody(bb[i], detached);
	}
	return q;
}

Vector4 DelaunayMesh::getNodeAverage(const Buffer<uint> &n) const {
	Vector4 sum(0,0,0,0);
	for(uint i=0; i<n.size(); i++) sum += getNodePosition(n[i]);
	return sum / double(n.size());
}

double DelaunayMesh::getLengthSq(const Vector4 &p, const uint node) const {
	switch(m_dim) {
	case 1: {
		const double r = p.x - getNodePosition1(node);
		return r * getTransformed1(r, 0);
	}
	case 2: {
		const Vector2 r = p.toVector2() - getNodePosition2(node);
		return r.dot(getTransformed2(r, 0));
	}
	case 3: {
		const Vector3 r = p.toVector3() - getNodePosition3(node);
		return r.dot(getTransformed3(r, 0));
	}
	default: {
		const Vector4 r = p - getNodePosition4(node);
		return r.dot(getTransformed4(r, 0));
	}
	}
}

double DelaunayMesh::getRadiusSq(const Vector4 &p, const uint node) const {
	if(node < m_w.size()) return getLengthSq(p, node) + m_w[node];
	return getLengthSq(p, node);
}

double DelaunayMesh::getRadiusSq(const Vector4 &p, const Buffer<uint> &n) const {
	double sq = 0.0;
	for(uint i=0; i<n.size(); i++) sq += getRadiusSq(p, n[i]);
	return sq / double(n.size());
}

bool DelaunayMesh::isInsideSphere(const Vector4 &p, const double sq, const uint node) const {
	return getRadiusSq(p, node) <= RADIUSSCALE * sq;
}

bool DelaunayMesh::isInsideSphere(const Vector4 &p, const double sq, const Buffer<uint> &n) const {
	const double safeSq = RADIUSSCALE * sq;
	for(uint i=0; i<n.size(); i++) {
		if(getRadiusSq(p, n[i]) > safeSq) return false;
	}
	return true;
}

bool DelaunayMesh::isOutsideSphere(const Vector4 &p, const double sq, const Buffer<uint> &n) const {
	const double safeSq = sq / RADIUSSCALE;
	for(uint i=0; i<n.size(); i++) {
		if(getRadiusSq(p, n[i]) < safeSq) return false;
	}
	return true;
}

uint DelaunayMesh::eraseNodesInside(const Vector4 &p, const double sq, uint near) {
	while(true) {
		near = searchNode(p, near);
		if(near == NONE || !isInsideSphere(p, sq, near)) return near;

		// erase near and set neighbor to near
		const uint toremove = near;
		const Buffer<uint> &e = getNodeEdges(near);
		near = (e.empty() ? 0 : getEdgeOtherNode(e[0], near));
		eraseNode(toremove);
	}
	return near;
}


