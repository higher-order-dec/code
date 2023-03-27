#include "BuilderMesh.hpp"
#include "../Types/Matrix.hpp"
#include <iostream>

using namespace gfd;
using namespace std;

BuilderMesh::BuilderMesh(const uint dim)
: DelaunayMesh(dim)
{
}

void BuilderMesh::clear() {
	DelaunayMesh::clear();
}

void BuilderMesh::createCopy(const Mesh &mesh) {
	clear();
	combine(mesh);
}

void BuilderMesh::createGrid(const Vector4 &minp, const Vector4 &maxp, const double h) {
	clear();
	const Vector4 d = maxp - minp;
	const uint nx = uint(std::abs(d.x) / h + 0.5);
	const uint ny = uint(std::abs(d.y) / h + 0.5);
	const uint nz = uint(std::abs(d.z) / h + 0.5);
	const uint nt = uint(std::abs(d.t) / h + 0.5);
	addNode(minp);
	if(nx > 0) stretchLinear(Vector4(d.x,0,0,0), nx, 0, 0);
	if(ny > 0) stretchLinear(Vector4(0,d.y,0,0), ny, 0, 0);
	if(nz > 0) stretchLinear(Vector4(0,0,d.z,0), nz, 0, 0);
	if(nt > 0) stretchLinear(Vector4(0,0,0,d.t), nt, 0, 0);
}

void BuilderMesh::createGrid(const Vector4 &minp, const Vector4 &maxp, const Vector4 &h) {
	clear();
	const Vector4 d = maxp - minp;
	const uint nx = uint(std::abs(d.x) / h.x + 0.5);
	const uint ny = uint(std::abs(d.y) / h.y + 0.5);
	const uint nz = uint(std::abs(d.z) / h.z + 0.5);
	const uint nt = uint(std::abs(d.t) / h.t + 0.5);
	addNode(minp);
	if(nx > 0) stretchLinear(Vector4(d.x,0,0,0), nx, 0, 0);
	if(ny > 0) stretchLinear(Vector4(0,d.y,0,0), ny, 0, 0);
	if(nz > 0) stretchLinear(Vector4(0,0,d.z,0), nz, 0, 0);
	if(nt > 0) stretchLinear(Vector4(0,0,0,d.t), nt, 0, 0);
}

void BuilderMesh::createTriangleGrid(const Vector2 &minp, const Vector2 &maxp, const double h, const bool rect) {
	clear();
	Vector2 d = maxp - minp;
	const uint xsize = uint(std::abs(d.x) / h + 0.5);
	if(xsize == 0) return createGrid(Vector4(minp.x,minp.y,0,0), Vector4(minp.x,maxp.y,0,0), std::sqrt(0.75)*h);
	const uint ysize = uint(std::abs(d.y) / (std::sqrt(0.75) * h) + 0.5);
	const uint xmax = xsize + 1;
	const uint ymax = ysize + 1;

	resizeNodeBuffer((2 * xmax + 1) * ymax / 2);
	resizeEdgeBuffer((xsize + xmax) * ymax / 2 + 2 * xmax * ysize);
	resizeFaceBuffer((xsize + xmax) * ysize);

	// create nodes
	uint xi, yi, i, j;
	for(yi=0; yi<=ysize; yi++) {
		const uint pair = (yi % 2);
		const double py = minp.y + yi * d.y / double(ysize);
		if(pair == 1) {
			addNode(Vector4((rect ? minp.x : minp.x - 0.5 * d.x / double(xsize)),py,0,0));
		}
		for(xi=0; xi<=xsize; xi++) {
			addNode(Vector4(minp.x + (rect && xi == xsize ? xi : xi + 0.5 * pair) * d.x / double(xsize),py,0,0));
		}
	}

	// create x-edges
	i = 0;
	for(yi=0; yi<=ysize; yi++) {
		const uint pair = (yi % 2);
		i++;
		for(xi=1-pair; xi<=xsize; xi++) {
			addEdge(i, i-1);
			i++;
		}
	}

	// create y-edges
	i = xmax;
	for(yi=1; yi<=ysize; yi++) {
		const uint pair = (yi % 2);
		if(pair == 1) i++;
		for(xi=0; xi<=xsize; xi++) {
			// create y-edges
			addEdge(i-pair, i-xmax-1);
			addEdge(i, i-xmax-pair);
			i++;
		}
	}

	// create faces
	i = xsize;
	j = (xsize + xmax) * ymax / 2;
	for(yi=1; yi<=ysize; yi++) {
		const uint pair = (yi % 2);
		if(pair > 0) i++;
		for(xi=0; xi<=xsize; xi++) {
			Buffer<uint> e(3);
			e[0] = j;
			if(xi > 0) {
				e[1] = j-1;
				e[2] = i-pair*xmax;
				addFace(e);
				i++;
			}
			e[1] = j+1;
			e[2] = i-xmax+pair*(xmax-1);
			addFace(e);
			j += 2;
		}
	}
}

void BuilderMesh::createHexagonGrid(const Vector2 &minp, const Vector2 &maxp, const double h) {
	const Vector2 d(h, h / std::sqrt(3.0));
	const uint xsize = uint((maxp.x - minp.x) / d.x + 0.999);
	const uint ysize = uint((maxp.y - minp.y) / d.y + 0.999);

	const uint xs = (xsize < 3 ? xsize : 3);
	const uint ys = (ysize < 3 ? ysize : 3);

	const Vector2 p0 = 0.5 * (maxp + minp - Vector2(xsize * d.x, ysize * d.y));
	const Vector2 p1 = p0 + Vector2(xs * d.x, ys * d.y);

	createTriangleGrid(p0, p1, h / 3.0, true);

	// generate nodes
	uint xi, yi;
	uint node = 0;
	for(yi=0; yi<ys; yi++) {
		const double py = p0.y + d.y * yi;
		for(xi=0; xi<xs; xi++) {
			const double px = p0.x + d.x * xi;
			if(yi > 0 && xi > 0) {
				node = Mesh::findNode(Vector4(px, py,0,0), 1e-13, node);
				if(node != NONE) eraseNode(node);
			}
			node = Mesh::findNode(Vector4(px+0.5*d.x, py+0.5*d.y,0,0), 1e-13, node);
			if(node != NONE) eraseNode(node);
		}
	}
	const Vector4 pp(p0 + 0.99999*d, 0.0, 0.0);
	if(xsize > xs) repeatMiddle(pp, Vector4(d.x,0,0,0), xsize - xs);
	if(ysize > ys) repeatMiddle(pp, Vector4(0,d.y,0,0), ysize - ys);
}

void BuilderMesh::createSnubSquareGrid(const Vector2 &minp, const Vector2 &maxp, const double h) {
	const uint xsize = uint((maxp.x - minp.x) / h + 0.5);
	const uint ysize = uint((maxp.y - minp.y) / h + 0.5);

	const uint xs = (xsize < 3 ? xsize : 3);
	const uint ys = (ysize < 3 ? ysize : 3);

//	const double SQ3 = std::sqrt(3.0);
	const double l = h / 3.0; // / (SQ3 + 1.0);// / 2.7725;
	const Vector2 p0 = 0.5 * (minp + maxp - Vector2(xsize * h, ysize * h));
	const Vector2 p1 = p0 + Vector2(xs * h, ys * h);

	createGrid(Vector4(p0,0,0), Vector4(p1,0,0), h);

	// generate nodes
	uint node = 0;
	uint xi, yi;
	for(yi=0; yi<=ys; yi++) {
		const double py = p0.y + h * yi;
		for(xi=0; xi<=xs; xi++) {
			const double px = p0.x + h * xi;
			if(yi > 0 && (xi == 0 || xi == xs)) {
				node = insertNode(Vector4(px, py - l, 0,0), 0.0, node, true);
				node = insertNode(Vector4(px, py - h + l, 0,0), 0.0, node, true);
			}

			if(xi > 0) {
				if(yi > 0) {
					node = insertNode(Vector4(px - 0.5 * h, py - 0.5 * h, 0,0), 0.0, node, true);
					if(xi > 1) node = insertNode(Vector4(px - 0.5 * h-l, py - 0.5 * h, 0,0), 0.0, node, true);
					if(yi < ys) node = insertNode(Vector4(px - 0.5 * (h + l), py - 0.5 * l, 0,0), 0.0, node, true);
					if(yi > 1) node = insertNode(Vector4(px - 0.5 * (h + l), py - h + 0.5 * l, 0,0), 0.0, node, true);
					if(xi < xs) {
						node = insertNode(Vector4(px - 0.5 * l, py - 0.5 * (h - l), 0,0), 0.0, node, true);
						node = insertNode(Vector4(px - 0.5 * l, py - 0.5 * (h + l), 0,0), 0.0, node, true);
					}
				}
				node = insertNode(Vector4(px - l, py, 0,0), 0.0, node, true);
				if(yi == 0 || yi == ys) node = insertNode(Vector4(px - 0.5 * (h + l), py, 0,0), 0.0, node, true);
			}
		}
	}

	const Vector4 pp(p0 + 0.99999*Vector2(h,h), 0.0, 0.0);
	if(xsize > xs) repeatMiddle(pp, Vector4(h,0,0,0), xsize - xs);
	if(ysize > ys) repeatMiddle(pp, Vector4(0,h,0,0), ysize - ys);
}

void BuilderMesh::createTetrilleGrid(const Vector2 &minp, const Vector2 &maxp, const double h) {
	const Vector2 d(h, h * std::sqrt(3.0));
	const uint xsize = uint((maxp.x - minp.x) / d.x + 0.999);
	const uint ysize = uint((maxp.y - minp.y) / d.y + 0.999);

	const uint xs = (xsize < 3 ? xsize : 3);
	const uint ys = (ysize < 3 ? ysize : 3);

	const Vector2 p0 = 0.5 * (maxp + minp - Vector2(xsize * d.x, ysize * d.y));
	const Vector2 p1 = p0 + Vector2(xs * d.x, ys * d.y);
	createTriangleGrid(p0, p1, 0.5 * h, true);

	// generate nodes
	uint node = 0;
	uint xi, yi;
	for(yi=0; yi<ys; yi++) {
		const double py = p0.y + d.y * yi;
		for(xi=0; xi<xs; xi++) {
			const double px = p0.x + d.x * xi;
			if(xi > 0) {
				node = insertNode(Vector4(px, py + d.y / 6.0, 0,0), 0.0, node, true);
				node = insertNode(Vector4(px, py + d.y * 5.0 / 6.0, 0,0), 0.0, node, true);
			}
			node = insertNode(Vector4(px + 0.5 * d.x, py + d.y / 3.0, 0,0), 0.0, node, true);
			node = insertNode(Vector4(px + 0.5 * d.x, py + d.y * 2.0 / 3.0, 0,0), 0.0, node, true);
		}
	}

	const Vector4 pp(p0 + 0.99999*d, 0.0, 0.0);
	if(xsize > xs) repeatMiddle(pp, Vector4(d.x,0,0,0), xsize - xs);
	if(ysize > ys) repeatMiddle(pp, Vector4(0,d.y,0,0), ysize - ys);
}

void BuilderMesh::createCircleBoundaryGrid(const double r, const double h) {
	clear();

	uint isize = uint(r * PIx2 / h + 0.5);
	if(isize < 3) isize = 3;

	uint i;
	for(i=0; i<isize; i++) {
		const double ang = i * PIx2 / double(isize);
		addNode(Vector4(r * std::cos(ang), r * std::sin(ang), 0.0, 0.0));
	}
	for(i=1; i<isize; i++) addEdge(i, i-1);
	addEdge(0, isize - 1);
}

void BuilderMesh::createFccGrid(const Vector3 &minp, const Vector3 &maxp, const double h) {
	const uint xsize = uint((maxp.x - minp.x) / h - 1e-8) + 1;
	const uint ysize = uint((maxp.y - minp.y) / h - 1e-8) + 1;
	const uint zsize = uint((maxp.z - minp.z) / h - 1e-8) + 1;

	const uint xs = (xsize < 3 ? xsize : 3);
	const uint ys = (ysize < 3 ? ysize : 3);
	const uint zs = (zsize < 3 ? zsize : 3);

	const Vector3 p0 = 0.5 * (maxp + minp - Vector3(xsize * h, ysize * h, zsize * h));
	const Vector3 p1 = p0 + Vector3(xs * h, ys * h, zs * h);

	createGrid(Vector4(p0,0), Vector4(p1,0), h);

	uint node = 0;
	uint xi, yi, zi;
	for(zi=0; zi<=zs; zi++) {
		const double pz = p0.z + h * zi;
		for(yi=0; yi<=ys; yi++) {
			const double py = p0.y + h * yi;
			for(xi=0; xi<=xs; xi++) {
				const double px = p0.x + h * xi;
				if(xi < xs && yi < ys) node = insertNode(Vector4(px + 0.5 * h, py + 0.5 * h, pz, 0), 0.0, node, true);
				if(xi < xs && zi < zs) node = insertNode(Vector4(px + 0.5 * h, py, pz + 0.5 * h, 0), 0.0, node, true);
				if(yi < ys && zi < zs) node = insertNode(Vector4(px, py + 0.5 * h, pz + 0.5 * h, 0), 0.0, node, true);
			}
		}
	}

	const Vector4 pp(p0 + 0.99999 * Vector3(h, h, h), 0.0);
	if(xsize > xs) repeatMiddle(pp, Vector4(h,0,0,0), xsize - xs);
	if(ysize > ys) repeatMiddle(pp, Vector4(0,h,0,0), ysize - ys);
	if(zsize > zs) repeatMiddle(pp, Vector4(0,0,h,0), zsize - zs);
}

void BuilderMesh::createBccGrid(const Vector3 &minp, const Vector3 &maxp, const double h) {
	const uint xsize = uint((maxp.x - minp.x) / h - 1e-8) + 1;
	const uint ysize = uint((maxp.y - minp.y) / h - 1e-8) + 1;
	const uint zsize = uint((maxp.z - minp.z) / h - 1e-8) + 1;

	const uint xs = (xsize < 3 ? xsize : 3);
	const uint ys = (ysize < 3 ? ysize : 3);
	const uint zs = (zsize < 3 ? zsize : 3);

	const Vector3 p0 = 0.5 * (maxp + minp - Vector3(xsize * h, ysize * h, zsize * h));
	const Vector3 p1 = p0 + Vector3(xs * h, ys * h, zs * h);

	createGrid(Vector4(p0,0), Vector4(p1,0), h);

	uint node = 0;
	uint xi, yi, zi;
	for(zi=0; zi<=zs; zi++) {
		const double pz = p0.z + h * zi;
		for(yi=0; yi<=ys; yi++) {
			const double py = p0.y + h * yi;
			for(xi=0; xi<=xs; xi++) {
				const double px = p0.x + h * xi;
				if(xi < xs && yi < ys && zi < zs) node = insertNode(Vector4(px + 0.5 * h, py + 0.5 * h, pz + 0.5 * h, 0), 0.0, node, true);
				if(xi < xs && yi < ys && (zi == 0 || zi == zs)) node = insertNode(Vector4(px + 0.5 * h, py + 0.5 * h, pz, 0), 0.0, node, true);
				if(xi < xs && zi < zs && (yi == 0 || yi == ys)) node = insertNode(Vector4(px + 0.5 * h, py, pz + 0.5 * h, 0), 0.0, node, true);
				if(zi < zs && yi < ys && (xi == 0 || xi == xs)) node = insertNode(Vector4(px, py + 0.5 * h, pz + 0.5 * h, 0), 0.0, node, true);
				if(xi < xs && (yi == 0 || yi == ys || zi == 0 || zi == zs)) node = insertNode(Vector4(px + 0.5 * h, py, pz, 0), 0.0, node, true);
				if(yi < ys && (xi == 0 || xi == xs || zi == 0 || zi == zs)) node = insertNode(Vector4(px, py + 0.5 * h, pz, 0), 0.0, node, true);
				if(zi < zs && (yi == 0 || yi == ys || xi == 0 || xi == xs)) node = insertNode(Vector4(px, py, pz + 0.5 * h, 0), 0.0, node, true);
			}
		}
	}
	const Vector4 pp(p0 + 0.99999 * Vector3(h, h, h), 0.0);
	if(xsize > xs) repeatMiddle(pp, Vector4(h,0,0,0), xsize - xs);
	if(ysize > ys) repeatMiddle(pp, Vector4(0,h,0,0), ysize - ys);
	if(zsize > zs) repeatMiddle(pp, Vector4(0,0,h,0), zsize - zs);
}

void BuilderMesh::createTruncatedOctahedraGrid(const Vector3 &minp, const Vector3 &maxp, const double h) {
	const uint xsize = uint((maxp.x - minp.x) / h - 1e-8) + 1;
	const uint ysize = uint((maxp.y - minp.y) / h - 1e-8) + 1;
	const uint zsize = uint((maxp.z - minp.z) / h - 1e-8) + 1;

	const uint xs = (xsize < 3 ? xsize : 3);
	const uint ys = (ysize < 3 ? ysize : 3);
	const uint zs = (zsize < 3 ? zsize : 3);

	const Vector3 p0 = 0.5 * (maxp + minp - Vector3(xsize * h, ysize * h, zsize * h));
	const Vector3 p1 = p0 + Vector3(xs * h, ys * h, zs * h);

	createGrid(Vector4(p0,0), Vector4(p1,0), h);

	uint node = 0;
	uint xi, yi, zi;
	for(zi=0; zi<=zs; zi++) {
		const double pz = p0.z + h * zi;
		for(yi=0; yi<=ys; yi++) {
			const double py = p0.y + h * yi;
			for(xi=0; xi<=xs; xi++) {
				const double px = p0.x + h * xi;
				if(xi > 0 && yi > 0 && zi > 0) {
					node = insertNode(Vector4(px - 0.5 * h, py - 0.5 * h, pz - 0.5 * h, 0), 0.0, node, true);
					node = insertNode(Vector4(px - 0.25 * h, py - 0.5 * h, pz - 0.75 * h, 0), 0.0, node, true);
					node = insertNode(Vector4(px - 0.75 * h, py - 0.75 * h, pz - 0.5 * h, 0), 0.0, node, true);
					node = insertNode(Vector4(px - 0.75 * h, py - 0.25 * h, pz - 0.5 * h, 0), 0.0, node, true);
					node = insertNode(Vector4(px - 0.25 * h, py - 0.5 * h, pz - 0.25 * h, 0), 0.0, node, true);
				}
				if(yi > 0 && zi > 0) node = insertNode(Vector4(px, py - 0.5 * h, pz - 0.5 * h, 0), 0.0, node, true);
				if(xi > 0 && yi > 0) {
					node = insertNode(Vector4(px - 0.25 * h, py - 0.25 * h, pz, 0), 0.0, node, true);
					node = insertNode(Vector4(px - 0.25 * h, py - 0.75 * h, pz, 0), 0.0, node, true);
				}
				if(xi > 0 && zi > 0) {
					node = insertNode(Vector4(px - 0.75 * h, py, pz - 0.75 * h, 0), 0.0, node, true);
					node = insertNode(Vector4(px - 0.75 * h, py, pz - 0.25 * h, 0), 0.0, node, true);
				}
				if(xi > 0) node = insertNode(Vector4(px - 0.5 * h, py, pz, 0), 0.0, node, true);
			}
		}
	}

	const Vector4 pp(p0 + 0.99999 * Vector3(h, h, h), 0.0);
	if(xsize > xs) repeatMiddle(pp, Vector4(h,0,0,0), xsize - xs);
	if(ysize > ys) repeatMiddle(pp, Vector4(0,h,0,0), ysize - ys);
	if(zsize > zs) repeatMiddle(pp, Vector4(0,0,h,0), zsize - zs);
}

void BuilderMesh::createA15Grid(const Vector3 &minp, const Vector3 &maxp, const double h) {
	const uint xsize = uint((maxp.x - minp.x) / h - 1e-8) + 1;
	const uint ysize = uint((maxp.y - minp.y) / h - 1e-8) + 1;
	const uint zsize = uint((maxp.z - minp.z) / h - 1e-8) + 1;

	const uint xs = (xsize < 3 ? xsize : 3);
	const uint ys = (ysize < 3 ? ysize : 3);
	const uint zs = (zsize < 3 ? zsize : 3);

	const Vector3 p0 = 0.5 * (maxp + minp - Vector3(xsize * h, ysize * h, zsize * h));
	const Vector3 p1 = p0 + Vector3(xs * h, ys * h, zs * h);

	createGrid(Vector4(p0,0), Vector4(p1,0), h);

	const double bias = 0.25;
	const double a = 0.5 - bias;
	const double b = 0.5 + bias;

	const double weight = -0.020292 * h * h;

	uint node = 0;
	uint xi, yi, zi;
	for(zi=0; zi<=zs; zi++) {
		const double pz = p0.z + h * zi;
		for(yi=0; yi<=ys; yi++) {
			const double py = p0.y + h * yi;
			for(xi=0; xi<=xs; xi++) {
				const double px = p0.x + h * xi;
				if(xi < xs && yi < ys && zi < zs) node = insertNode(Vector4(px + 0.5 * h, py + 0.5 * h, pz + 0.5 * h, 0), 0.0, node, true);

				if(xi < xs && yi < ys) {
					if(zi == 0 || zi == zs) {
						node = insertNode(Vector4(px + h / 3.0, py + h / 3.0, pz, 0), 0.0, node, true);
						node = insertNode(Vector4(px + 2*h / 3.0, py + h / 3.0, pz, 0), 0.0, node, true);
						node = insertNode(Vector4(px + h / 3.0, py + 2*h / 3.0, pz, 0), 0.0, node, true);
						node = insertNode(Vector4(px + 2*h / 3.0, py + 2*h / 3.0, pz, 0), 0.0, node, true);
					}
					else {
						node = insertNode(Vector4(px + 0.5 * h, py + a * h, pz, 0), weight, node, true);
						node = insertNode(Vector4(px + 0.5 * h, py + b * h, pz, 0), weight, node, true);
					}
				}
				if(xi < xs && zi < zs) {
					if(yi == 0 || yi == ys) {
						node = insertNode(Vector4(px + h / 3.0, py, pz + h / 3.0, 0), 0.0, node, true);
						node = insertNode(Vector4(px + 2*h / 3.0, py, pz + h / 3.0, 0), 0.0, node, true);
						node = insertNode(Vector4(px + h / 3.0, py, pz + 2*h / 3.0, 0), 0.0, node, true);
						node = insertNode(Vector4(px + 2*h / 3.0, py, pz + 2*h / 3.0, 0), 0.0, node, true);
					}
					else {
						node = insertNode(Vector4(px + a * h, py, pz + 0.5 * h, 0), weight, node, true);
						node = insertNode(Vector4(px + b * h, py, pz + 0.5 * h, 0), weight, node, true);
					}
				}
				if(zi < zs && yi < ys) {
					if(xi == 0 || xi == xs) {
						node = insertNode(Vector4(px, py + h / 3.0, pz + h / 3.0, 0), 0.0, node, true);
						node = insertNode(Vector4(px, py + 2*h / 3.0, pz + h / 3.0, 0), 0.0, node, true);
						node = insertNode(Vector4(px, py + h / 3.0, pz + 2*h / 3.0, 0), 0.0, node, true);
						node = insertNode(Vector4(px, py + 2*h / 3.0, pz + 2*h / 3.0, 0), 0.0, node, true);
					}
					else {
						node = insertNode(Vector4(px, py + 0.5 * h, pz + a * h, 0), weight, node, true);
						node = insertNode(Vector4(px, py + 0.5 * h, pz + b * h, 0), weight, node, true);
					}

				}
				if(xi < xs && (yi == 0 || yi == ys || zi == 0 || zi == zs)) {
					node = insertNode(Vector4(px + h / 3.0, py, pz, 0), 0.0, node, true);
					node = insertNode(Vector4(px + 2*h / 3.0, py, pz, 0), 0.0, node, true);
				}
				if(yi < ys && (xi == 0 || xi == xs || zi == 0 || zi == zs)) {
					node = insertNode(Vector4(px, py + h / 3.0, pz, 0), 0.0, node, true);
					node = insertNode(Vector4(px, py + 2*h / 3.0, pz, 0), 0.0, node, true);
				}
				if(zi < zs && (yi == 0 || yi == ys || xi == 0 || xi == xs)) {
					node = insertNode(Vector4(px, py, pz + h / 3.0, 0), 0.0, node, true);
					node = insertNode(Vector4(px, py, pz + 2*h / 3.0, 0), 0.0, node, true);
				}
			}
		}
	}

	const Vector4 pp(p0 + 0.99999 * Vector3(h, h, h), 0.0);
	if(xsize > xs) repeatMiddle(pp, Vector4(h,0,0,0), xsize - xs);
	if(ysize > ys) repeatMiddle(pp, Vector4(0,h,0,0), ysize - ys);
	if(zsize > zs) repeatMiddle(pp, Vector4(0,0,h,0), zsize - zs);
}

void BuilderMesh::createC15Grid(const Vector3 &minp, const Vector3 &maxp, const double h) {
	const uint xsize = uint((maxp.x - minp.x) / h - 1e-8) + 1;
	const uint ysize = uint((maxp.y - minp.y) / h - 1e-8) + 1;
	const uint zsize = uint((maxp.z - minp.z) / h - 1e-8) + 1;

	const uint xs = (xsize < 3 ? xsize : 3);
	const uint ys = (ysize < 3 ? ysize : 3);
	const uint zs = (zsize < 3 ? zsize : 3);

	const Vector3 p0 = 0.5 * (maxp + minp - Vector3(xsize * h, ysize * h, zsize * h));
	const Vector3 p1 = p0 + Vector3(xs * h, ys * h, zs * h);

	//clear();
	createGrid(Vector4(p0,0), Vector4(p1,0), h);

	const double weight = 0.017137 * h * h;
	uint node = 0;
	uint xi, yi, zi;
	for(zi=0; zi<=zs; zi++) {
		const double pz = p0.z + zi * h;
		for(yi=0; yi<=ys; yi++) {
			const double py = p0.y + yi * h;
			for(xi=0; xi<=xs; xi++) {
				const double px = p0.x + xi * h;
				if(yi < ys && zi < zs) {
					node = insertNode(Vector4(px, py + 0.5 * h, pz + 0.5 * h, 0), 0.0, node, true);
					if(xi == 0) {
						node = insertNode(Vector4(px, py + 0.25 * h, pz + 0.25 * h, 0), 0.0, node, true);
						node = insertNode(Vector4(px, py + 0.75 * h, pz + 0.75 * h, 0), 0.0, node, true);

						if(yi > 0) node = insertNode(Vector4(px, py + 0.125 * h, pz + 0.625 * h, 0), 0.0, node, true);
						if(zi > 0) node = insertNode(Vector4(px, py + 0.625 * h, pz + 0.125 * h, 0), 0.0, node, true);
						if(yi + 1 < ys) node = insertNode(Vector4(px, py + 0.875 * h, pz + 0.375 * h, 0), 0.0, node, true);
						if(zi + 1 < zs) node = insertNode(Vector4(px, py + 0.375 * h, pz + 0.875 * h, 0), 0.0, node, true);
					}
					else if(xi == xs) {
						node = insertNode(Vector4(px, py + 0.25 * h, pz + 0.75 * h, 0), 0.0, node, true);
						node = insertNode(Vector4(px, py + 0.75 * h, pz + 0.25 * h, 0), 0.0, node, true);

						if(yi > 0) node = insertNode(Vector4(px, py + 0.125 * h, pz + 0.375 * h, 0), 0.0, node, true);
						if(zi > 0) node = insertNode(Vector4(px, py + 0.375 * h, pz + 0.125 * h, 0), 0.0, node, true);
						if(yi + 1 < ys) node = insertNode(Vector4(px, py + 0.875 * h, pz + 0.625 * h, 0), 0.0, node, true);
						if(zi + 1 < zs) node = insertNode(Vector4(px, py + 0.625 * h, pz + 0.875 * h, 0), 0.0, node, true);
					}
				}
				if(xi < xs && zi < zs) {
					node = insertNode(Vector4(px + 0.5 * h, py, pz + 0.5 * h, 0), 0.0, node, true);
					if(yi == 0) {
						node = insertNode(Vector4(px + 0.25 * h, py, pz + 0.25 * h, 0), 0.0, node, true);
						node = insertNode(Vector4(px + 0.75 * h, py, pz + 0.75 * h, 0), 0.0, node, true);

						if(xi > 0) node = insertNode(Vector4(px + 0.125 * h, py, pz + 0.625 * h, 0), 0.0, node, true);
						if(zi > 0) node = insertNode(Vector4(px + 0.625 * h, py, pz + 0.125 * h, 0), 0.0, node, true);
						if(xi + 1 < xs) node = insertNode(Vector4(px + 0.875 * h, py, pz + 0.375 * h, 0), 0.0, node, true);
						if(zi + 1 < zs) node = insertNode(Vector4(px + 0.375 * h, py, pz + 0.875 * h, 0), 0.0, node, true);
					}
					else if(yi == ys) {
						node = insertNode(Vector4(px + 0.25 * h, py, pz + 0.75 * h, 0), 0.0, node, true);
						node = insertNode(Vector4(px + 0.75 * h, py, pz + 0.25 * h, 0), 0.0, node, true);

						if(xi > 0) node = insertNode(Vector4(px + 0.125 * h, py, pz + 0.375 * h, 0), 0.0, node, true);
						if(zi > 0) node = insertNode(Vector4(px + 0.375 * h, py, pz + 0.125 * h, 0), 0.0, node, true);
						if(xi + 1 < xs) node = insertNode(Vector4(px + 0.875 * h, py, pz + 0.625 * h, 0), 0.0, node, true);
						if(zi + 1 < zs) node = insertNode(Vector4(px + 0.625 * h, py, pz + 0.875 * h, 0), 0.0, node, true);
					}
				}
				if(xi < xs && yi < ys) {
					node = insertNode(Vector4(px + 0.5 * h, py + 0.5 * h, pz, 0), 0.0, node, true);
					if(zi == 0) {
						node = insertNode(Vector4(px + 0.25 * h, py + 0.25 * h, pz, 0), 0.0, node, true);
						node = insertNode(Vector4(px + 0.75 * h, py + 0.75 * h, pz, 0), 0.0, node, true);

						if(xi > 0) node = insertNode(Vector4(px + 0.125 * h, py + 0.625 * h, pz, 0), 0.0, node, true);
						if(yi > 0) node = insertNode(Vector4(px + 0.625 * h, py + 0.125 * h, pz, 0), 0.0, node, true);
						if(xi + 1 < xs) node = insertNode(Vector4(px + 0.875 * h, py + 0.375 * h, pz, 0), 0.0, node, true);
						if(yi + 1 < ys) node = insertNode(Vector4(px + 0.375 * h, py + 0.875 * h, pz, 0), 0.0, node, true);
					}
					else if(zi == zs) {
						node = insertNode(Vector4(px + 0.25 * h, py + 0.75 * h, pz, 0), 0.0, node, true);
						node = insertNode(Vector4(px + 0.75 * h, py + 0.25 * h, pz, 0), 0.0, node, true);

						if(xi > 0) node = insertNode(Vector4(px + 0.125 * h, py + 0.375 * h, pz, 0), 0.0, node, true);
						if(yi > 0) node = insertNode(Vector4(px + 0.375 * h, py + 0.125 * h, pz, 0), 0.0, node, true);
						if(xi + 1 < xs) node = insertNode(Vector4(px + 0.875 * h, py + 0.625 * h, pz, 0), 0.0, node, true);
						if(yi + 1 < ys) node = insertNode(Vector4(px + 0.625 * h, py + 0.875 * h, pz, 0), 0.0, node, true);
					}
				}
				if(xi < xs && (yi == 0 || yi == ys) && (zi == 0 || zi == zs)) {
					node = insertNode(Vector4(px + 0.333333 * h, py, pz, 0), 0.0, node, true);
					node = insertNode(Vector4(px + 0.666667 * h, py, pz, 0), 0.0, node, true);
				}
				if(yi < ys && (xi == 0 || xi == xs) && (zi == 0 || zi == zs)) {
					node = insertNode(Vector4(px, py + 0.333333 * h, pz, 0), 0.0, node, true);
					node = insertNode(Vector4(px, py + 0.666667 * h, pz, 0), 0.0, node, true);
				}
				if(zi < zs && (yi == 0 || yi == ys) && (xi == 0 || xi == xs)) {
					node = insertNode(Vector4(px, py, pz + 0.333333 * h, 0), 0.0, node, true);
					node = insertNode(Vector4(px, py, pz + 0.666667 * h, 0), 0.0, node, true);
				}
				if(xi < xs && yi < ys && zi < zs) {
					const double w = 0.0;
					if(xi > 0 && yi + 1 < ys && zi + 1 < zs) node = insertNode(Vector4(px + 0.25 * h, py + 0.75 * h, pz + 0.75 * h, 0), 0.0, node, true);
					else node = insertNode(Vector4(px + 0.25 * h, py + 0.75 * h, pz + 0.75 * h, 0), w, node, true);
					if(yi > 0 && xi + 1 < xs && zi + 1 < zs) node = insertNode(Vector4(px + 0.75 * h, py + 0.25 * h, pz + 0.75 * h, 0), 0.0, node, true);
					else node = insertNode(Vector4(px + 0.75 * h, py + 0.25 * h, pz + 0.75 * h, 0), w, node, true);
					if(zi > 0 && xi + 1 < xs && yi + 1 < ys) node = insertNode(Vector4(px + 0.75 * h, py + 0.75 * h, pz + 0.25 * h, 0), 0.0, node, true);
					else node = insertNode(Vector4(px + 0.75 * h, py + 0.75 * h, pz + 0.25 * h, 0), w, node, true);
					if(xi > 0 && yi > 0 && zi > 0) node = insertNode(Vector4(px + 0.25 * h, py + 0.25 * h, pz + 0.25 * h, 0), 0.0, node, true);
					else node = insertNode(Vector4(px + 0.25 * h, py + 0.25 * h, pz + 0.25 * h, 0), w, node, true);

					if(xi > 0 && zi < zs-1) node = insertNode(Vector4(px + 0.125 * h, py + 0.375 * h, pz + 0.875 * h, 0), weight, node, true);
					if(xi > 0 && yi > 0) node = insertNode(Vector4(px + 0.125 * h, py + 0.125 * h, pz + 0.625 * h, 0), weight, node, true);
					if(yi > 0 && zi < zs-1) node = insertNode(Vector4(px + 0.375 * h, py + 0.125 * h, pz + 0.875 * h, 0), weight, node, true);
					node = insertNode(Vector4(px + 0.375 * h, py + 0.375 * h, pz + 0.625 * h, 0), weight, node, true);

					if(xi > 0 && yi < ys-1) node = insertNode(Vector4(px + 0.125 * h, py + 0.875 * h, pz + 0.375 * h, 0), weight, node, true);
					if(zi > 0 && yi < ys-1) node = insertNode(Vector4(px + 0.375 * h, py + 0.875 * h, pz + 0.125 * h, 0), weight, node, true);
					node = insertNode(Vector4(px + 0.375 * h, py + 0.625 * h, pz + 0.375 * h, 0), weight, node, true);
					if(xi > 0 && zi > 0) node = insertNode(Vector4(px + 0.125 * h, py + 0.625 * h, pz + 0.125 * h, 0), weight, node, true);

					if(yi > 0 && zi > 0) node = insertNode(Vector4(px + 0.625 * h, py + 0.125 * h, pz + 0.125 * h, 0), weight, node, true);
					if(xi < xs-1 && zi > 0) node = insertNode(Vector4(px + 0.875 * h, py + 0.375 * h, pz + 0.125 * h, 0), weight, node, true);
					node = insertNode(Vector4(px + 0.625 * h, py + 0.375 * h, pz + 0.375 * h, 0), weight, node, true);
					if(yi > 0 && xi < xs-1) node = insertNode(Vector4(px + 0.875 * h, py + 0.125 * h, pz + 0.375 * h, 0), weight, node, true);

					node = insertNode(Vector4(px + 0.625 * h, py + 0.625 * h, pz + 0.625 * h, 0), weight, node, true);
					if(xi < xs-1 && zi < zs-1) node = insertNode(Vector4(px + 0.875 * h, py + 0.625 * h, pz + 0.875 * h, 0), weight, node, true);
					if(xi < xs-1 && yi < ys-1) node = insertNode(Vector4(px + 0.875 * h, py + 0.875 * h, pz + 0.625 * h, 0), weight, node, true);
					if(yi < ys-1 && zi < zs-1) node = insertNode(Vector4(px + 0.625 * h, py + 0.875 * h, pz + 0.875 * h, 0), weight, node, true);
				}
			}
		}
	}

	const Vector4 pp(p0 + 0.99999 * Vector3(h, h, h), 0.0);
	if(xsize > xs) repeatMiddle(pp, Vector4(h,0,0,0), xsize - xs);
	if(ysize > ys) repeatMiddle(pp, Vector4(0,h,0,0), ysize - ys);
	if(zsize > zs) repeatMiddle(pp, Vector4(0,0,h,0), zsize - zs);
}

void BuilderMesh::createZGrid(const Vector3 &minp, const Vector3 &maxp, const double h) {
	const Vector3 d(h, h * std::sqrt(3.0), h);
	const uint xsize = uint((maxp.x - minp.x) / d.x + 0.99999);
	const uint ysize = uint((maxp.y - minp.y) / d.y + 0.99999);
	const uint zsize = uint((maxp.z - minp.z) / d.z + 0.99999);

	const uint xs = (xsize < 3 ? xsize : 3);
	const uint ys = (ysize < 3 ? ysize : 3);
	const uint zs = (zsize < 3 ? zsize : 3);

	const Vector3 p0 = 0.5 * (maxp + minp - Vector3(xsize * d.x, ysize * d.y, zsize * d.z));
	const Vector3 p1 = p0 + Vector3(xs * d.x, ys * d.y, zs * d.z);

	createTriangleGrid(Vector2(p0.x, p0.y), Vector2(p1.x, p1.y), h, true);
	move(Vector4(0,0,p0.z,0));
	stretchLinear(Vector4(0,0,p1.z-p0.z,0), 2*zs);

	uint node = 0;
	uint xi, yi, zi;
	for(zi=0; zi<zs; zi++) {
		const double pz = p0.z + d.z * zi;
		for(yi=0; yi<ys; yi++) {
			const double py = p0.y + d.y * yi;
			for(xi=0; xi<xs; xi++) {
				const double px = p0.x + d.x * xi;
				node = insertNode(Vector4(px + 0.25 * d.x, py + 0.25 * d.y, pz + 0.75 * d.z, 0), 0.0, node, true);
				node = insertNode(Vector4(px + 0.25 * d.x, py + 0.75 * d.y, pz + 0.75 * d.z, 0), 0.0, node, true);
				node = insertNode(Vector4(px + 0.75 * d.x, py + 0.25 * d.y, pz + 0.75 * d.z, 0), 0.0, node, true);
				node = insertNode(Vector4(px + 0.75 * d.x, py + 0.75 * d.y, pz + 0.75 * d.z, 0), 0.0, node, true);
				node = insertNode(Vector4(px + 0.5 * d.x, py + 1.0 / 6.0 * d.y, pz + 0.25 * d.z, 0), 0.0, node, true);
				node = insertNode(Vector4(px + 0.5 * d.x, py + 5.0 / 6.0 * d.y, pz + 0.25 * d.z, 0), 0.0, node, true);
				if(xi > 0) {
                    node = insertNode(Vector4(px, py + 0.5 * d.y, pz + 0.75 * d.z, 0), 0.0, node, true);
                    node = insertNode(Vector4(px, py + 1.0 / 3.0 * d.y, pz + 0.25 * d.z, 0), 0.0, node, true);
					node = insertNode(Vector4(px, py + 2.0 / 3.0 * d.y, pz + 0.25 * d.z, 0), 0.0, node, true);
				}
				if(yi > 0) {
                    node = insertNode(Vector4(px + 0.5 * d.x, py, pz + 0.75 * d.z, 0), 0.0, node, true);
				}
			}
		}
	}

	const Vector4 pp(p0 + 0.99999 * d, 0.0);
	if(xsize > xs) repeatMiddle(pp, Vector4(d.x,0,0,0), xsize - xs);
	if(ysize > ys) repeatMiddle(pp, Vector4(0,d.y,0,0), ysize - ys);
	if(zsize > zs) repeatMiddle(pp, Vector4(0,0,d.z,0), zsize - zs);
}

void BuilderMesh::createSphereBoundaryGrid(const double r, const double h, const uint optimizationSteps, const bool integerDivision) {
	clear();

	uint i, j, k;
	double jang = 0.0;
	uint jmax = 1;

	double jang0 = 0.0;
	uint jmax0 = 0;
	uint nodei0 = 0;
	uint edgei0 = 0;

	Buffer<uint> e(3);
	addNode(Vector4(r, 0.0, 0.0, 0.0));
	addNode(Vector4(-r, 0.0, 0.0, 0.0));

	const double itmax = 0.5 * r * PI / h;
	uint imax = uint(itmax + 0.5);
	if(imax < 1) imax = 1;
	for(i=1; i<=imax; i++) {
		double ang = 0.5 * i * PI / double(imax);
		const double icos = r * std::cos(ang);
		const double isin = r * std::sin(ang);

		const uint nodei1 = getNodeSize();
		const uint edgei1 = getEdgeSize();
		const double jtmax = isin * PIx2 / h;
		if(integerDivision) jmax *= uint(jtmax / double(jmax) + 0.5);
		else jmax = uint(jtmax + 0.5);
		if(i < imax) {
			for(j=0; j<jmax; j++) {
				ang = jang + 2.0 * j * PI / double(jmax);
				const double jcos = isin * std::cos(ang);
				const double jsin = isin * std::sin(ang);
				addNode(Vector4( icos, jcos, jsin, 0.0));
				addNode(Vector4(-icos, jcos, jsin, 0.0));
				if(j > 0) {
					addEdge(getNodeSize() - 2, getNodeSize() - 4);
					addEdge(getNodeSize() - 1, getNodeSize() - 3);
				}
			}
			addEdge(getNodeSize() - 2 * jmax, getNodeSize() - 2);
			addEdge(getNodeSize() - 2 * jmax + 1, getNodeSize() - 1);

			uint edge1 = NONE;
			uint edge2 = NONE;
			for(j=0, k=0; k<jmax; ) {
				addEdge(nodei0 + 2 * (j<jmax0?j:j-jmax0), nodei1 + 2 * k);
				addEdge(nodei0 + 2 * (j<jmax0?j:j-jmax0) + 1, nodei1 + 2 * k + 1);

				if(edge1 != NONE) {
					e[0] = getEdgeSize() - 2;
					e[1] = getEdgeSize() - 4;
					e[2] = edge1;
					addFace(e);
				}
				if(edge2 != NONE) {
					e[0] = getEdgeSize() - 1;
					e[1] = getEdgeSize() - 3;
					e[2] = edge2;
					addFace(e);
				}

				if(((j + 0.5) * PIx2 * jmax) < ((jang - jang0) * jmax + (k + 0.5) * PIx2) * jmax0) {
					edge1 = edgei0 + 2 * j;
					edge2 = edgei0 + 2 * j + 1;
					j++;
				}
				else {
					edge1 = edgei1 + 2 * k;
					edge2 = edgei1 + 2 * k + 1;
					k++;
				}
			}

			if(edge1 != NONE) {
				e[0] = edgei1 + 2 * jmax;
				e[1] = getEdgeSize() - 2;
				e[2] = edge1;
				addFace(e);
			}
			if(edge2 != NONE) {
				e[0] = edgei1 + 2 * jmax + 1;
				e[1] = getEdgeSize() - 1;
				e[2] = edge2;
				addFace(e);
			}
		}
		else {
			for(j=0; j<jmax; j++) {
				ang = jang + 2.0 * j * PI / double(jmax);
				const double jcos = isin * std::cos(ang);
				const double jsin = isin * std::sin(ang);
				addNode(Vector4(icos, jcos, jsin, 0.0));
				if(j > 0) addEdge(getNodeSize() - 1, getNodeSize() - 2);
			}
			addEdge(getNodeSize() - jmax, getNodeSize() - 1);

			uint edge1 = NONE;
			uint edge2 = NONE;
			for(j=0, k=0; k<jmax; ) {
				addEdge(nodei0 + 2 * (j<jmax0?j:j-jmax0), nodei1 + k);
				addEdge(nodei0 + 2 * (j<jmax0?j:j-jmax0) + 1, nodei1 + k);

				if(edge1 != NONE) {
					e[0] = getEdgeSize() - 2;
					e[1] = getEdgeSize() - 4;
					e[2] = edge1;
					addFace(e);
				}
				if(edge2 != NONE) {
					e[0] = getEdgeSize() - 1;
					e[1] = getEdgeSize() - 3;
					e[2] = edge2;
					addFace(e);
				}

				if(((j + 0.5) * PIx2 * jmax) < ((jang - jang0) * jmax + (k + 0.5) * PIx2) * jmax0) {
					edge1 = edgei0 + 2 * j;
					edge2 = edgei0 + 2 * j + 1;
					j++;
				}
				else {
					edge1 = edgei1 + k;
					edge2 = edgei1 + k;
					k++;
				}
			}

			if(edge1 != NONE) {
				e[0] = edgei1 + jmax;
				e[1] = getEdgeSize() - 2;
				e[2] = edge1;
				addFace(e);
			}
			if(edge2 != NONE) {
				e[0] = edgei1 + jmax + 1;
				e[1] = getEdgeSize() - 1;
				e[2] = edge2;
				addFace(e);
			}
		}

		jang0 = jang;
		jmax0 = jmax;
		nodei0 = nodei1;
		edgei0 = edgei1;
		jang += 0.5 * jtmax * h / (jmax * isin);
	}

	// optimize elements
	for(uint k=0; k<optimizationSteps; k++) {
		for(i=0; i<getNodeSize(); i++) {
			improveNodeByHodge(i, true, false);
			setNodePosition(i, r * getNodePosition(i).unit());
		}
	}
}

void BuilderMesh::createSphereBoundaryPentaGrid(const double r, const double h) {
	uint i;

	// create base mesh
	const double icos = std::sqrt(0.2) * r;
	const double isin = std::sqrt(r * r - icos * icos);
	double ang = 0.0;
	for(i=0; i<5; i++) {
		addNode(Vector4( icos, isin * std::cos(ang), isin * std::sin(ang), 0.0));
		ang += 0.2 * PI;
		addNode(Vector4(-icos, isin * std::cos(ang), isin * std::sin(ang), 0.0));
		ang += 0.2 * PI;
	}
	addNode(Vector4(r, 0.0, 0.0, 0.0));
	addNode(Vector4(-r, 0.0, 0.0, 0.0));

	for(i=0; i<5; i++) addEdge(10, 2*i);
	for(i=0; i<5; i++) addEdge(2*i, (2*i+2)%10);
	for(i=0; i<5; i++) addEdge(2*i, 2*i+1);
	for(i=0; i<5; i++) addEdge(2*i+1, (2*i+2)%10);
	for(i=0; i<5; i++) addEdge(2*i+1, (2*i+3)%10);
	for(i=0; i<5; i++) addEdge(11, 2*i+1);

	Buffer<uint> e(3);
	for(i=0; i<5; i++) { e[0] = (i+1)%5; e[1] = i+5; e[2] = i; addFace(e); }
	for(i=0; i<5; i++) { e[0] = i+5; e[1] = i+15; e[2] = i+10; addFace(e); }
	for(i=0; i<5; i++) { e[0] = i+20; e[1] = i+15; e[2] = ((i+1)%5)+10; addFace(e); }
	for(i=0; i<5; i++) { e[0] = i+25; e[1] = i+20; e[2] = ((i+1)%5)+25; addFace(e); }

	//double len = r * std::acos(sqrt(0.2)); 
	double len = r * 1.107148717794090503017; //acos(sqrt(0.2))=1.107148717794090503017...
	while(len > h) {
		Mesh mesh;
		swap(mesh);
		createFaceSplit(mesh, 2);
		for(i=0; i<getNodeSize(); i++) setNodePosition(i, r * getNodePosition(i).unit());
		len *= 0.5;
	}
}

void BuilderMesh::createFaceSplit(const Mesh &mesh, const uint div) {
	clear();

	uint i, j, k;

	// create nodes
	for(i=0; i<mesh.getNodeSize(); i++) {
		addNode(mesh.getNodePosition(i));
	}
	for(i=0; i<mesh.getEdgeSize(); i++) {
		const Buffer<uint> &n = mesh.getEdgeNodes(i);
		const Vector4 p = mesh.getNodePosition(n[0]);
		const Vector4 d = (mesh.getNodePosition(n[1]) - p) / double(div);
		for(j=1; j<div; j++) addNode(p + j * d);
	}
	Buffer<uint> n((div + 1) * (div + 2) / 2);
	const uint i1 = div*(div+1)/2;
	for(i=0; i<mesh.getFaceSize(); i++) {
		const Buffer<uint> &e = mesh.getFaceEdges(i);
		n[0] = mesh.getEdgeIntersection(e[2], e[0]);
		n[i1] = mesh.getEdgeIntersection(e[0], e[1]);
		n[i1 + div] = mesh.getEdgeIntersection(e[1], e[2]);

		uint en = mesh.getNodeSize() + (div - 1) * e[0];
		if(n[0] == mesh.getEdgeNodes(e[0])[0]) {
			for(j=1; j<div; j++) n[j*(j+1)/2] = en + j - 1;
		}
		else {
			for(j=1; j<div; j++) n[j*(j+1)/2] = en + div - j - 1;
		}
		en = mesh.getNodeSize() + (div - 1) * e[1];
		if(n[i1] == mesh.getEdgeNodes(e[1])[0]) {
			for(j=1; j<div; j++) n[i1+j] = en + j - 1;
		}
		else {
			for(j=1; j<div; j++) n[i1+j] = en + div - j - 1;
		}
		en = mesh.getNodeSize() + (div - 1) * e[2];
		if(n[0] == mesh.getEdgeNodes(e[2])[0]) {
			for(j=1; j<div; j++) n[j*(j+3)/2] = en + j - 1;
		}
		else {
			for(j=1; j<div; j++) n[j*(j+3)/2] = en + div - j - 1;
		}

		for(j=2; j<div; j++) {
			const uint j1 = j*(j+1)/2;
			const Vector4 p = getNodePosition(n[j1]);
			const Vector4 d = (getNodePosition(n[j1+j]) - p) / double(j);
			for(k=1; k<j; k++) n[j1+k] = addNode(p + k * d);
		}

		// create edges and faces
		Buffer<uint> ee(3);
		for(j=1; j<=div; j++) {
			for(k=1; k<=j; k++) {
				const uint node = j*(j+1)/2+k;
				ee[0] = addEdge(n[node-1], n[node-j-1]);
				if(k > 1) {
					ee[2] = findEdge(n[node-j-1], n[node-j-2]);
					addFace(ee);
				}
				ee[1] = addEdge(n[node], n[node-j-1]);
				ee[2] = addEdge(n[node-1], n[node]);
				addFace(ee);
			}
		}
	}
}

void BuilderMesh::createQccGrid(const Vector4 &minp, const Vector4 &maxp, const Vector4 &h) {
	const uint xsize = uint((maxp.x - minp.x) / h.x - 1e-8) + 1;
	const uint ysize = uint((maxp.y - minp.y) / h.y - 1e-8) + 1;
	const uint zsize = uint((maxp.z - minp.z) / h.z - 1e-8) + 1;
	const uint tsize = uint((maxp.t - minp.t) / h.t - 1e-8) + 1;

	const uint xs = (xsize < 3 ? xsize : 3);
	const uint ys = (ysize < 3 ? ysize : 3);
	const uint zs = (zsize < 3 ? zsize : 3);
	const uint ts = (tsize < 3 ? tsize : 3);

	const Vector4 p0 = 0.5 * (maxp + minp - Vector4(xsize * h.x, ysize * h.y, zsize * h.z, tsize * h.t));
	const Vector4 p1 = p0 + Vector4(xs * h.x, ys * h.y, zs * h.z, ts * h.t);

	createGrid(p0, p1, h);

	uint node = 0;
	uint xi, yi, zi, ti;
	for(ti=0; ti<ts; ti++) {
		const double pt = p0.t + h.t * (ti + 0.5);
		for(zi=0; zi<zs; zi++) {
			const double pz = p0.z + h.z * (zi + 0.5);
			for(yi=0; yi<ys; yi++) {
				const double py = p0.y + h.y * (yi + 0.5);
				for(xi=0; xi<xs; xi++) {
					const double px = p0.x + h.x * (xi + 0.5);
					node = insertNode(Vector4(px, py, pz, pt), 0.0, node, true);
				}
			}
		}
	}

	const Vector4 pp(p0 + 0.99999 * h);
	if(xsize > xs) repeatMiddle(pp, Vector4(h.x,0,0,0), xsize - xs);
	if(ysize > ys) repeatMiddle(pp, Vector4(0,h.y,0,0), ysize - ys);
	if(zsize > zs) repeatMiddle(pp, Vector4(0,0,h.z,0), zsize - zs);
	if(tsize > ts) repeatMiddle(pp, Vector4(0,0,0,h.t), tsize - ts);
}

void BuilderMesh::createAsyncQccGrid(const Vector4 &minp, const Vector4 &maxp, const Vector4 &h) {
	const uint xsize = uint((maxp.x - minp.x) / h.x - 1e-8) + 1;
	const uint ysize = uint((maxp.y - minp.y) / h.y - 1e-8) + 1;
	const uint zsize = uint((maxp.z - minp.z) / h.z - 1e-8) + 1;
	const uint tsize = uint((maxp.t - minp.t) / h.t - 1e-8) + 1;

	const uint xs = (xsize < 3 ? xsize : 3);
	const uint ys = (ysize < 3 ? ysize : 3);
	const uint zs = (zsize < 3 ? zsize : 3);
	const uint ts = (tsize < 3 ? tsize : 3);

	const Vector4 p0 = 0.5 * (maxp + minp - Vector4(xsize * h.x, ysize * h.y, zsize * h.z, tsize * h.t));
	const Vector4 p1 = p0 + Vector4(xs * h.x, ys * h.y, zs * h.z, ts * h.t);

	createGrid(p0, p1, h);

	uint node = 0;
	uint xi, yi, zi, ti;
	for(ti=0; ti<=ts; ti++) {
		const double pt = p0.t + h.t * ti;
		for(zi=0; zi<=zs; zi++) {
			const double pz = p0.z + h.z * zi;
			for(yi=0; yi<=ys; yi++) {
				const double py = p0.y + h.y * yi;
				for(xi=0; xi<=xs; xi++) {
					const double px = p0.x + h.x * xi;
					if(xi<xs && yi<ys && zi<zs && ti<ts) {
						node = insertNode(Vector4(px + 0.5 * h.x, py + 0.5 * h.y, pz + 0.5 * h.z, pt + 0.5 * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.25 * h.x, py + 0.25 * h.y, pz + 0.25 * h.z, pt + 0.25 * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.75 * h.x, py + 0.25 * h.y, pz + 0.25 * h.z, pt + 0.75 * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.25 * h.x, py + 0.75 * h.y, pz + 0.25 * h.z, pt + 0.75 * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.25 * h.x, py + 0.25 * h.y, pz + 0.75 * h.z, pt + 0.75 * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.25 * h.x, py + 0.75 * h.y, pz + 0.75 * h.z, pt + 0.25 * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.75 * h.x, py + 0.25 * h.y, pz + 0.75 * h.z, pt + 0.25 * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.75 * h.x, py + 0.75 * h.y, pz + 0.25 * h.z, pt + 0.25 * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.75 * h.x, py + 0.75 * h.y, pz + 0.75 * h.z, pt + 0.75 * h.t), 0.0, node, true);
					}

					if(xi<xs && ti<ts) node = insertNode(Vector4(px + 0.5 * h.x, py, pz, pt + 0.5 * h.t), 0.0, node, true);
					if(yi<ys && ti<ts) node = insertNode(Vector4(px, py + 0.5 * h.y, pz, pt + 0.5 * h.t), 0.0, node, true);
					if(zi<zs && ti<ts) node = insertNode(Vector4(px, py, pz + 0.5 * h.z, pt + 0.5 * h.t), 0.0, node, true);
					if(yi<ys && zi<zs) node = insertNode(Vector4(px, py + 0.5 * h.y, pz + 0.5 * h.z, pt), 0.0, node, true);
					if(xi<xs && zi<zs) node = insertNode(Vector4(px + 0.5 * h.x, py, pz + 0.5 * h.z, pt), 0.0, node, true);
					if(xi<xs && yi<ys) node = insertNode(Vector4(px + 0.5 * h.x, py + 0.5 * h.y, pz, pt), 0.0, node, true);
				}
			}
		}
	}

	const Vector4 pp(p0 + 0.99999 * h);
	if(xsize > xs) repeatMiddle(pp, Vector4(h.x,0,0,0), xsize - xs);
	if(ysize > ys) repeatMiddle(pp, Vector4(0,h.y,0,0), ysize - ys);
	if(zsize > zs) repeatMiddle(pp, Vector4(0,0,h.z,0), zsize - zs);
	if(tsize > ts) repeatMiddle(pp, Vector4(0,0,0,h.t), tsize - ts);
}

void BuilderMesh::createAsyncQccGrid2(const Vector4 &minp, const Vector4 &maxp, const Vector4 &h) {
	const uint xsize = uint((maxp.x - minp.x) / h.x - 1e-8) + 1;
	const uint ysize = uint((maxp.y - minp.y) / h.y - 1e-8) + 1;
	const uint zsize = uint((maxp.z - minp.z) / h.z - 1e-8) + 1;
	const uint tsize = uint((maxp.t - minp.t) / h.t - 1e-8) + 1;

	const uint xs = (xsize < 3 ? xsize : 3);
	const uint ys = (ysize < 3 ? ysize : 3);
	const uint zs = (zsize < 3 ? zsize : 3);
	const uint ts = (tsize < 3 ? tsize : 3);

	const Vector4 p0 = 0.5 * (maxp + minp - Vector4(xsize * h.x, ysize * h.y, zsize * h.z, tsize * h.t));
	const Vector4 p1 = p0 + Vector4(xs * h.x, ys * h.y, zs * h.z, ts * h.t);

	createGrid(p0, p1, h);

	const double twist = 0.25; //(std::sqrt(3.0) - 1.0) / (2.0 + std::sqrt(3.0));//1.0 / (2.0 + std::sqrt(3.0)); // or 0.196;
	uint node = 0;
	uint xi, yi, zi, ti;
	for(ti=0; ti<=ts; ti++) {
		const double pt = p0.t + h.t * ti;
		for(zi=0; zi<=zs; zi++) {
			const double pz = p0.z + h.z * zi;
			for(yi=0; yi<=ys; yi++) {
				const double py = p0.y + h.y * yi;
				for(xi=0; xi<=xs; xi++) {
					const double px = p0.x + h.x * xi;
					if(xi<xs && yi<ys && zi<zs && ti<ts) {
						node = insertNode(Vector4(px + 0.5 * h.x, py + 0.5 * h.y, pz + 0.5 * h.z, pt + twist * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.25 * h.x, py + 0.25 * h.y, pz + 0.25 * h.z, pt + 0.5 * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.75 * h.x, py + 0.25 * h.y, pz + 0.25 * h.z, pt + (0.5 + twist) * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.25 * h.x, py + 0.75 * h.y, pz + 0.25 * h.z, pt + (0.5 + twist) * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.25 * h.x, py + 0.25 * h.y, pz + 0.75 * h.z, pt + (0.5 + twist) * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.25 * h.x, py + 0.75 * h.y, pz + 0.75 * h.z, pt + 0.5 * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.75 * h.x, py + 0.25 * h.y, pz + 0.75 * h.z, pt + 0.5 * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.75 * h.x, py + 0.75 * h.y, pz + 0.25 * h.z, pt + 0.5 * h.t), 0.0, node, true);
						node = insertNode(Vector4(px + 0.75 * h.x, py + 0.75 * h.y, pz + 0.75 * h.z, pt + (0.5 + twist) * h.t), 0.0, node, true);
					}

					if(xi<xs && ti<ts) node = insertNode(Vector4(px + 0.5 * h.x, py, pz, pt + twist * h.t), 0.0, node, true);
					if(yi<ys && ti<ts) node = insertNode(Vector4(px, py + 0.5 * h.y, pz, pt + twist * h.t), 0.0, node, true);
					if(zi<zs && ti<ts) node = insertNode(Vector4(px, py, pz + 0.5 * h.z, pt + twist * h.t), 0.0, node, true);
					if(yi<ys && zi<zs) node = insertNode(Vector4(px, py + 0.5 * h.y, pz + 0.5 * h.z, pt), 0.0, node, true);
					if(xi<xs && zi<zs) node = insertNode(Vector4(px + 0.5 * h.x, py, pz + 0.5 * h.z, pt), 0.0, node, true);
					if(xi<xs && yi<ys) node = insertNode(Vector4(px + 0.5 * h.x, py + 0.5 * h.y, pz, pt), 0.0, node, true);
				}
			}
		}
	}

	const Vector4 pp(p0 + 0.99999 * h);
	if(xsize > xs) repeatMiddle(pp, Vector4(h.x,0,0,0), xsize - xs);
	if(ysize > ys) repeatMiddle(pp, Vector4(0,h.y,0,0), ysize - ys);
	if(zsize > zs) repeatMiddle(pp, Vector4(0,0,h.z,0), zsize - zs);
	if(tsize > ts) repeatMiddle(pp, Vector4(0,0,0,h.t), tsize - ts);
}

bool BuilderMesh::createRepeated(const Mesh &mesh, const Vector4 &pos, const Vector4 &dir, const Vector4 &step, const uint steps) {
	uint i, j, k;
	const double dirsq = std::abs(dir.dot(step));

	// find a slot for each existing cell
	const uint qsize = mesh.getQuadSize();
	Buffer<uint> qslot(qsize);
	Buffer<uint> qind(qsize);
	Buffer<uint> qs(4, 0);
	for(i=0; i<qsize; i++) {
		const Vector4 p = mesh.getQuadPosition(i) - pos;
		const double dot = p.dot(dir);
		if(dot < 0.0) qslot[i] = 0;
		else if(dot < dirsq) qslot[i] = 2;
		else qslot[i] = 3;
		qind[i] = qs[qslot[i]]++;
	}
	const uint bsize = mesh.getBodySize();
	Buffer<uint> bslot(bsize);
	Buffer<uint> bind(bsize);
	Buffer<uint> bs(4, 0);
	for(i=0; i<bsize; i++) {
		const Buffer<uint> &par = mesh.getBodyQuads(i);
		if(par.empty()) {
			const Vector4 p = mesh.getBodyPosition(i) - pos;
			const double dot = p.dot(dir);
			if(dot < 0.0) bslot[i] = 0;
			else if(dot < dirsq) bslot[i] = 2;
			else bslot[i] = 3;
		}
		else {
			bslot[i] = qslot[par[0]];
			uint maxslot = bslot[i];
			for(j=1; j<par.size(); j++) {
				const uint jslot = qslot[par[j]];
				if(jslot < bslot[i]) bslot[i] = jslot;
				else if(jslot > maxslot) maxslot = jslot;
			}
			if(bslot[i] == 0 && maxslot != 0) bslot[i] = 1;
		}
		bind[i] = bs[bslot[i]]++;
	}
	const uint fsize = mesh.getFaceSize();
	Buffer<uint> fslot(fsize);
	Buffer<uint> find(fsize);
	Buffer<uint> fs(4, 0);
	for(i=0; i<fsize; i++) {
		const Buffer<uint> &par = mesh.getFaceBodies(i);
		if(par.empty()) {
			const Vector4 p = mesh.getFacePosition(i) - pos;
			const double dot = p.dot(dir);
			if(dot < 0.0) fslot[i] = 0;
			else if(dot < dirsq) fslot[i] = 2;
			else fslot[i] = 3;
		}
		else {
			fslot[i] = bslot[par[0]];
			uint maxslot = fslot[i];
			for(j=1; j<par.size(); j++) {
				const uint jslot = bslot[par[j]];
				if(jslot < fslot[i]) fslot[i] = jslot;
				else if(jslot > maxslot) maxslot = jslot;
			}
			if(fslot[i] == 0 && maxslot != 0) fslot[i] = 1;
		}
		find[i] = fs[fslot[i]]++;
	}
	const uint esize = mesh.getEdgeSize();
	Buffer<uint> eslot(esize);
	Buffer<uint> eind(esize);
	Buffer<uint> es(4, 0);
	for(i=0; i<esize; i++) {
		const Buffer<uint> &par = mesh.getEdgeFaces(i);
		if(par.empty()) {
			const Vector4 p = mesh.getEdgePosition(i) - pos;
			const double dot = p.dot(dir);
			if(dot < 0.0) eslot[i] = 0;
			else if(dot < dirsq) eslot[i] = 2;
			else eslot[i] = 3;
		}
		else {
			eslot[i] = fslot[par[0]];
			uint maxslot = eslot[i];
			for(j=1; j<par.size(); j++) {
				const uint jslot = fslot[par[j]];
				if(jslot < eslot[i]) eslot[i] = jslot;
				else if(jslot > maxslot) maxslot = jslot;
			}
			if(eslot[i] == 0 && maxslot != 0) eslot[i] = 1;
		}
		eind[i] = es[eslot[i]]++;
	}
	const uint nsize = mesh.getNodeSize();
	Buffer<uint> nslot(nsize);
	Buffer<uint> nind(nsize);
	Buffer<uint> ns(4, 0);
	for(i=0; i<nsize; i++) {
		const Buffer<uint> &par = mesh.getNodeEdges(i);
		if(par.empty()) {
			const Vector4 p = getNodePosition(i) - pos;
			const double dot = p.dot(dir);
			if(dot < 0.0) nslot[i] = 0;
			else if(dot < dirsq) nslot[i] = 2;
			else nslot[i] = 3;
		}
		else {
			nslot[i] = eslot[par[0]];
			uint maxslot = nslot[i];
			for(j=1; j<par.size(); j++) {
				const uint jslot = eslot[par[j]];
				if(jslot < nslot[i]) nslot[i] = jslot;
				else if(jslot > maxslot) maxslot = jslot;
			}
			if(nslot[i] == 0 && maxslot != 0) nslot[i] = 1;
		}
		nind[i] = ns[nslot[i]]++;
	}

	// compute links for repetition
	uint cell = 0;
	Buffer<uint> nlink(ns[1]);
	for(i=0; i<nsize; i++) {
		if(nslot[i] != 1) continue;
		cell = mesh.findNode(mesh.getNodePosition(i) + step, 1e-13, cell);
		if(cell == NONE) return false;
		nlink[nind[i]] = cell;
	}
	Buffer<uint> elink(es[1]);
	for(i=0; i<esize; i++) {
		if(eslot[i] != 1) continue;
		const Buffer<uint> bou = mesh.getEdgeNodes(i);
		for(j=0; j<bou.size(); j++) bou[j] = nlink[nind[bou[j]]];
		cell = mesh.findEdge(bou[0], bou[1]);
		if(cell == NONE) return false;
		elink[eind[i]] = cell;
	}
	Buffer<uint> flink(fs[1]);
	for(i=0; i<fsize; i++) {
		if(fslot[i] != 1) continue;
		const Buffer<uint> bou = mesh.getFaceEdges(i);
		for(j=0; j<bou.size(); j++) bou[j] = elink[eind[bou[j]]];
		cell = mesh.findFace(bou);
		if(cell == NONE) return false;
		flink[find[i]] = cell;
	}
	Buffer<uint> blink(bs[1]);
	for(i=0; i<bsize; i++) {
		if(bslot[i] != 1) continue;
		const Buffer<uint> bou = mesh.getBodyFaces(i);
		for(j=0; j<bou.size(); j++) bou[j] = flink[find[bou[j]]];
		cell = mesh.findBody(bou);
		if(cell == NONE) return false;
		blink[bind[i]] = cell;
	}
	Buffer<uint> qlink(qs[1]); // is empty

	// create mesh
	clear();

	// create nodes
	for(i=mesh.getMetricSize(); i-->0; ) setMetric(mesh.getMetric(i), i);
	m_nsize = nsize + steps * ns[2];
	resizeNodeBuffer(m_nsize);
	for(i=0; i<nsize; i++) {
		uint ii = getIndex(i, 0, nind, nslot, nlink, ns);
		if(nslot[i] <= 1) {
			setNodePosition(ii, mesh.getNodePosition(i));
			setNodeWeight(ii, mesh.getNodeWeight(i));
			setNodeFlag(ii, mesh.getNodeFlag(i));
		}
		else if(nslot[i] == 2) {
			for(j=0; j<=steps; j++) {
				setNodePosition(ii, mesh.getNodePosition(i) + j * step);
				setNodeWeight(ii, mesh.getNodeWeight(i));
				setNodeFlag(ii, mesh.getNodeFlag(i));
				ii += ns[2];
			}
		}
		else {
			ii += steps * ns[2];
			setNodePosition(ii, mesh.getNodePosition(i) + steps * step);
			setNodeWeight(ii, mesh.getNodeWeight(i));
			setNodeFlag(ii, mesh.getNodeFlag(i));
		}
	}

	// create edges
	m_esize = esize + steps * es[2];
	resizeEdgeBuffer(m_esize);
	for(i=0; i<esize; i++) {
		uint ii = getIndex(i, 0, eind, eslot, elink, es);
		if(eslot[i] <= 1) {
			Buffer<uint> &bou = m_e[ii].n;
			bou = mesh.getEdgeNodes(i);
			for(k=0; k<bou.size(); k++) {
				bou[k] = getIndex(bou[k], 0, nind, nslot, nlink, ns);
				m_n[bou[k]].e.push_back(ii);
			}
			setEdgeFlag(ii, mesh.getEdgeFlag(i));
		}
		else if(eslot[i] == 2) {
			for(j=0; j<=steps; j++) {
				Buffer<uint> &bou = m_e[ii].n;
				bou = mesh.getEdgeNodes(i);
				for(k=0; k<bou.size(); k++) {
					bou[k] = getIndex(bou[k], j, nind, nslot, nlink, ns);
					m_n[bou[k]].e.push_back(ii);
				}
				setEdgeFlag(ii, mesh.getEdgeFlag(i));
				ii += es[2];
			}
		}
		else {
			ii += steps * es[2];
			Buffer<uint> &bou = m_e[ii].n;
			bou = mesh.getEdgeNodes(i);
			for(k=0; k<bou.size(); k++) {
				bou[k] = getIndex(bou[k], steps, nind, nslot, nlink, ns);
				m_n[bou[k]].e.push_back(ii);
			}
			setEdgeFlag(ii, mesh.getEdgeFlag(i));
		}
	}
	// create faces
	m_fsize = fsize + steps * fs[2];
	resizeFaceBuffer(m_fsize);
	for(i=0; i<fsize; i++) {
		uint ii = getIndex(i, 0, find, fslot, flink, fs);
		if(fslot[i] <= 1) {
			Buffer<uint> &bou = m_f[ii].e;
			bou = mesh.getFaceEdges(i);
			for(k=0; k<bou.size(); k++) {
				bou[k] = getIndex(bou[k], 0, eind, eslot, elink, es);
				m_e[bou[k]].f.push_back(ii);
			}
			setFaceFlag(ii, mesh.getFaceFlag(i));
		}
		else if(fslot[i] == 2) {
			for(j=0; j<=steps; j++) {
				Buffer<uint> &bou = m_f[ii].e;
				bou = mesh.getFaceEdges(i);
				for(k=0; k<bou.size(); k++) {
					bou[k] = getIndex(bou[k], j, eind, eslot, elink, es);
					m_e[bou[k]].f.push_back(ii);
				}
				setFaceFlag(ii, mesh.getFaceFlag(i));
				ii += fs[2];
			}
		}
		else {
			ii += steps * fs[2];
			Buffer<uint> &bou = m_f[ii].e;
			bou = mesh.getFaceEdges(i);
			for(k=0; k<bou.size(); k++) {
				bou[k] = getIndex(bou[k], steps, eind, eslot, elink, es);
				m_e[bou[k]].f.push_back(ii);
			}
			setFaceFlag(ii, mesh.getFaceFlag(i));
		}
	}
	// create bodies
	m_bsize = bsize + steps * bs[2];
	resizeBodyBuffer(m_bsize);
	for(i=0; i<bsize; i++) {
		uint ii = getIndex(i, 0, bind, bslot, blink, bs);
		if(bslot[i] <= 1) {
			Buffer<uint> &bou = m_b[ii].f;
			bou = mesh.getBodyFaces(i);
			for(k=0; k<bou.size(); k++) {
				bou[k] = getIndex(bou[k], 0, find, fslot, flink, fs);
				m_f[bou[k]].b.push_back(ii);
			}
			setBodyFlag(ii, mesh.getBodyFlag(i));
		}
		else if(bslot[i] == 2) {
			for(j=0; j<=steps; j++) {
				Buffer<uint> &bou = m_b[ii].f;
				bou = mesh.getBodyFaces(i);
				for(k=0; k<bou.size(); k++) {
					bou[k] = getIndex(bou[k], j, find, fslot, flink, fs);
					m_f[bou[k]].b.push_back(ii);
				}
				setBodyFlag(ii, mesh.getBodyFlag(i));
				ii += bs[2];
			}
		}
		else {
			ii += steps * bs[2];
			Buffer<uint> &bou = m_b[ii].f;
			bou = mesh.getBodyFaces(i);
			for(k=0; k<bou.size(); k++) {
				bou[k] = getIndex(bou[k], steps, find, fslot, flink, fs);
				m_f[bou[k]].b.push_back(ii);
			}
			setBodyFlag(ii, mesh.getBodyFlag(i));
		}
	}
	// create quads
	m_qsize = qsize + steps * qs[2];
	resizeQuadBuffer(m_qsize);
	for(i=0; i<qsize; i++) {
		uint ii = getIndex(i, 0, qind, qslot, qlink, qs);
		if(qslot[i] <= 1) {
			Buffer<uint> &bou = m_q[ii].b;
			bou = mesh.getQuadBodies(i);
			for(k=0; k<bou.size(); k++) {
				bou[k] = getIndex(bou[k], 0, bind, bslot, blink, bs);
				m_b[bou[k]].q.push_back(ii);
			}
			setQuadFlag(ii, mesh.getQuadFlag(i));
		}
		else if(qslot[i] == 2) {
			for(j=0; j<=steps; j++) {
				Buffer<uint> &bou = m_q[ii].b;
				bou = mesh.getQuadBodies(i);
				for(k=0; k<bou.size(); k++) {
					bou[k] = getIndex(bou[k], j, bind, bslot, blink, bs);
					m_b[bou[k]].q.push_back(ii);
				}
				setQuadFlag(ii, mesh.getQuadFlag(i));
				ii += qs[2];
			}
		}
		else {
			ii += steps * qs[2];
			Buffer<uint> &bou = m_q[ii].b;
			bou = mesh.getQuadBodies(i);
			for(k=0; k<bou.size(); k++) {
				bou[k] = getIndex(bou[k], steps, bind, bslot, blink, bs);
				m_b[bou[k]].q.push_back(ii);
			}
			setQuadFlag(ii, mesh.getQuadFlag(i));
		}
	}
	return true;
}

uint BuilderMesh::getIndex(uint i, uint num, const Buffer<uint> &ind, const Buffer<uint> &slot, const Buffer<uint> &link, const Buffer<uint> &s) {
	if(slot[i] == 0) return ind[i];
	while(slot[i] == 1) {
		if(num == 0) return s[0] + ind[i];
		i = link[ind[i]];
		num--;
	}
	if(slot[i] == 2) return s[0] + s[1] + s[2] * num + ind[i];
	return s[0] + s[1] + s[2] * (num + 1) + ind[i];
}

void BuilderMesh::createIntersection(const Mesh &mesh, const Vector4 &v, const double dot) {
	uint i, j;
	clear();

	// calculate side (of the plane) for each node
	Buffer<bool> nside(mesh.getNodeSize());
	for(i=0; i<nside.size(); i++) {
		nside[i] = (v.dot(mesh.getNodePosition(i)) > dot);
	}

	// find edges which are intersect by the plane
	uint esize = 0;
	Buffer<uint> e(mesh.getEdgeSize());
	Buffer<uint> enew(e.size(), NONE);
	for(i=0; i<e.size(); i++) {
		const Buffer<uint> &en = mesh.getEdgeNodes(i);
		if(nside[en[0]] != nside[en[1]]) {
			enew[i] = esize;
			e[esize++] = i;
		}
	}

	// add nodes at intersections of each edge and the plane
	resizeNodeBuffer(esize);
	for(i=0; i<esize; i++) {
		const Buffer<uint> &en = mesh.getEdgeNodes(e[i]);
		Vector4 p = mesh.getNodePosition(en[0]);
		const Vector4 d = mesh.getNodePosition(en[1]) - p;
		const double pdotv = dot - p.dot(v);
		const double ddotv = d.dot(v);
		if(pdotv * ddotv >= ddotv * ddotv) p += d;
		else if(pdotv * ddotv > 0.0) p += (dot - p.dot(v)) / ddotv * d;
		setNodeFlag(addNode(p), mesh.getEdgeFlag(e[i]));
	}

	// find faces which are intersect by the plane
	uint fsize = 0;
	Buffer<uint> f(mesh.getFaceSize());
	Buffer<uint> fnew(f.size(), NONE);
	for(i=0; i<f.size(); i++) {
		const Buffer<uint> &fe = mesh.getFaceEdges(i);
		for(j=0; j<fe.size(); j++) {
			if(enew[fe[j]] != NONE) {
				fnew[i] = fsize;
				f[fsize++] = i;
				break;
			}
		}
	}

	// add edges
	resizeEdgeBuffer(fsize);
	for(i=0; i<fsize; i++) {
		const Buffer<uint> &fe = mesh.getFaceEdges(f[i]);
		uint ens = 0;
		Buffer<uint> en(fe.size());
		for(j=0; j<fe.size(); j++) {
			if(enew[fe[j]] != NONE) en[ens++] = enew[fe[j]];
		}
		setEdgeFlag(addEdge(en[0], en[1]), mesh.getFaceFlag(f[i]));
	}

	// find bodies which are intersect by the plane
	uint bsize = 0;
	Buffer<uint> b(mesh.getBodySize());
	Buffer<uint> bnew(b.size(), NONE);
	for(i=0; i<b.size(); i++) {
		const Buffer<uint> &bf = mesh.getBodyFaces(i);
		for(j=0; j<bf.size(); j++) {
			if(fnew[bf[j]] != NONE) {
				bnew[i] = bsize;
				b[bsize++] = i;
				break;
			}
		}
	}

	// add faces
	resizeFaceBuffer(bsize);
	for(i=0; i<bsize; i++) {
		const Buffer<uint> &bf = mesh.getBodyFaces(b[i]);
		uint fes = 0;
		Buffer<uint> fe(bf.size());
		for(j=0; j<bf.size(); j++) {
			if(fnew[bf[j]] != NONE) fe[fes++] = fnew[bf[j]];
		}
		fe.resize(fes);
		setFaceFlag(addFace(fe), mesh.getBodyFlag(b[i]));
	}

	// find quads which are intersect by the plane
	uint qsize = 0;
	Buffer<uint> q(mesh.getQuadSize());
	for(i=0; i<q.size(); i++) {
		const Buffer<uint> &qb = mesh.getQuadBodies(i);
		for(j=0; j<qb.size(); j++) {
			if(bnew[qb[j]] != NONE) {
				q[qsize++] = i;
				break;
			}
		}
	}

	// add bodies
	resizeBodyBuffer(qsize);
	for(i=0; i<qsize; i++) {
		const Buffer<uint> &qb = mesh.getQuadBodies(q[i]);
		uint bfs = 0;
		Buffer<uint> bf(qb.size());
		for(j=0; j<qb.size(); j++) {
			if(bnew[qb[j]] != NONE) bf[bfs++] = bnew[qb[j]];
		}
		bf.resize(bfs);
		setBodyFlag(addBody(bf), mesh.getQuadFlag(q[i]));
	}
}

void BuilderMesh::createDualMesh(const Mesh &mesh) {
	clear();
	if(mesh.getEdgeSize() == 0) {
		createCopy(mesh);
		return;
	}
	uint i, j;
	if(mesh.getFaceSize() == 0) {
		for(i=0; i<mesh.getEdgeSize(); i++) {
			addNode(mesh.getEdgePosition(i));
		}
		for(i=0; i<mesh.getNodeSize(); i++) {
			const Buffer<uint> &n = mesh.getNodeEdges(i);
			if(n.size() != 2) continue;
			addEdge(n[0], n[1]);
		}
		return;
	}
	if(mesh.getBodySize() == 0) {
		for(i=0; i<mesh.getFaceSize(); i++) {
			addNode(mesh.getFacePosition(i));
		}
		Buffer<uint> me(mesh.getEdgeSize(), NONE);
		for(i=0; i<me.size(); i++) {
			const Buffer<uint> &n = mesh.getEdgeFaces(i);
			if(n.size() != 2) continue;
			me[i] = addEdge(n[0], n[1]);
		}
		for(i=0; i<mesh.getNodeSize(); i++) {
			Buffer<uint> e = mesh.getNodeEdges(i);
			if(e.size() < 3) continue;
			for(j=0; j<e.size(); j++) {
				if(me[e[j]] == NONE) break;
				e[j] = me[e[j]];
			}
			if(j < e.size()) continue;
			addFace(e);
		}
		return;
	}
	if(mesh.getQuadSize() == 0) {
		for(i=0; i<mesh.getBodySize(); i++) {
			addNode(mesh.getBodyPosition(i));
		}
		Buffer<uint> me(mesh.getFaceSize(), NONE);
		for(i=0; i<me.size(); i++) {
			const Buffer<uint> &n = mesh.getFaceBodies(i);
			if(n.size() != 2) continue;
			me[i] = addEdge(n[0], n[1]);
		}
		Buffer<uint> mf(mesh.getEdgeSize(), NONE);
		for(i=0; i<mf.size(); i++) {
			Buffer<uint> e = mesh.getEdgeFaces(i);
			if(e.size() < 3) continue;
			for(j=0; j<e.size(); j++) {
				if(me[e[j]] == NONE) break;
				e[j] = me[e[j]];
			}
			if(j < e.size()) continue;
			mf[i] = addFace(e);
		}
		for(i=0; i<mesh.getNodeSize(); i++) {
			Buffer<uint> f = mesh.getNodeEdges(i);
			if(f.size() < 4) continue;
			for(j=0; j<f.size(); j++) {
				if(mf[f[j]] == NONE) break;
				f[j] = mf[f[j]];
			}
			if(j < f.size()) continue;
			addBody(f);
		}
		return;
	}
}

void BuilderMesh::clearFlags() {
	clearNodeFlags();
	clearEdgeFlags();
	clearFaceFlags();
	clearBodyFlags();
	clearQuadFlags();
}
void BuilderMesh::fillNodeFlags(const uint flag, const UintSet &oldflag) {
	for(uint i=0; i<m_nsize; i++) {
		if(oldflag.includes(getNodeFlag(i))) setNodeFlag(i, flag);
	}
}
void BuilderMesh::fillEdgeFlags(const uint flag, const UintSet &oldflag) {
	for(uint i=0; i<m_esize; i++) {
		if(oldflag.includes(getEdgeFlag(i))) setEdgeFlag(i, flag);
	}
}
void BuilderMesh::fillFaceFlags(const uint flag, const UintSet &oldflag) {
	for(uint i=0; i<m_fsize; i++) {
		if(oldflag.includes(getFaceFlag(i))) setFaceFlag(i, flag);
	}
}
void BuilderMesh::fillBodyFlags(const uint flag, const UintSet &oldflag) {
	for(uint i=0; i<m_bsize; i++) {
		if(oldflag.includes(getBodyFlag(i))) setBodyFlag(i, flag);
	}
}
void BuilderMesh::fillQuadFlags(const uint flag, const UintSet &oldflag) {
	for(uint i=0; i<m_qsize; i++) {
		if(oldflag.includes(getQuadFlag(i))) setQuadFlag(i, flag);
	}
}
void BuilderMesh::fillFlags(const uint flag, const UintSet &oldflag) {
	fillNodeFlags(flag, oldflag);
	fillEdgeFlags(flag, oldflag);
	fillFaceFlags(flag, oldflag);
	fillBodyFlags(flag, oldflag);
	fillQuadFlags(flag, oldflag);
}

void BuilderMesh::fillBoundaryFlags(const uint flag, const UintSet &oldflag) {
	uint i, j, k, l;
	if(m_fsize == 0) { // 1d
		for(i=0; i<m_nsize; i++) {
			if(!oldflag.includes(getNodeFlag(i))) continue;
			if(getNodeEdges(i).size() >= 2) continue;
			setNodeFlag(i, flag);
		}
		return;
	}
	if(m_bsize == 0) { // 2d
		for(i=0; i<m_esize; i++) {
			if(!oldflag.includes(getEdgeFlag(i))) continue;
			if(getEdgeFaces(i).size() >= 2) continue;
			setEdgeFlag(i, flag);
			const Buffer<uint> &n = getEdgeNodes(i);
			for(j=0; j<n.size(); j++)
			{
				if(oldflag.includes(getNodeFlag(n[j]))) setNodeFlag(n[j], flag);
			}
		}
		return;
	}
	if(m_qsize == 0) { // 3d
		for(i=0; i<m_fsize; i++) {
			if(!oldflag.includes(getFaceFlag(i))) continue;
			if(getFaceBodies(i).size() >= 2) continue;
			setFaceFlag(i, flag);
			const Buffer<uint> &e = getFaceEdges(i);
			for(j=0; j<e.size(); j++) {
				if(!oldflag.includes(getEdgeFlag(e[j]))) continue;
				setEdgeFlag(e[j], flag);
				const Buffer<uint> &n = getEdgeNodes(e[j]);
				for(k=0; k<n.size(); k++) {
					if(oldflag.includes(getNodeFlag(n[k]))) setNodeFlag(n[k], flag);
				}
			}
		}
		return;
	}
	// 4d
	for(i=0; i<m_bsize; i++) {
		if(!oldflag.includes(getBodyFlag(i))) continue;
		if(getBodyQuads(i).size() >= 2) continue;
		setBodyFlag(i, flag);
		const Buffer<uint> &f = getBodyFaces(i);
		for(j=0; j<f.size(); j++) {
			if(!oldflag.includes(getFaceFlag(f[j]))) continue;
			setFaceFlag(f[j], flag);
			const Buffer<uint> &e = getFaceEdges(f[j]);
			for(k=0; k<e.size(); k++) {
				if(!oldflag.includes(getEdgeFlag(e[k]))) continue;
				setEdgeFlag(e[k], flag);
				const Buffer<uint> &n = getEdgeNodes(e[k]);
				for(l=0; l<n.size(); l++) {
					if(oldflag.includes(getNodeFlag(n[l]))) setNodeFlag(n[l], flag);
				}
			}
		}
	}
}

void BuilderMesh::fillNodeRectangleFlags(const Vector4 &minp, const Vector4 &maxp, const uint flag, const UintSet &oldflag) {
	for(uint i=0; i<m_nsize; i++) {
		if(!oldflag.includes(getNodeFlag(i))) continue;
		const Vector4 p = getNodePosition(i);
		if(p.x < minp.x || maxp.x <= p.x) continue;
		if(p.y < minp.y || maxp.y <= p.y) continue;
		if(p.z < minp.z || maxp.z <= p.z) continue;
		if(p.t < minp.t || maxp.t <= p.t) continue;
		setNodeFlag(i, flag);
	}
}
void BuilderMesh::fillEdgeRectangleFlags(const Vector4 &minp, const Vector4 &maxp, const uint flag, const UintSet &oldflag) {
	for(uint i=0; i<m_esize; i++) {
		if(!oldflag.includes(getEdgeFlag(i))) continue;
		const Vector4 p = getEdgeAverage(i);
		if(p.x < minp.x || maxp.x <= p.x) continue;
		if(p.y < minp.y || maxp.y <= p.y) continue;
		if(p.z < minp.z || maxp.z <= p.z) continue;
		if(p.t < minp.t || maxp.t <= p.t) continue;
		setEdgeFlag(i, flag);
	}
}
void BuilderMesh::fillFaceRectangleFlags(const Vector4 &minp, const Vector4 &maxp, const uint flag, const UintSet &oldflag) {
	for(uint i=0; i<m_fsize; i++) {
		if(!oldflag.includes(getFaceFlag(i))) continue;
		const Vector4 p = getFaceAverage(i);
		if(p.x < minp.x || maxp.x <= p.x) continue;
		if(p.y < minp.y || maxp.y <= p.y) continue;
		if(p.z < minp.z || maxp.z <= p.z) continue;
		if(p.t < minp.t || maxp.t <= p.t) continue;
		setFaceFlag(i, flag);
	}
}
void BuilderMesh::fillBodyRectangleFlags(const Vector4 &minp, const Vector4 &maxp, const uint flag, const UintSet &oldflag) {
	for(uint i=0; i<m_bsize; i++) {
		if(!oldflag.includes(getBodyFlag(i))) continue;
		const Vector4 p = getBodyAverage(i);
		if(p.x < minp.x || maxp.x <= p.x) continue;
		if(p.y < minp.y || maxp.y <= p.y) continue;
		if(p.z < minp.z || maxp.z <= p.z) continue;
		if(p.t < minp.t || maxp.t <= p.t) continue;
		setBodyFlag(i, flag);
	}
}
void BuilderMesh::fillQuadRectangleFlags(const Vector4 &minp, const Vector4 &maxp, const uint flag, const UintSet &oldflag) {
	for(uint i=0; i<m_qsize; i++) {
		if(!oldflag.includes(getQuadFlag(i))) continue;
		const Vector4 p = getQuadAverage(i);
		if(p.x < minp.x || maxp.x <= p.x) continue;
		if(p.y < minp.y || maxp.y <= p.y) continue;
		if(p.z < minp.z || maxp.z <= p.z) continue;
		if(p.t < minp.t || maxp.t <= p.t) continue;
		setQuadFlag(i, flag);
	}
}
void BuilderMesh::fillRectangleFlags(const Vector4 &minp, const Vector4 &maxp, const uint flag, const UintSet &oldflag) {
	fillNodeRectangleFlags(minp, maxp, flag, oldflag);
	fillEdgeRectangleFlags(minp, maxp, flag, oldflag);
	fillFaceRectangleFlags(minp, maxp, flag, oldflag);
	fillBodyRectangleFlags(minp, maxp, flag, oldflag);
	fillQuadRectangleFlags(minp, maxp, flag, oldflag);
}

void BuilderMesh::expandFlags(const uint flag, const UintSet &oldflag, uint layers) {
	uint i;
	uint qs = 0;
	Buffer<uint> q;
	for(i=0; i<m_qflag.size(); i++) {
		if(m_qflag[i] == flag) q.gather(i, qs);
	}
	uint bs = 0;
	Buffer<uint> b;
	for(i=0; i<m_bflag.size(); i++) {
		if(m_bflag[i] == flag) b.gather(i, bs);
	}
	uint fs = 0;
	Buffer<uint> f;
	for(i=0; i<m_fflag.size(); i++) {
		if(m_fflag[i] == flag) f.gather(i, fs);
	}
	uint es = 0;
	Buffer<uint> e;
	for(i=0; i<m_eflag.size(); i++) {
		if(m_eflag[i] == flag) e.gather(i, es);
	}
	uint ns = 0;
	Buffer<uint> n;
	for(i=0; i<m_nflag.size(); i++) {
		if(m_nflag[i] == flag) n.gather(i, ns);
	}

	while(ns + es + fs + bs + qs > 0) {
		while(qs > 0) {
			const Buffer<uint> &c = getQuadBodies(q[--qs]);
			for(i=0; i<c.size(); i++) {
				const uint iflag = getBodyFlag(c[i]);
				if(iflag == flag || !oldflag.includes(iflag)) continue;
				setBodyFlag(c[i], flag);
				b.gather(c[i], bs);
			}
		}
		while(bs > 0) {
			const Buffer<uint> &c = getBodyFaces(b[--bs]);
			for(i=0; i<c.size(); i++) {
				const uint iflag = getFaceFlag(c[i]);
				if(iflag == flag || !oldflag.includes(iflag)) continue;
				setFaceFlag(c[i], flag);
				f.gather(c[i], fs);
			}

			if(layers == 0) continue;
			const Buffer<uint> &cc = getBodyQuads(b[bs]);
			for(i=0; i<cc.size(); i++) {
				const uint iflag = getQuadFlag(cc[i]);
				if(iflag == flag || !oldflag.includes(iflag)) continue;
				setQuadFlag(cc[i], flag);
				q.gather(cc[i], qs);
			}
		}
		while(fs > 0) {
			const Buffer<uint> &c = getFaceEdges(f[--fs]);
			for(i=0; i<c.size(); i++) {
				const uint iflag = getEdgeFlag(c[i]);
				if(iflag == flag || !oldflag.includes(iflag)) continue;
				setEdgeFlag(c[i], flag);
				e.gather(c[i], es);
			}

			if(layers == 0) continue;
			const Buffer<uint> &cc = getFaceBodies(f[fs]);
			for(i=0; i<cc.size(); i++) {
				const uint iflag = getBodyFlag(cc[i]);
				if(iflag == flag || !oldflag.includes(iflag)) continue;
				setBodyFlag(cc[i], flag);
				b.gather(cc[i], bs);
			}
		}
		while(es > 0) {
			const Buffer<uint> &c = getEdgeNodes(e[--es]);
			for(i=0; i<c.size(); i++) {
				const uint iflag = getNodeFlag(c[i]);
				if(iflag == flag || !oldflag.includes(iflag)) continue;
				setNodeFlag(c[i], flag);
				n.gather(c[i], ns);
			}

			if(layers == 0) continue;
			const Buffer<uint> &cc = getEdgeFaces(e[es]);
			for(i=0; i<cc.size(); i++) {
				const uint iflag = getFaceFlag(cc[i]);
				if(iflag == flag || !oldflag.includes(iflag)) continue;
				setFaceFlag(cc[i], flag);
				f.gather(cc[i], fs);
			}
		}

		if(layers == 0) break;
		while(ns > 0) {
			const Buffer<uint> &cc = getNodeEdges(n[--ns]);
			for(i=0; i<cc.size(); i++) {
				const uint iflag = getEdgeFlag(cc[i]);
				if(iflag == flag || !oldflag.includes(iflag)) continue;
				setEdgeFlag(cc[i], flag);
				e.gather(cc[i], es);
			}
		}
		--layers;
	}
}

void BuilderMesh::removeByFlags(const UintSet &flag) {
	uint i;
	for(i=m_qsize; i-->0; ) {
		if(flag.includes(getQuadFlag(i))) removeQuad(i);
	}
	for(i=m_bsize; i-->0; ) {
		if(flag.includes(getBodyFlag(i))) removeBody(i);
	}
	for(i=m_fsize; i-->0; ) {
		if(flag.includes(getFaceFlag(i))) removeFace(i);
	}
	for(i=m_esize; i-->0; ) {
		if(flag.includes(getEdgeFlag(i))) removeEdge(i);
	}
	for(i=m_nsize; i-->0; ) {
		if(flag.includes(getNodeFlag(i))) removeNode(i);
	}
}

bool BuilderMesh::repeatMiddle(const Vector4 &pos, const Vector4 &dir, const Vector4 &step, const uint steps) {
	if(steps == 0) return true;
	BuilderMesh mesh(m_dim);
	if(mesh.createRepeated(*this, pos, dir, step, steps)) {
		swap(mesh);
		return true;
	}
	return false;
}

void BuilderMesh::combine(const Mesh &mesh) {
	uint i, j;

	const uint ns = getNodeSize();
	const uint nsize = mesh.getNodeSize();
	resizeNodeBuffer(ns + nsize);
	for(i=0; i<nsize; i++) {
		const uint n = addNode(mesh.getNodePosition(i));
		setNodeWeight(n, mesh.getNodeWeight(i));
		setNodeFlag(n, mesh.getNodeFlag(i));
	}
	const uint es = getEdgeSize();
	const uint esize = mesh.getEdgeSize();
	resizeEdgeBuffer(es + esize);
	for(i=0; i<esize; i++) {
		const Buffer<uint> &par = mesh.getEdgeNodes(i);
		const uint e = addEdge(ns + par[0], ns + par[1]);
		setEdgeFlag(e, mesh.getEdgeFlag(i));
	}
	const uint fs = getFaceSize();
	const uint fsize = mesh.getFaceSize();
	resizeFaceBuffer(fs + fsize);
	for(i=0; i<fsize; i++) {
		Buffer<uint> par = mesh.getFaceEdges(i);
		for(j=0; j<par.size(); j++) par[j] += es;
		const uint f = addFace(par);
		setFaceFlag(f, mesh.getFaceFlag(i));
	}
	const uint bs = getBodySize();
	const uint bsize = mesh.getBodySize();
	resizeBodyBuffer(bs + bsize);
	for(i=0; i<bsize; i++) {
		Buffer<uint> par = mesh.getBodyFaces(i);
		for(j=0; j<par.size(); j++) par[j] += fs;
		const uint b = addBody(par);
		setBodyFlag(b, mesh.getBodyFlag(i));
	}
	const uint qs = getQuadSize();
	const uint qsize = mesh.getQuadSize();
	resizeQuadBuffer(qs + qsize);
	for(i=0; i<qsize; i++) {
		Buffer<uint> par = mesh.getQuadBodies(i);
		for(j=0; j<par.size(); j++) par[j] += bs;
		const uint q = addQuad(par);
		setQuadFlag(q, mesh.getQuadFlag(i));
	}
}

void BuilderMesh::stretchLinear(const Vector4 &v, const uint steps, const UintSet &flag, const uint flagMiddle, const uint flagEnd) {
	if(steps == 0) return;
	const Vector4 d = v / double(steps);

	// find nodes to duplicate
	uint i, j;
	uint nnsize = 0;
	Buffer<uint> nn(m_nsize);
	for(i=0; i<m_nsize; i++) {
		if(flag.includes(getNodeFlag(i))) nn[nnsize++] = i;
	}
	if(nnsize == 0) return;
	nn.resize(nnsize);

	// create new nodes
	const uint nsize = m_nsize;
	m_nsize += steps * nnsize;
	resizeNodeBuffer(m_nsize);
	for(i=0; i<nnsize; i++) {
		const Vector4 p = getNodePosition(nn[i]);
		const double w = getNodeWeight(nn[i]);
		for(j=0; j<steps; j++) {
			const uint ii = nsize + nnsize * j + i;
			setNodePosition(ii, p + (j + 1) * d);
			setNodeWeight(ii, w);
		}
	}
	stretch(nn, steps, flag, flagMiddle, flagEnd);
}

void BuilderMesh::stretch(const Buffer<uint> &n, const uint steps, const UintSet &flag, const uint flagMiddle, const uint flagEnd) {
	uint i, j, k;
	const uint nnsize = n.size();
	const uint nsize = m_nsize - steps * nnsize;
	const uint esize = m_esize;
	const uint fsize = m_fsize;
	const uint bsize = m_bsize;
	const uint qsize = m_qsize;

	// find edges to stretch
	uint eesize = 0;
	Buffer<uint> ee(esize);
	for(i=0; i<esize; i++) {
		if(flag.includes(getEdgeFlag(i))) ee[eesize++] = i;
	}

	// find faces to stretch
	uint ffsize = 0;
	Buffer<uint> ff(fsize);
	for(i=0; i<fsize; i++) {
		if(flag.includes(getFaceFlag(i))) ff[ffsize++] = i;
	}

	// find bodies to stretch
	uint bbsize = 0;
	Buffer<uint> bb(bsize);
	for(i=0; i<bsize; i++) {
		if(flag.includes(getBodyFlag(i))) bb[bbsize++] = i;
	}

	// resize cell buffers
	m_esize = esize + steps * (nnsize + eesize);
	resizeEdgeBuffer(m_esize);
	m_fsize = fsize + steps * (eesize + ffsize);
	resizeFaceBuffer(m_fsize);
	m_bsize = bsize + steps * (ffsize + bbsize);
	resizeBodyBuffer(m_bsize);
	m_qsize = qsize + steps * bbsize;
	resizeQuadBuffer(m_qsize);

	// stretch nodes
	for(i=0; i<nnsize; i++) {
		const uint flag0 = getNodeFlag(n[i]);
		for(j=0; j<steps; j++) {
			const uint n1 = nsize + j * nnsize + i;
			const uint n0 = (j > 0 ? n1 - nnsize : n[i]);
			setNodeFlag(n1, flag0 + (j + 1 < steps ? flagMiddle : flagEnd));

			const uint ei = esize + j * nnsize + i;
			setEdgeFlag(ei, flag0 + flagMiddle);
			m_e[ei].n.resize(2);
			m_e[ei].n[0] = n0;
			m_e[ei].n[1] = n1;
			m_n[n0].e.push_back(ei);
			m_n[n1].e.push_back(ei);
		}
	}

	// stretch edges
	for(i=0; i<eesize; i++) {
		const uint flag0 = getEdgeFlag(ee[i]);
		for(j=0; j<steps; j++) {
			const uint e1 = esize + steps * nnsize + j * eesize + i;
			const uint e0 = (j > 0 ? e1 - eesize : ee[i]);
			Buffer<uint> en = getEdgeNodes(e0);
			Buffer<uint> fe(4);
			fe[1] = e0;
			fe[3] = e1;
			for(k=0; k<en.size(); k++) {
				if(j == 0) fe[2 * k] = getNodeEdges(en[k]).back();
				else fe[2 * k] = getNodeEdges(en[k])[1];
				en[k] = getEdgeNodes(fe[2 * k])[1];
			}
			setEdgeFlag(e1, flag0 + (j + 1 < steps ? flagMiddle : flagEnd));
			m_e[e1].n = en;
			m_n[en[0]].e.push_back(e1);
			m_n[en[1]].e.push_back(e1);

			const uint fi = fsize + j * eesize + i;
			setFaceFlag(fi, flag0 + flagMiddle);
			m_f[fi].e = fe;
			for(k=0; k<fe.size(); k++) m_e[fe[k]].f.push_back(fi);
		}
	}

	// stretch faces
	for(i=0; i<ffsize; i++) {
		const uint flag0 = getFaceFlag(ff[i]);
		for(j=0; j<steps; j++) {
			const uint f1 = fsize + steps * eesize + j * ffsize + i;
			const uint f0 = (j > 0 ? f1 - ffsize : ff[i]);
			Buffer<uint> fe = getFaceEdges(f0);
			Buffer<uint> bf(2 + fe.size());
			bf[0] = f0;
			bf.back() = f1;
			for(k=0; k<fe.size(); k++) {
				if(j == 0) bf[1 + k] = getEdgeFaces(fe[k]).back();
				else bf[1 + k] = getEdgeFaces(fe[k])[1];
				fe[k] = getFaceEdges(bf[1 + k])[3];
			}
			setFaceFlag(f1, flag0 + (j + 1 < steps ? flagMiddle : flagEnd));
			m_f[f1].e = fe;
			for(k=0; k<fe.size(); k++) m_e[fe[k]].f.push_back(f1);

			const uint bi = bsize + j * ffsize + i;
			setBodyFlag(bi, flag0 + flagMiddle);
			m_b[bi].f = bf;
			for(k=0; k<bf.size(); k++) m_f[bf[k]].b.push_back(bi);
		}
	}

	// stretch bodies
	for(i=0; i<bbsize; i++) {
		const uint flag0 = getBodyFlag(bb[i]);
		for(j=0; j<steps; j++) {
			const uint b1 = bsize + steps * ffsize + j * bbsize + i;
			const uint b0 = (j > 0 ? b1 - bbsize : bb[i]);
			Buffer<uint> bf = getBodyFaces(b0);
			Buffer<uint> qb(2 + bf.size());
			qb[0] = b0;
			qb.back() = b1;
			for(k=0; k<bf.size(); k++) {
				if(j == 0) qb[1 + k] = getFaceBodies(bf[k]).back();
				else qb[1 + k] = getFaceBodies(bf[k])[1];
				bf[k] = getBodyFaces(qb[1 + k]).back();
			}
			setBodyFlag(b1, flag0 + (j + 1 < steps ? flagMiddle : flagEnd));
			m_b[b1].f = bf;
			for(k=0; k<bf.size(); k++) m_f[bf[k]].b.push_back(b1);

			const uint qi = qsize + j * bbsize + i;
			setQuadFlag(qi, flag0 + flagMiddle);
			m_q[qi].b = qb;
			for(k=0; k<qb.size(); k++) m_b[qb[k]].q.push_back(qi);
		}
	}
}

void BuilderMesh::improveNodeByHodge(const uint n, const bool position, const bool weight) {
	// changes node position and weight to improve Hodge quality
	uint i, j;
	Buffer<double> sq;
	Buffer<Vector4> v;
	Buffer<Vector4> d;
	double fac;

	if(getEdgeSize() == 0) return;
	if(getFaceSize() == 0) { // 1-dimensional mesh
		const Buffer<uint> &e = getNodeEdges(n);
		sq.resize(e.size());
		v.resize(e.size());
		d.resize(e.size());
		for(i=0; i<e.size(); i++) {
			const uint opp = getEdgeOtherNode(e[i], n);
			const Vector4 p = getNodePosition(opp);
			sq[i] = getNodeWeight(opp);
			v[i] = getNodePosition(n) - p;
			d[i] = v[i];
		}
		fac = 1.0;
	}
	else if(getBodySize() == 0) { // 2-dimensional mesh
		const Buffer<uint> f = getNodeFaces(n);
		const Buffer<uint> &e = getNodeEdges(n);
		sq.resize(f.size());
		v.resize(f.size());
		d.resize(f.size());
		for(i=0; i<f.size(); i++) {
			const Buffer<uint> &opp = getFaceEdges(f[i]);
			if(opp.size() != 3) return;
			for(j=0; e.includes(opp[j]); j++); // find the opposite cell
			const Vector4 p = getEdgePosition(opp[j]);
			sq[i] = getRadiusSq(p, getEdgeAnyNode(opp[j]));
			v[i] = getNodePosition(n) - p;
			d[i] = getEdgeDeviation(opp[j], v[i]);
		}
		fac = 1.5;
	}
	else if(getQuadSize() == 0) { // 3-dimensional mesh
		const Buffer<uint> b = getNodeBodies(n);
		const Buffer<uint> f = getNodeFaces(n);
		sq.resize(b.size());
		v.resize(b.size());
		d.resize(b.size());
		for(i=0; i<b.size(); i++) {
			const Buffer<uint> &opp = getBodyFaces(b[i]);
			if(opp.size() != 4) return;
			for(j=0; f.includes(opp[j]); j++); // find the opposite cell
			const Vector4 p = getFacePosition(opp[j]);
			sq[i] = getRadiusSq(p, getFaceAnyNode(opp[j]));
			v[i] = getNodePosition(n) - p;
			d[i] = getFaceDeviation(opp[j], v[i]);
		}
		fac = 2.0;
	}
	else {
		const Buffer<uint> q = getNodeQuads(n);
		const Buffer<uint> b = getNodeBodies(n);
		sq.resize(q.size());
		v.resize(q.size());
		d.resize(q.size());
		for(i=0; i<q.size(); i++) {
			const Buffer<uint> &opp = getQuadBodies(q[i]);
			if(opp.size() != 5) return;
			for(j=0; b.includes(opp[j]); j++); // find the opposite cell
			const Vector4 p = getBodyPosition(opp[j]);
			sq[i] = getRadiusSq(p, getBodyAnyNode(opp[j]));
			v[i] = getNodePosition(n) - p;
			d[i] = getBodyDeviation(opp[j], v[i]);
		}
		fac = 2.5;
	}

	if(position) {
		double sumcost = 0.0; // relative value of cost function
		Vector4 sumgrad(0,0,0,0); // one sixth of the gradient
		SymMatrix4 sumdiv(ZEROSYMMATRIX4); // a matrix to compute the zero of gradient
		for(i=0; i<v.size(); i++) {
			const Vector4 tv = getTransformed(v[i]);
			const Vector4 td = getTransformed(d[i]);
			const double dsq = 1.0 / d[i].dot(td);
			const double rsq = 0.5 * dsq * (v[i].dot(tv) + getNodeWeight(n) - sq[i]);
			sumcost += rsq * rsq / dsq;
			const Vector4 grad = (tv - rsq * td);
			sumgrad += rsq * grad;
			sumdiv += dsq * grad.outerProduct() + rsq * getMetric();
		}

		const Vector4 dp = sumgrad * sumdiv.inverse();
		setNodePosition(n, getNodePosition(n) - dp);
		for(i=0; i<v.size(); i++) {
			const Vector4 td = getTransformed(d[i]);
			v[i] -= dp;
			d[i] *= v[i].dot(td) / d[i].dot(td);
		}
	}
	if(weight) {
		double sumdiv = 0.0;
		double sumw = 0.0;
		for(i=0; i<v.size(); i++)
		{
			const double dsq = fac / d[i].dot(getTransformed(d[i]));
			sumw += dsq * (sq[i] - v[i].dot(getTransformed(v[i]))) + 1.0;
			sumdiv += dsq;
		}
		setNodeWeight(n, sumw / sumdiv);
	}
}

void BuilderMesh::optimizeNodes(const UintSet &flag, const uint iterations, const bool position, const bool weight) {
	// select nodes
	uint ns = 0;
	Buffer<uint> n(getNodeSize());
	for(uint i=0; i<n.size(); i++) {
		if(flag.includes(getNodeFlag(i))) n[ns++] = i;
	}
	n.resize(ns);

	// optimize
	bool improves = true;
	for(uint i=0; i<iterations && improves; i++) {
		improves = false;
		if(weight && optimizeNodesIteration(n, false, true)) improves = true;
		if(position && optimizeNodesIteration(n, true, false)) improves = true;
	}
}

bool BuilderMesh::optimizeNodesIteration(const Buffer<uint> &n, const bool position, const bool weight) {
	uint i, j;
	Buffer<bool> chk; // necessary cells
	if(m_fsize == 0) chk.resize(getEdgeSize());
	else if(m_bsize == 0) chk.resize(getFaceSize());
	else if(m_qsize == 0) chk.resize(getBodySize());
	else chk.resize(getQuadSize());
	chk.fill(false);

	// compute gradient and sum of radius squares
	const double step = 1e-5;
	double fcurr = 0.0;
	double grsq = 0.0; // gradient square
	Buffer<double> grw((weight ? n.size() : 0), 0.0); // weight gradient
	Buffer<Vector4> grp((position ? n.size() : 0), Vector4(0,0,0,0)); // position gradient
	for(i=0; i<n.size(); i++) {
		Buffer<uint> ele;
		if(m_fsize == 0) ele = getNodeEdges(n[i]);
		else if(m_bsize == 0) ele = getNodeFaces(n[i]);
		else if(m_qsize == 0) ele = getNodeBodies(n[i]);
		else ele = getNodeQuads(n[i]);

		double r0 = 0.0;
		for(j=0; j<ele.size(); j++) {
			const double r = getCellLengthSq(ele[j]);
			r0 += r;
			if(chk[ele[j]]) continue;
			chk[ele[j]] = true;
			fcurr += r;
		}
		if(weight) {
			const double w = getNodeWeight(n[i]);
			setNodeWeight(n[i], w + step);
			double rw = 0.0;
			for(j=0; j<ele.size(); j++) rw += getCellLengthSq(ele[j]);
			grw[i] = (rw - r0) / step;
			grsq += grw[i] * grw[i];				
			setNodeWeight(n[i], w);
		}
		if(position) {
			const Vector4 p = getNodePosition(n[i]);
			setNodePosition(n[i], p + Vector4(step,0,0,0));
			double rx = 0.0;
			for(j=0; j<ele.size(); j++) rx += getCellLengthSq(ele[j]);
			grp[i].x = (rx - r0) / step;
			if(m_dim >= 2) {
				setNodePosition(n[i], p + Vector4(0,step,0,0));
				double ry = 0.0;
				for(j=0; j<ele.size(); j++) ry += getCellLengthSq(ele[j]);
				grp[i].y = (ry - r0) / step;
			}
			if(m_dim >= 3) {
				setNodePosition(n[i], p + Vector4(0,0,step,0));
				double rz = 0.0;
				for(j=0; j<ele.size(); j++) rz += getCellLengthSq(ele[j]);
				grp[i].z = (rz - r0) / step;
			}
			if(m_dim == 4) {
				setNodePosition(n[i], p + Vector4(0,0,0,step));
				double rt = 0.0;
				for(j=0; j<ele.size(); j++) rt += getCellLengthSq(ele[j]);
				grp[i].t = (rt - r0) / step;
			}
			grsq += grp[i].lensq();
			setNodePosition(n[i], p);
		}
	}
	cout << "Minimize " << fcurr << endl;
	if(grsq < 1e-20) return false;

	double xcurr = 0.0;
	double xprev = -step / std::sqrt(grsq);
	double fprev = fcurr - xprev * grsq;
	for(uint iter=0; iter<1000; iter++) {
		const double xnext = 3.0 * xcurr - 2.0 * xprev;
		for(i=0; i<n.size(); i++) {
			if(weight) setNodeWeight(n[i], getNodeWeight(n[i]) - (xnext - xcurr) * grw[i]);
			if(position) setNodePosition(n[i], getNodePosition(n[i]) - (xnext - xcurr) * grp[i]);
		}
		double fnext = 0.0;
		for(i=0; i<chk.size(); i++) {
			if(chk[i]) fnext += getCellLengthSq(i);
		}
		if(fnext > fcurr) {
			const double xdiff = (xnext - xcurr) * (fprev - 2.25 * fcurr + 1.25 * fnext) / (2.0 * fprev - 3.0 * fcurr + fnext);
			for(i=0; i<n.size(); i++) {
				if(weight) setNodeWeight(n[i], getNodeWeight(n[i]) + xdiff * grw[i]);
				if(position) setNodePosition(n[i], getNodePosition(n[i]) + xdiff * grp[i]);
			}
			return (iter > 0);
		}
		xprev = xcurr;
		fprev = fcurr;
		xcurr = xnext;
		fcurr = fnext;
	}
	cout << "BuilderMesh::optimizeNodesIteration -> Iteration limit exceeded." << endl;
	return true;
}

double BuilderMesh::getCellLengthSq(const uint i) const {
	if(m_fsize == 0) return getLengthSq(getEdgePosition(i), getEdgeAnyNode(i));
	if(m_bsize == 0) return getLengthSq(getFacePosition(i), getFaceAnyNode(i));
	if(m_qsize == 0) return getLengthSq(getBodyPosition(i), getBodyAnyNode(i));
	return getLengthSq(getQuadPosition(i), getQuadAnyNode(i));
}
