/*
 * QuadTreeNode.h
 *
 *  Created on: 6 Jan 2023
 *      Author: mac
 */

#ifndef QUADTREENODE_H_
#define QUADTREENODE_H_

#include <iostream>
#include <stddef.h>
#include <vector>

#include "../utility/appexception.h"
#include "../utility/point.h"

//#include "QuadTree.h"
//template<typename T>
//class QuadTree;

//#include "AcceptPlacement.h"

namespace manet {

//class BoundingBox {
//protected:
////	double x0, x1, y0, y1;
//public:
//	point a, b;
//	BoundingBox() : a(0.0, 0.0), b(1.0,1.0) {}
//	BoundingBox(double nux0, double nux1, double nuy0, double nuy1) : a(nux0, nux1), b(nuy0, nuy1) {}
//
//	template<typename T>
////	bool contains(std::forward_iterator& begin, std::forward_iterator& end) { return false; }
//
//	bool contains(double x, double y) { return contains(point(x, y)); }
//	bool contains(point pt) {
////		std::cout << "\tboundary = " << *this << "; point is " << pt << std::endl;
//		bool containsX = (a.x <= pt.x) && (pt.x <= b.x);
////		std::cout << "\t\ta.x <= pt.x ? " << (a.x <= pt.x) << std::endl;
////		std::cout << "\t\tpt.x <= b.x ? " << (pt.x <= b.x) << std::endl;
////		std::cout << "\tcontainsX = " << containsX << std::endl;
//		bool containsY = (a.y <= pt.y) && (pt.y <= b.y);
////		std::cout << "\t\ta.y <= pt.y ? " << (a.y <= pt.y) << std::endl;
////		std::cout << "\t\tpt.y <= b.y ? " << (pt.y <= b.y) << std::endl;
////		std::cout << "\tcontainsY = " << containsY << std::endl;
//		return containsX && containsY;
//	}
//	bool contains(BoundingBox box) {
//		return (contains(box.a) && contains(box.b));
//	}
//	friend std::ostream& operator<<(std::ostream& os, const BoundingBox& box);
//};
//
//template<typename T>
//class AcceptPlacement {
//	QuadTree<T>* QT;
//	explicit AcceptPlacement() : QT(nullptr) {}
//public:
//	virtual ~AcceptPlacement() {}
//
//	virtual void setQuadTree(QuadTree<T>* qt) { QT = qt; }
//
//	bool accept(point pt);
//	bool operator()(point pt) const { return false; }	// override this in derived classes for different acceptance functions
//};
//
//template<typename T>
//class QuadTree {
//protected:
//	QuadTreeNode<T> root;
//public:
//	QuadTree<T>(int cap=10) : root(cap) {}
//	virtual ~QuadTree<T>() {}
//
//	QuadTreeNode<T>& getRoot() { return root; }
//	const QuadTreeNode<T>& getRoot() const { return root; }
//
//	bool insert(point p, T d) {
//		return root.insert(p, d);
//	}
//
//	friend std::ostream& operator<<(std::ostream& os, QuadTree<T>& Q) {
//		os << Q.root;
//		return os;
//	}
//
//	void printAll() { root.printAll(); }
//
//	void putPointsInRadius(std::vector<std::pair<point, T> >& ptsInRange, point p, double r) {
//		root.putPointsInRadius(ptsInRange, p, r);
//	}
//	void putPointsInRange(std::vector<std::pair<point, T> >& ptsInRange, BoundingBox range) {
//		root.putPointsInRange(ptsInRange, range);
//	}
//
//	void traverse(void* func()) { root.traverse(func); }
//};
//
//template<typename T>
//class QuadTreeNode {
//protected:
////	QuadTree<T>* QT;
//	QuadTreeNode *NW, *NE, *SW, *SE;
//	std::vector<std::pair<point, T> > data;
//	BoundingBox boundary;
//	size_t capacity;
//public:
//	explicit QuadTreeNode(int cap=10) : NW(nullptr), NE(nullptr), SW(nullptr), SE(nullptr), capacity(cap) {}
//	QuadTreeNode(double nux0, double nux1, double nuy0, double nuy1, size_t cap=10) :
//		boundary(nux0, nux1, nuy0, nuy1), capacity(cap) {
////		std::cout << "New QNode at " << boundary << std::endl;
//		NE = NW = SE = SW = nullptr;
//	}
//	virtual ~QuadTreeNode() {}
//
//	bool insert(std::pair<point, T> pr) {
////		std::cout << "Inserting point " << pr.first << " by calling insert(point, T data)" << std::endl;
//		return insert(pr.first, pr.second, nullptr);
//	}
//	bool insert(point pt, T d, AcceptPlacement<T>* acc = nullptr) {
//		if (!boundary.contains(pt)) {
//			return false;
//		}
//		if (capacity > data.size()) {
//			// here, test if the new placement will be accepted:
//			if (acc == nullptr) {
//				data.push_back(std::pair<point, T>(pt, d));
//				return true;
//			} else if (acc->accept(pt)) {
//				data.push_back(std::pair<point, T>(pt, d));
//			} else {
//				return false;
//			}
////			std::cout << "Inserting data " << d << " at point " << pt << std::endl;
//		}
//
//		if (isLeaf()) {
//			split();
//		}
//		if (NW->insert(pt, d, acc)) {
//			return true;
//		}
//		if (NE->insert(pt, d, acc)) {
//			return true;
//		}
//		if (SW->insert(pt, d, acc)) {
//			return true;
//		}
//		if (SE->insert(pt, d, acc)) {
//			return true;
//		}
//		return false;
//	}
//
//	bool intersects(BoundingBox& box) {
////		std::cout << "This boundary = " << boundary << "; box = " << box << std::endl;
//		bool intersectX = (boundary.a.x <= box.b.x) && (box.a.x <= boundary.b.x);
////		std::cout << "\tintersectX = " << intersectX << std::endl;
//		bool intersectY = (boundary.a.y <= box.b.y) && (box.a.y <= boundary.b.y);
////		std::cout << "\tintersectY = " << intersectY << std::endl;
//		return (intersectX && intersectY);
//	}
//
//	int intersects(point p, double r) {
//		// return true iff this box overlaps a circle centered at p with radius r.
//		// First check if the boundary's lie completely outside the
//		if (p.x + r < boundary.a.x) return 0;
//		if (p.x - r > boundary.b.x) return 0;
//		if (p.y + r < boundary.a.y) return 0;
//		if (p.y - r > boundary.b.y) return 0;
//
//		// next check if corners are within range r of the centre:
//		int cornersInRange(0);
//		double r2 = r*r;
//		cornersInRange += (sqrDist(p, boundary.a) <= r2) ? 1 : 0;
//		cornersInRange += (sqrDist(p, boundary.b) <= r2) ? 1 : 0;
//		cornersInRange += (sqrDist(p.x, p.y, boundary.a.x, boundary.b.y) <= r2) ? 1 : 0;
//		cornersInRange += (sqrDist(p.x, p.y, boundary.b.x, boundary.a.y) <= r2) ? 1 : 0;
//		if (cornersInRange > 0) {
//			return cornersInRange;	// XXX Still need to describe full or partial overlap
//		}
//		/*
//		 *  remaining possibility for a non-empty intersection is that 2 or 4 corners of the square lie
//		 *  *outside* the circle but *inside* the minimal square containing it:
//		 */
//		if ((boundary.a.x <= p.x) && (p.x <= boundary.b.x)) return 1;
//		if ((boundary.a.y <= p.y) && (p.y <= boundary.b.y)) return 1;
//		return 0;	// XXX STUB
//	}
//
//	inline bool isLeaf() const { return (NW == nullptr); }
//
//	friend std::ostream& operator<<(std::ostream& os, QuadTreeNode<T>& QN) {
//		if (QN.isLeaf()) {
//			for (auto pr : QN.data) {
//				os << pr.first << ',' << pr.second << ',' << QN.boundary << std::endl;
//			}
////			os << QN.boundary << "\t{ ";
////			for (auto pr : QN.data) {
////				os << pr.second << "@" << pr.first << ' ';
////			}
////			os << '}' << std::endl;
//			return os;
//		}
//		os << *(QN.NW);
//		os << *(QN.NE);
//		os << *(QN.SW);
//		os << *(QN.SE);
//		return os;
//	}
//
//	void printAll() {
//		if (isLeaf()) {
//			printData();
//			return;
//		}
//		NW->printAll();
//		NE->printAll();
//		SW->printAll();
//		SE->printAll();
//	}
//
//	void printData() {
////		std::cout << boundary << std::endl; //"[(" << boundary.a.x << ',' << boundary.a.y << "),(" << boundary.b.x << ',' << boundary.b.y << ")]" << std::endl;
//		for (auto pr : data) {
//			std::cout << pr.second << " @ " << pr.first << " in " << boundary << std::endl;
//		}
//	}
//
//	void putAllPoints(std::vector<std::pair<point, T>>& pts) {
//		if (isLeaf()) {
//			pts.insert(pts.begin(), data.begin(), data.end());
////			std::cout << "adding all " << data.size() << " points at leaf" << std::endl;
//			return;
//		}
//		NW->putAllPoints(pts);
//		NE->putAllPoints(pts);
//		SW->putAllPoints(pts);
//		SE->putAllPoints(pts);
//	}
//
//	void putPointsInRadius(std::vector<std::pair<point, T> > & ptsInRadius, point p, double r) {
//		switch (intersects(p, r)) {
////			BoundingBox box(p.x-r, p.y-r, p.x+r, p.y+r);
//			case 0: {
//				return;
//			case 1:
//			case 2:
//			case 3:
//				if (isLeaf()) {
//					// test EACH point to see if it's within r of p
//					for (auto pr : data) {
//						if (sqrDist(p, pr.first) <= r*r) {
//							ptsInRadius.push_back(pr);
////							std::cout << "*";
//						}
//					}
////					std::cout << std::endl;
//				} else {
//					NW->putPointsInRadius(ptsInRadius, p, r);
//					NE->putPointsInRadius(ptsInRadius, p, r);
//					SW->putPointsInRadius(ptsInRadius, p, r);
//					SE->putPointsInRadius(ptsInRadius, p, r);
//				}
//				break;
//			}
//			case 4:
//				putAllPoints(ptsInRadius);
//				break;
//			default:
//				break;
//		}
//	}
//
//	void putPointsInRange(std::vector<std::pair <point, T>>& ptsInRange, BoundingBox& box) {
//		if (!intersects(box)) {
//			return;
//		}
//		if (isLeaf()) {
////			printData();
////			std::cout << "Adding " << data.size() << " points to the collection." << std::endl;
//			if (box.contains(boundary) || boundary.contains(box)) {
//				ptsInRange.insert(ptsInRange.begin(), data.begin(), data.end());	// XXX NOT correct just yet
//			} else {
//				for (auto pr : data) {
//					if (box.contains(pr.first)) {
//						ptsInRange.push_back(pr);
//					}
//				}
//			}
//			return;
//		}
//		NW->putPointsInRange(ptsInRange, box);
//		NE->putPointsInRange(ptsInRange, box);
//		SW->putPointsInRange(ptsInRange, box);
//		SE->putPointsInRange(ptsInRange, box);
//	}
//
////	void setQuadTree(	QuadTree<T>* qt) { QT = qt; }
//
//	void split() {
////		std::cout << "Splitting up the node with boundary " << boundary << std::endl;
//		double xMid(0.5*(boundary.a.x + boundary.b.x));
//		double yMid(0.5*(boundary.a.y + boundary.b.y));
//		NW = new QuadTreeNode(boundary.a.x, yMid, xMid, boundary.b.y, capacity);
//		NE = new QuadTreeNode(xMid, yMid, boundary.b.x, boundary.b.y, capacity);
//		SW = new QuadTreeNode(boundary.a.x, boundary.a.y, xMid, yMid, capacity);
//		SE = new QuadTreeNode(xMid, boundary.a.y, boundary.b.x, yMid, capacity);
//		// put all the data in this node into the leaves:
//		for (auto pr : data) {
////			std::cout << "Copying data " << pr.first << ',' << pr.second << std::endl;
//			if (NW->insert(pr)) {
////				std::cout << "Successfully put in NW quadrant" << std::endl;
//				continue;
//			}
//			if (NE->insert(pr)) {
////				std::cout << "Successfully put in NE quadrant" << std::endl;
//				continue;
//			}
//			if (SW->insert(pr)) {
////				std::cout << "Successfully put in SW quadrant" << std::endl;
//				continue;
//			}
//			if (SE->insert(pr)) {
////				std::cout << "Successfully put in SE quadrant" << std::endl;
//				continue;
//			}
//			throw new app_exception("Failed to move data into a leaf. WEIRD.");
//		}
//		data.clear();
//	}
//
//
//	void traverse(void* func()) {
//		if (isLeaf()) {
//			return;
//		}
//		NW->traverse(func);
//		NE->traverse(func);
//		SW->traverse(func);
//		SE->traverse(func);
//	}
//
//};
//
//template<typename T>
//bool AcceptPlacement<T>::accept(point pt) {
//	return true;
//}

} /* namespace manet */

#endif /* QUADTREENODE_H_ */
