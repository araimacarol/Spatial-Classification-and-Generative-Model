#include <iostream>
#include <stdio.h>
#include <vector>

#include "../utility/appexception.h"
#include "../utility/myrandom.h"

#include "QuadTree.h"
#include "AcceptSparseNeighbourhood.h"

using namespace std;
using namespace manet;

bool _debugging(false);

void testInsertingData() {
	cout << "Testing basic insert" << endl;
	int cap(1024);
	int ss(1000);
	QuadTree<int> QT(cap);
	cout << "#samples = " << ss << endl;
	cout << "capacity / cell " << cap << endl;
	int numInserted(0);
	for (int i = 0; i < ss; ++i) {
		point p(dran(), dran());
		if (QT.insert(p, i)) {
			++numInserted;
		} else {
			cout << "Item " << i << " at position " << p << " could not be inserted." << endl;
			throw new app_exception("Insertion failure!");
		}
	}
	QT.printAll();
	cout << "#numInserted = " << numInserted << endl;
}

void testDataInRange() {
	cout << "Hello World" << endl;
	int cap(1);
	int ss(1000);
	QuadTree<int> QT(cap);
	int numInserted(0);
	for (int i = 0; i < ss; ++i) {
		point p(dran(), dran());
		if (QT.insert(p, i)) {
			++numInserted;
		} else {
			cout << "Item " << i << " at position " << p << " could not be inserted." << endl;
			throw new app_exception("Insertion failure!");
		}
	}
	cout << "#numInserted = " << numInserted << endl;
	vector<std::pair<point, int>> D;
	BoundingBox box(0.4,0.4,0.6,0.6);
	QT.putPointsInRange(D, box);
	cout << "Selected points in " << box << ':' << endl;
	for (auto pr : D) {
		std::cout << pr.second << " @ " << pr.first << std::endl;
	}
	QT.printAll();
	cout << "Number added to D = " << D.size() << endl;
}

void testIntersectionCases() {
	point p(0.5,0.5);
	BoundingBox B(0.1,0.1,0.2,0.2);
	vector< QuadTreeNode<int> > Box;
	Box.push_back(QuadTreeNode<int>(0.75,0.15,0.95,0.35));	// should NOT intersect (case 1)
	Box.push_back(QuadTreeNode<int>(0.7,0.51,0.9,0.71));		// SHOULD intersect (case 2)
	Box.push_back(QuadTreeNode<int>(0.45,0.4,0.65,0.6));		// SHOULD intersect (case 3)
	Box.push_back(QuadTreeNode<int>(0.0,0.4,0.3,0.7));			// SHOULD intersect (case 4(a))
	Box.push_back(QuadTreeNode<int>(0.31,0.31,0.69,0.69));	// SHOULD intersect (case 4(b))
	Box.push_back(QuadTreeNode<int>(0.2,0.05,0.9,0.75));		// SHOULD intersect (case 5)
	for (auto b : Box) {
		cout << "intersects " << b.intersects(p,0.22) << endl;
	}
}

void testGettingPointsInRadius() {
	time_t start = clock();

	point p(0.5,0.5);
	int cap(256);
	int ss(1000000);
	QuadTree<int> QT(cap);
	for (int i = 0; i < ss; ++i) {
		point p(dran(), dran());
		if (!QT.insert(p, i)) {
			cout << "Value " << i << " at point " << p << " could not be added!" << endl;
		}
	}
	cout << "Inserting " << ss << " items into QuadTree: " << static_cast<float>(clock() - start) / CLOCKS_PER_SEC
			<< " seconds" << endl;
	start = clock();

	vector<std::pair<point, int>> D;
	for (double r(0.1); r < 0.5; r+= 0.1) {
		for (int i = 0; i < 5; ++i) {
			D.clear();
			p.x = dran(1.0-2*r) + r;
			p.y = dran(1.0-2*r) + r;
			QT.putPointsInRadius(D,p,r);
			cout << "Finding all " << D.size() << " points distance " << r << " from p=" << p << ": "  << static_cast<float>(clock() - start) / CLOCKS_PER_SEC
					<< " seconds" << endl;
			start = clock();
		}
	}
//	for (auto pr : D) {
//		cout << pr.first << ',' << pr.second << std::endl;
//	}
//	cout << std::endl;
//	cout << QT;
}

void testAddingSpacedNeighbours() {
	QuadTree<int> QT(16);
	int n(1000);
	int i(0);
	AcceptPlacement<int>* acc = new AcceptSparseNeighbourhood<int>(0.1, 20);
	acc->setQuadTree(&QT);
	int failures(0);
	while (i < n) {
		point pt(dran(), dran());
		if (QT.insert(pt, i, acc)) {
			++i;
			failures = 0;
		} else {
			++failures;
		}
		if (failures > 100) {
			break;
		}
	}
	QT.printAll();
	cout << "number of vertices added = " << i << endl;
}

int main(int argn, char** argc) {
//	testInsertingData();
//	testDataInRange();
//	testGettingPointsInRadius();
	testAddingSpacedNeighbours();
	return 0;
}
