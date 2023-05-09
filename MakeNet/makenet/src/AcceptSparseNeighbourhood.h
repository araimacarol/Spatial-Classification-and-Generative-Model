/*
 * AcceptSparseNeighbourhood.h
 *
 *  Created on: 10 Jan 2023
 *      Author: mac
 */

#ifndef ACCEPTSPARSENEIGHBOURHOOD_H_
#define ACCEPTSPARSENEIGHBOURHOOD_H_

#include <stdio.h>
#include <string>
#include <vector>

#include "QuadTree.h"
//#include "AcceptPlacement.h"

namespace manet {

template<typename T>
class AcceptSparseNeighbourhood : public AcceptPlacement<T> {
protected:
	double radius;
	int maxNeighbours;
	std::vector<std::pair<point, int>> D;
public:
	AcceptSparseNeighbourhood(double r, int n) : radius(r), maxNeighbours(n) {}
	virtual ~AcceptSparseNeighbourhood() {}

	bool operator()(point pt) {
		D.clear();
		AcceptPlacement<T>::QT->putPointsInRadius(D, pt, radius);
		std::cout << "D.size() = " << D.size() << std::endl;
		return (D.size() <= maxNeighbours);
	}	// override this in derived classes for different acceptance functions

	std::string name() { return "AcceptSparseNeighbourhood"; }

	bool accept(point pt) {
		D.clear();
		AcceptPlacement<T>::QT->putPointsInRadius(D, pt, radius);
		std::cout << "D.size() = " << D.size() << std::endl;
		if (D.size() > maxNeighbours) {
			std::cout << "\tToo many neighbours!" << std::endl;
		} else {
			std::cout << "\tThis place is nice and quiet." << std::endl;
		}
		return (D.size() <= maxNeighbours);
	}

};

}

#endif /* ACCEPTSPARSENEIGHBOURHOOD_H_ */
