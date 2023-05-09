/*
 * point.h
 *
 *  Created on: 10 Jan 2023
 *      Author: mac
 */

#ifndef UTILITY_POINT_H_
#define UTILITY_POINT_H_

#include <iostream>

namespace manet {

class point {
public:
	double x, y;
	point(double d, double e) : x(d), y(e) {}
	friend std::ostream& operator<<(std::ostream& os, const point& pt);
};

//double sqrDist(point a, point b);
//double sqrDist(double x0, double y0, double x1, double y1);
//

double sqrDist(point a, point b) {
	return (a.x - b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y);
}
double sqrDist(double x0, double y0, double x1, double y1) {
	return (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1);
}

std::ostream& operator<<(std::ostream& os, const point& pt) {
	os << '(' << pt.x << ',' << pt.y << ')';
	return os;
}
};

#endif /* UTILITY_POINT_H_ */
