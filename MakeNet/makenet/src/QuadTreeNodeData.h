/*
 * QuadTreeNodeData.h
 *
 *  Created on: 6 Jan 2023
 *      Author: mac
 */

#ifndef QUADTREENODEDATA_H_
#define QUADTREENODEDATA_H_

#include <stddef.h>
#include <vector>

namespace manet {

template <typename T>
class QuadTreeNodeData {
protected:
	double x, y;
	T* data;
	explicit QuadTreeNodeData() : x(0.0), y(0.0), data(nullptr) {}
public:
	QuadTreeNodeData(double nux, double nuy, T* d) : x(nux), y(nuy), data(d) {}
	virtual ~QuadTreeNodeData() {}
};


class BoundingBox {
protected:
	double x0, x1, y0, y1;
public:
	BoundingBox() : x0(0.0), x1(1.0), y0(0.0), y1(1.0) {}
	BoundingBox(double nux0, double nux1, double nuy0, double nuy1) : x0(nux0), x1(nux1), y0(nuy0), y1(nuy1) {}

	template<typename T>
	bool contains(const QuadTreeNodeData<T>& qd) { return false; }

	bool contains(double x, double y) { return false; }
};

} /* namespace manet */

#endif /* QUADTREENODEDATA_H_ */
