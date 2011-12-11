#ifndef MAP_H
#define MAP_H

#include "State.h"

template<typename T>
struct Map {
	typedef T Type;
	typedef std::vector<Type> Vector;

	Map(uint r = 0, uint c = 0)
	{
		resize(r, c);
	}

	void resize(uint r, uint c)
	{
		rows = r;
		cols = c;
		m_data.resize(rows * cols);
	}

	void fill(const Type & t)
	{
		m_data.assign(rows * cols, t);
	}

	typename Vector::reference operator[](const Location & pos)
	{
		return m_data[pos.row * cols + pos.col];
	}

	typename Vector::const_reference operator[](const Location & pos) const
	{
		return m_data[pos.row * cols + pos.col];
	}

	uint rows;
	uint cols;
	std::vector<Type> m_data;
};


#endif // MAP_H
