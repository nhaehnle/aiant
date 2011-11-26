#ifndef MAP_H
#define MAP_H

#include "State.h"

template<typename T>
struct Map {
	typedef T Type;

	Map(const State & state) :
		m_state(state)
	{
		m_data.resize(state.rows * state.cols);
	}

	void fill(const Type & t)
	{
		m_data.assign(m_state.rows * m_state.cols, t);
	}

	Type & operator[](const Location & pos)
	{
		return m_data[pos.row * m_state.cols + pos.col];
	}

	const Type & operator[](const Location & pos) const
	{
		return m_data[pos.row * m_state.cols + pos.col];
	}

private:
	const State & m_state;
	std::vector<Type> m_data;
};

#endif // MAP_H
