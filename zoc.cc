
#include "zoc.h"

#include <cassert>
#include <limits>

#include "State.h"

using namespace std;

Zoc::Zoc(State & state_) :
	state(state_)
{
}

Zoc::~Zoc()
{
}

void Zoc::init()
{
	m_me.resize(state.rows, state.cols);
	m_enemy.resize(state.rows, state.cols);
	m_enemy.fill(0);
}

void Zoc::update()
{
	std::vector<Location> oldenemies;

	Location cur;
	for (cur.row = 0; cur.row < state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < state.cols; ++cur.col) {
			Square & sq = state.grid[cur.row][cur.col];
			if (!sq.isVisible && m_enemy[cur] == 0)
				oldenemies.push_back(cur);
		}
	}

	m_enemy.fill(std::numeric_limits<uint>::max());

	std::queue<Location> open;
	for (uint idx = 0; idx < oldenemies.size(); ++idx) {
		cur = oldenemies[idx];
		m_enemy[cur] = 0;
		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n = state.getLocation(cur, dir);
			Square & sq = state.grid[n.row][n.col];
			if (sq.isWater)
				continue;
			if (!sq.isVisible && m_enemy[n] != 0) {
				m_enemy[n] = 0;
				open.push(n);
			}
		}
	}

	for (uint idx = 0; idx < state.enemyAnts.size(); ++idx) {
		cur = state.enemyAnts[idx];
		if (m_enemy[cur] > 0) {
			m_enemy[cur] = 0;
			open.push(cur);
		}
	}

	for (uint idx = 0; idx < state.enemyHills.size(); ++idx) {
		cur = state.enemyHills[idx];
		if (m_enemy[cur] > 0) {
			m_enemy[cur] = 0;
			open.push(cur);
		}
	}

	while (!open.empty()) {
		cur = open.front();
		open.pop();

		uint dist = m_enemy[cur];

		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n = state.getLocation(cur, dir);
			if (state.grid[n.row][n.col].isWater)
				continue;
			if (m_enemy[n] > dist+1) {
				m_enemy[n] = dist+1;
				open.push(n);
			}
		}
	}

	m_me.fill(std::numeric_limits<uint>::max());

	for (uint idx = 0; idx < state.myAnts.size(); ++idx) {
		cur = state.myAnts[idx];
		if (m_me[cur] > 0) {
			m_me[cur] = 0;
			open.push(cur);
		}
	}

	while (!open.empty()) {
		cur = open.front();
		open.pop();

		uint dist = m_me[cur];

		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n = state.getLocation(cur, dir);
			if (state.grid[n.row][n.col].isWater)
				continue;
			if (m_me[n] > dist+1) {
				m_me[n] = dist+1;
				open.push(n);
			}
		}
	}

	//print();
}

void Zoc::print()
{
	state.bug << "ZOC data" << endl;

	Location cur;
	for (cur.row = 0; cur.row < state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < state.cols; ++cur.col) {
			if (state.grid[cur.row][cur.col].isWater) {
				state.bug << '%';
				continue;
			}
			if (m_enemy[cur] == 0) {
				state.bug << '@';
			} else {
				int diff = m_enemy[cur] - m_me[cur];
				if (diff > 5)
					state.bug << ' ';
				else if (diff > 0)
					state.bug << '.';
				else if (diff > -2)
					state.bug << '+';
				else
					state.bug << '*';
			}
		}
		state.bug << endl;
	}
}
