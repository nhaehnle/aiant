
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

ostream & operator<<(ostream & out, const Zoc & zoc)
{
	Location cur;
	for (cur.row = 0; cur.row < zoc.state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < zoc.state.cols; ++cur.col) {
			const Square & sq = zoc.state.grid[cur.row][cur.col];
			if (sq.isWater)
				out << '%';
			else if (sq.hillPlayer >= 0)
				out << (char)('A' + sq.hillPlayer);
			else if (sq.ant >= 0)
				out << (char)('a' + sq.ant);
			else if (sq.isFood)
				out << '*';
			else if (zoc.m_enemy[cur] == 0)
				out << '@';
			else {
				int diff = zoc.m_enemy[cur] - zoc.m_me[cur];
				if (diff > 3)
					out << ' ';
				else if (diff > 0)
					out << '.';
				else
					out << '+';
			}
		}
		out << endl;
	}

	return out;
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
		if (m_enemy[cur] > 0) {
			m_enemy[cur] = 0;
			open.push(cur);
		}
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

	state.bug << "ZOC data turn " << state.turn << endl;
	state.bug << *this;
}
