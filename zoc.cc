/*
 * Copyright 2011 Nicolai Hähnle
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#include "zoc.h"
#include <limits>

#include "State.h"

using namespace std;

struct Zoc::Data {
	vector<Location> oldenemies;
	vector<Location> queue;
};

Zoc::Zoc(State & state_) :
	state(state_),
	d(*new Data)
{
}

Zoc::~Zoc()
{
	delete &d;
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
	d.oldenemies.clear();
	d.queue.clear();

	Location cur;
	for (cur.row = 0; cur.row < state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < state.cols; ++cur.col) {
			Square & sq = state.grid[cur.row][cur.col];
			if (!sq.isVisible && m_enemy[cur] == 0)
				d.oldenemies.push_back(cur);
		}
	}

	m_enemy.fill(std::numeric_limits<uint>::max());

	for (uint idx = 0; idx < d.oldenemies.size(); ++idx) {
		cur = d.oldenemies[idx];
		if (m_enemy[cur] > 0) {
			m_enemy[cur] = 0;
			d.queue.push_back(cur);
		}
		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n = state.getLocation(cur, dir);
			Square & sq = state.grid[n.row][n.col];
			if (sq.isWater)
				continue;
			if (!sq.isVisible && m_enemy[n] != 0) {
				m_enemy[n] = 0;
				d.queue.push_back(n);
			}
		}
	}

	for (uint idx = 0; idx < state.enemyAnts.size(); ++idx) {
		cur = state.enemyAnts[idx];
		if (m_enemy[cur] > 0) {
			m_enemy[cur] = 0;
			d.queue.push_back(cur);
		}
	}

	for (uint idx = 0; idx < state.enemyHills.size(); ++idx) {
		cur = state.enemyHills[idx];
		if (m_enemy[cur] > 0) {
			m_enemy[cur] = 0;
			d.queue.push_back(cur);
		}
	}

	uint queue_head = 0;
	while (queue_head < d.queue.size()) {
		cur = d.queue[queue_head++];

		uint dist = m_enemy[cur];

		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n = state.getLocation(cur, dir);
			if (state.grid[n.row][n.col].isWater)
				continue;
			if (m_enemy[n] > dist+1) {
				m_enemy[n] = dist+1;
				d.queue.push_back(n);
			}
		}
	}

	m_me.fill(std::numeric_limits<uint>::max());
	d.queue.clear();
	queue_head = 0;

	for (uint idx = 0; idx < state.myAnts.size(); ++idx) {
		cur = state.myAnts[idx];
		if (m_me[cur] > 0) {
			m_me[cur] = 0;
			d.queue.push_back(cur);
		}
	}

	while (queue_head < d.queue.size()) {
		cur = d.queue[queue_head++];

		uint dist = m_me[cur];

		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n = state.getLocation(cur, dir);
			if (state.grid[n.row][n.col].isWater)
				continue;
			if (m_me[n] > dist+1) {
				m_me[n] = dist+1;
				d.queue.push_back(n);
			}
		}
	}

	state.bug << "ZOC data turn " << state.turn << endl;
	state.bug << *this;
}
