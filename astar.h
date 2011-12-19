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

#ifndef ASTAR_H
#define ASTAR_H

#include <cassert>

#include "State.h"
#include "map.h"

struct BaseAStar {
	BaseAStar(State & state) :
		m_state(state),
		m_map(state.rows, state.cols),
		m_queue(MyCompare(*this))
	{
	}

	void clear() {
		while (!m_queue.empty())
			m_queue.pop();
		m_map.fill(Field());
	}

	Location getsource(const Location & to) const {
		//assert(m_map[to].state == Field::WHITE);

		Location cur = to;

		int backlink;
		while ((backlink = m_map[cur].backlink) >= 0)
			cur = m_state.getLocation(cur, reversedir(backlink));

		return cur;
	}

	int getlaststep(const Location & to) const {
		//assert(m_map[to].state == Field::WHITE);
		return m_map[to].backlink;
	}

protected:
	struct Field {
		enum State {
			BLACK = 0,
			GREY,
			WHITE
		};

		State state;
		int32_t cost;
		int32_t estimate;
		int32_t backlink;

		Field() : state(BLACK), backlink(-1) {}
	};

	typedef std::pair<Location, int32_t> QueuePair;

	struct MyCompare {
		MyCompare(const BaseAStar & astar) : m_astar(astar) {}

		bool operator()(const QueuePair & a, const QueuePair & b) const {
			return a.second > b.second;
		}

		const BaseAStar & m_astar;
	};

	State & m_state;
	Map<Field> m_map;
	std::priority_queue<Location, std::vector<QueuePair>, MyCompare> m_queue;
};

struct LocationEvalZero {
	int32_t operator()(const Location & pos) const {return 0;}
};

struct LocationEvalSeek {
	LocationEvalSeek(State & s,const Location & t) : state(s), target(t) {}

	int32_t operator()(const Location & pos) const {
		return state.manhattanDistance(pos, target);
	}

	State & state;
	Location target;
};

struct StepEvalDefault {
	StepEvalDefault(const State & state) : m_state(state) {}

	bool operator()(const Location & from, const Location & to, int direction, int32_t * pstepcost) const {
		if (m_state.grid[to.row][to.col].isWater)
			return false;

		*pstepcost = 1;
		return true;
	}

	const State & m_state;
};

template<typename LocationEval, typename StepEval>
struct AStar : BaseAStar {
	AStar(State & state, const LocationEval & loceval, const StepEval & stepeval) :
		BaseAStar(state),
		m_loceval(loceval),
		m_stepeval(stepeval)
	{
	}

	void push(const Location & pos, int32_t cost = 0, int32_t backlink = -1)
	{
		Field & f = m_map[pos];
		if (f.state == Field::BLACK) {
			f.state = Field::GREY;
			f.cost = cost;
			f.estimate = m_loceval(pos);
			f.backlink = backlink;
			m_queue.push(std::make_pair(pos, f.cost + f.estimate));
		} else if (f.state == Field::GREY && cost < f.cost) {
			f.cost = cost;
			f.backlink = backlink;
			m_queue.push(std::make_pair(pos, f.cost + f.estimate));
		}
	}

	bool step(Location & pos, int32_t & cost)
	{
		do {
			if (m_queue.empty())
				return false;
			pos = m_queue.top().first;
			m_queue.pop();
		} while (m_map[pos].state != Field::GREY);

		Field & f = m_map[pos];
		cost = f.cost;
		f.state = Field::WHITE;

		const int * dirperm = getdirperm();
		for (int predir = 0; predir < TDIRECTIONS; ++predir) {
			int dir = dirperm[predir];
			Location neighbour = m_state.getLocation(pos, dir);
			int32_t stepcost;

			Field & fn = m_map[neighbour];
			if (fn.state == Field::BLACK || fn.state == Field::GREY) {
				if (m_stepeval(pos, neighbour, dir, &stepcost))
					push(neighbour, cost + stepcost, dir);
			}
		}

		return true;
	}

private:
	LocationEval m_loceval;
	StepEval m_stepeval;
};

#endif // ASTAR_H
