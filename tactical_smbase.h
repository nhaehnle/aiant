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

#ifndef TACTICAL_SMBASE_H
#define TACTICAL_SMBASE_H

#include <cassert>
#include <cstdlib>

#include "map.h"
#include "module.h"
#include "State.h"

struct Bot;
struct Submap;

static const uint MaxFoodDist = 16;

struct TacticalSmBase : Module {
	TacticalSmBase(Bot & b);
	virtual ~TacticalSmBase();

	virtual void init();

	void gensubmap(Submap & sm, const Location & center);
	void gensubmap_field(Submap & sm, const Location & local, const Location & global);
	void compute_fooddists();

	struct BaseData;

	Bot & bot;
	State & state;
	BaseData & bd;

	Map<uint32_t> fooddist;

	float ValueKill;
	float ValueLoss;
	float ValueHill;
	float ValueEat;
	float ValueEnemyDist;
	float ValueFoodDist;

	uint MinAggressionAnts;
	float AggressionThresholdShift;
	float AggressionScale;
};

template<typename T>
struct BaseSubmap {
	static const int Radius = 8;
	static const int Size = 2 * Radius + 1;

	T map[Size * Size];

	void fill(T c) {
		for (uint idx = 0; idx < Size * Size; ++idx)
			map[idx] = c;
	}
	static bool inside(const Location & pos) {
		return pos.row >= 0 && pos.row < Size && pos.col >= 0 && pos.col < Size;
	}
	T & operator[](const Location & pos) {
		//assert(inside(pos));
		return map[pos.row * Size + pos.col];
	}
	const T & operator[](const Location & pos) const {
		//assert(inside(pos));
		return map[pos.row * Size + pos.col];
	}
	static bool getneighbour(Location local, int direction, Location & out) {
		out.row = local.row + DIRECTIONS[direction][0];
		out.col = local.col + DIRECTIONS[direction][1];
		return inside(out);
	}
	static bool getneighbouropt(Location local, int direction, Location & out) {
		if (direction >= 0) {
			out.row = local.row + DIRECTIONS[direction][0];
			out.col = local.col + DIRECTIONS[direction][1];
		} else {
			out = local;
		}
		return inside(out);
	}
	static uint manhattandist(const Location & a, const Location & b) {
		return abs(a.row - b.row) + abs(a.col - b.col);
	}

	static uint eucliddist2(const Location & a, const Location & b) {
		int dr = a.row - b.row;
		int dc = a.col - b.col;
		return (dr * dr) + (dc * dc);
	}

	static uint infinitydist(const Location & a, const Location & b) {
		int dr = a.row - b.row;
		int dc = a.col - b.col;
		return std::max(abs(dr), abs(dc));
	}
};

struct Submap : BaseSubmap<uint8_t> {
	static const uint8_t Water = 0x20;
	static const uint8_t Food = 0x10;
	static const uint8_t Hill = 0x08;
	static const uint8_t Ant = 0x04;
	static const uint8_t Mine = 0x02; // must be the two least significant bits
	static const uint8_t Enemy = 0x01;

	void flipsides()
	{
		for (uint idx = 0; idx < Size * Size; ++idx) {
			uint8_t & field = map[idx];

			if ((field ^ (field >> 1)) & 1)
				field ^= Mine | Enemy;
		}
	}
};

std::ostream & operator<<(std::ostream & out, const Submap & sm);

// area in which enemies are counted as attackers
static const int AttackKernelRadius = 2;
static const int AttackKernelSize = 2 * AttackKernelRadius + 1;
static const bool AttackKernel[AttackKernelSize][AttackKernelSize] = {
	{ false, true,  true,  true,  false },
	{ true,  true,  true,  true,  true  },
	{ true,  true,  false, true,  true  },
	{ true,  true,  true,  true,  true  },
	{ false, true,  true,  true,  false },
};

static void add_enemy_kernel(Submap & sm, const Location & center, int delta)
{
	Location kernel;
	for (kernel.row = 0; kernel.row < AttackKernelSize; ++kernel.row) {
		for (kernel.col = 0; kernel.col < AttackKernelSize; ++kernel.col) {
			if (!AttackKernel[kernel.row][kernel.col])
				continue;

			Location other
				((center.row + kernel.row - AttackKernelRadius),
				 (center.col + kernel.col - AttackKernelRadius));
			if (!sm.inside(other))
				continue;

			sm[other] += delta;
		}
	}
}

struct PlayerMove {
	struct AntMove {
		Location pos;
		int8_t direction;
		uint8_t collided : 1; // death by collision
		uint8_t collider : 1; // we ran into the other guy; oops

		uint8_t killed : 1;
		uint8_t killer : 1;
		uint8_t sentoffense : 1;

		AntMove(const Location & p, int8_t dir) :
			pos(p), direction(dir),
			collided(0), collider(0), killed(0), sentoffense(0)
		{}
	};

	std::vector<AntMove> antmoves;
	uint hash;
	uint nrcollided;

	static const uint8_t AntsMask = 0xe0;
	static const uint8_t AntsShift = 5;
	static const uint8_t AttackMask = 0x1f;
	Submap map;

	struct VsOutcome {
		float value;
		int32_t improved;

		VsOutcome(float v, int32_t imp) : value(v), improved(imp) {}
	};

	float worstvalue;
	std::vector<VsOutcome> outcomes;
	float weight;
	float counterweight;

	PlayerMove() {
		reset();
	}

	void reset() {
		hash = 0;
		nrcollided = 0;
		antmoves.clear();
		worstvalue = std::numeric_limits<float>::max();
		outcomes.clear();
		weight = 1.0;
		counterweight = 1.0;
	}

	void init(const PlayerMove & o) {
		antmoves = o.antmoves;
		nrcollided = o.nrcollided;
		map = o.map;
		weight = o.weight;
		counterweight = o.counterweight;
	}

	void computehash() {
		hash = 0;
		for (uint idx = 0; idx < Submap::Size * Submap::Size; ++idx)
			hash = (hash * 312007) + map.map[idx];
	}

	void ant_mark(uint idx) {
		// assume: direction and pos already set
		AntMove & am = antmoves[idx];

		if (!map.inside(am.pos))
			return;

		uint prevnrants = map[am.pos] >> AntsShift;

		map[am.pos] += 1 << AntsShift;
		if (prevnrants == 0) {
			add_enemy_kernel(map, am.pos, 1);
			am.collided = 0;
			am.collider = 0;
		} else if (prevnrants == 1) {
			add_enemy_kernel(map, am.pos, -1);
			for (uint other = 0; other < antmoves.size(); ++other) {
				AntMove & otherm = antmoves[other];
				if (other == idx || am.pos != otherm.pos)
					continue;
				otherm.collided = 1;
			}
			am.collided = 1;
			am.collider = 1;
			nrcollided += 2;
		} else {
			am.collider = 1;
			am.collided = 1;
			nrcollided++;
		}
	}

	void ant_unmark(uint idx) {
		// assume: direction and pos still set
		AntMove & am = antmoves[idx];

		if (!map.inside(am.pos))
			return;

		uint prevnrants = map[am.pos] >> AntsShift;
		//assert(prevnrants >= 1);

		map[am.pos] -= 1 << AntsShift;
		if (prevnrants == 1) {
			add_enemy_kernel(map, am.pos, -1);
		} else if (prevnrants == 2) {
			for (uint other = 0; other < antmoves.size(); ++other) {
				AntMove & otherm = antmoves[other];
				if (other == idx || am.pos != otherm.pos)
					continue;

				add_enemy_kernel(map, am.pos, 1);
				otherm.collided = 0;
				otherm.collider = 0;
				nrcollided -= 2;
				break;
			}
		} else {
			if (!am.collider) {
				for (uint other = 0; other < antmoves.size(); ++other) {
					AntMove & otherm = antmoves[other];
					if (other == idx || am.pos != otherm.pos)
						continue;

					otherm.collider = 0;
					break;
				}
			}
			nrcollided--;
		}
	}

	void reset_killed()
	{
		for (uint idx = 0; idx < antmoves.size(); ++idx) {
			PlayerMove::AntMove & am = antmoves[idx];
			am.sentoffense = 0;
			am.killed = 0;
			am.killer = 0;
		}
	}
};


#endif // TACTICAL_SMBASE_H
