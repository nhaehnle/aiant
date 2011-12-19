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

#include "offense.h"

#include <algorithm>
#include <limits>

#include "Bot.h"
#include "hilldefense.h"
#include "symmetry.h"

using namespace std;

static const uint MinUnassignedAnts = 60;
static const uint AttackAntsMargin = 10;
static const uint AttackAntNom = 1;
static const uint AttackAntDenom = 2;

static const uint AttackReconsiderDist = 12;

struct EnemyHillData {
	Location pos;
	Map<uint> dist;
	uint distance_to_own_hill;
	bool needupdate;

	EnemyHillData(State & state, const Location & p) :
		pos(p),
		dist(state.rows, state.cols),
		distance_to_own_hill(numeric_limits<uint>::max()),
		needupdate(true)
	{
		dist.fill(numeric_limits<uint>::max());
	}
};

struct Offense::Data {
	vector<EnemyHillData *> hills;
	uint lastupdated;
	uint nrattackants;
	int attackhill;

	vector<Location> updatequeue;

	Data() :
		lastupdated(0),
		nrattackants(0),
		attackhill(-1)
	{
	}
	~Data() {
		while (!hills.empty()) {
			delete hills.back();
			hills.pop_back();
		}
	}
};

Offense::Offense(Bot & b) :
	bot(b),
	state(b.state),
	d(*new Data)
{
}

Offense::~Offense()
{
}

void Offense::init()
{

}

void Offense::update_hill_distance(uint hillidx)
{
	EnemyHillData * ehd = d.hills[hillidx];

	state.bug << "  update distance map for enemy hill " << hillidx << " at " << ehd->pos << endl;
	ehd->needupdate = false;

	ehd->dist.fill(numeric_limits<uint>::max());
	d.updatequeue.clear();

	uint queue_head = 0;

	ehd->dist[ehd->pos] = 0;
	d.updatequeue.push_back(ehd->pos);

	while (queue_head < d.updatequeue.size()) {
		const Location cur = d.updatequeue[queue_head++];
		uint dist = ehd->dist[cur] + 1;

		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n = state.getLocation(cur, dir);
			if (bot.m_symmetry.map[n] & SymmetryFinder::MapWater)
				continue;
			if (dist >= ehd->dist[n])
				continue;

			ehd->dist[n] = dist;
			d.updatequeue.push_back(n);
		}
	}

	ehd->distance_to_own_hill = numeric_limits<uint>::max();

	uint nrmyhills = bot.m_hilldefense.getnrhills();
	for (uint idx = 0; idx < nrmyhills; ++idx) {
		const Location & myhillpos = bot.m_hilldefense.gethill(idx);
		ehd->distance_to_own_hill = min(ehd->distance_to_own_hill, ehd->dist[myhillpos]);
	}
}


void Offense::update_hills()
{
	if (bot.m_symmetry.newwater || bot.m_hilldefense.hilldestroyed()) {
		for (uint hillidx = 0; hillidx < d.hills.size(); ++hillidx)
			d.hills[hillidx]->needupdate = true;
	}

	//
	for (uint hillidx = 0; hillidx < d.hills.size(); ++hillidx) {
		EnemyHillData * ehd = d.hills[hillidx];
		Square & sq = state.grid[ehd->pos.row][ehd->pos.col];

		if (sq.isVisible && sq.hillPlayer < 0) {
			state.bug << "  hill destroyed at " << ehd->pos << endl;

			if (d.attackhill == (int)hillidx) {
				d.attackhill = -1;
			} else if (d.attackhill == (int)(d.hills.size() - 1)) {
				d.attackhill = hillidx;
			}

			delete ehd;
			d.hills[hillidx] = d.hills.back();
			d.hills.pop_back();
			hillidx--;
		}
	}

	for (uint stidx = 0; stidx < bot.m_symmetry.enemy_hills.size(); ++stidx) {
		const Location & pos = bot.m_symmetry.enemy_hills[stidx];

		for (uint hillidx = 0; hillidx < d.hills.size(); ++hillidx) {
			if (d.hills[hillidx]->pos == pos)
				goto found;
		}

		state.bug << "  found new enemy hill at " << pos << endl;
		d.hills.push_back(new EnemyHillData(state, pos));

	found: ;
	}

	if (d.hills.empty())
		return;

	//
	if (d.lastupdated >= d.hills.size())
		d.lastupdated = 0;
	uint hillidx = d.lastupdated;
	do {
		hillidx++;
		if (hillidx >= d.hills.size())
			hillidx = 0;

		if (d.hills[hillidx]->needupdate) {
			update_hill_distance(hillidx);
			d.lastupdated = hillidx;
			break;
		}
	} while (hillidx != d.lastupdated);
}

uint Offense::closesthill()
{
	uint closest = 0;
	uint bestdist = numeric_limits<uint>::max();

	for (uint idx = 0; idx < d.hills.size(); ++idx) {
		if (d.hills[idx]->distance_to_own_hill <= bestdist) {
			bestdist = d.hills[idx]->distance_to_own_hill;
			closest = idx;
		}
	}

	return closest;
}

struct AttackCandidate {
	AttackCandidate(uint idx) : antidx(idx), dist(0) {}

	bool operator<(const AttackCandidate & o) const {return dist < o.dist;}

	uint antidx;
	uint dist;
};

void Offense::run()
{
	state.bug << "Offense turn " << state.turn << endl;

	update_hills();

	if (d.hills.empty())
		return;

	//
	vector<AttackCandidate> unassignedants;

	for (uint antidx = 0; antidx < bot.m_ants.size(); ++antidx) {
		if (!bot.m_ants[antidx].assigneddirection)
			unassignedants.push_back(AttackCandidate(antidx));
	}

	//
	uint minants = 0;
	uint maxants = 0;

	if (unassignedants.size() >= MinUnassignedAnts)
		minants = (unassignedants.size() * AttackAntNom) / AttackAntDenom;
	if (unassignedants.size() + AttackAntsMargin >= MinUnassignedAnts)
		maxants = ((unassignedants.size() + AttackAntsMargin) * AttackAntNom) / AttackAntDenom;

	d.nrattackants = min(maxants, max(minants, d.nrattackants));
	d.nrattackants = min(d.nrattackants, (uint)unassignedants.size());

	if (!d.nrattackants) {
		d.attackhill = -1;
		return;
	}

	//
	if (d.attackhill < 0) {
		d.attackhill = closesthill();
		state.bug << "  start attack on " << d.hills[d.attackhill]->pos << endl;
	}

	//
	EnemyHillData * ehd = d.hills[d.attackhill];
	for (uint idx = 0; idx < unassignedants.size(); ++idx) {
		Ant & ant = bot.m_ants[unassignedants[idx].antidx];
		unassignedants[idx].dist = ehd->dist[ant.where];
	}

	sort(unassignedants.begin(), unassignedants.end());

	if (unassignedants.front().dist >= AttackReconsiderDist) {
		uint closest = closesthill();
		if ((int)closest != d.attackhill) {
			state.bug << "  closest hill changed, rerouting attack" << endl;
			d.attackhill = closest;
			return; // until next turn
		}
	}

	//
	state.bug << "  attacking " << ehd->pos << " with " << d.nrattackants << " ants" << endl;
	for (uint idx = 0; idx < d.nrattackants; ++idx) {
		Ant & ant = bot.m_ants[unassignedants[idx].antidx];

		int bestdir = -1;
		uint bestdist = ehd->dist[ant.where];
		const int * dirperm = getdirperm();
		for (int predir = 0; predir < TDIRECTIONS; ++predir) {
			int dir = dirperm[predir];
			Location n = state.getLocation(ant.where, dir);
			uint dist = ehd->dist[n];

			if (dist < bestdist) {
				bestdir = dir;
				bestdist = dist;
			}
		}

		state.bug << "  ant " << unassignedants[idx].antidx << " at " << ant.where << " to " << cdir(bestdir) << endl;
		ant.direction = bestdir;
		ant.assigneddirection = true;
	}
}
