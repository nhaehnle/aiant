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

#include "hilldefense.h"

#include <algorithm>
#include <cassert>
#include <limits>

#include "Bot.h"
#include "zoc.h"

using namespace std;

static const uint DefenseEnemyDistance = 20;
static const uint DefensePanicRadius2 = 110;
static const uint MinDefenders = 3;
static const uint EquivalentDangerSlack = 4;

struct HillData {
	struct CandidateDefender {
		CandidateDefender(uint idx, uint d) : antidx(idx), dist(d) {}

		bool operator<(const CandidateDefender & o) const {return dist > o.dist;}

		uint antidx;
		uint dist;
	};

	Location pos;
	Map<uint> dist;
	bool needupdate;
	uint nrdefenders;
	vector<CandidateDefender> candidates;
	vector<uint> defenders;

	HillData(State & state, const Location & p) :
		pos(p),
		dist(state.rows, state.cols),
		needupdate(true)
	{
		dist.fill(0);
	}
};

struct HillDefense::Data {
	vector<HillData *> hills;
	uint lastupdated;
	bool hilldestroyed;

	vector<Location> updatequeue;

	Data() :
		lastupdated(0),
		hilldestroyed(false)
	{
	}

	~Data()
	{
		while (!hills.empty()) {
			delete hills.back();
			hills.pop_back();
		}
	}
};

HillDefense::HillDefense(Bot & b) :
	bot(b),
	d(*new Data),
	state(b.state)
{
}

HillDefense::~HillDefense()
{
	delete &d;
}

void HillDefense::init()
{
}

void HillDefense::update_hill_distances(uint hillidx)
{
	HillData * hd = d.hills[hillidx];

	state.bug << "Update distance map for hill " << hillidx << " at " << hd->pos << endl;

	hd->dist.fill(numeric_limits<uint>::max());

	d.updatequeue.clear();

	hd->dist[hd->pos] = 0;
	d.updatequeue.push_back(hd->pos);

	uint queue_head = 0;
	while (queue_head < d.updatequeue.size()) {
		Location cur = d.updatequeue[queue_head++];
		uint dist = hd->dist[cur];

		const int * dirperm = getdirperm();
		for (int predir = 0; predir < TDIRECTIONS; ++predir) {
			int dir = dirperm[predir];
			Location n = state.getLocation(cur, dir);

			if (state.grid[n.row][n.col].isWater)
				continue;

			if (dist + 1 < hd->dist[n]) {
				hd->dist[n] = dist + 1;
				d.updatequeue.push_back(n);
			}
		}
	}
}

void HillDefense::update_hills()
{
	d.hilldestroyed = false;

	if (state.newwater) {
		for (uint hillidx = 0; hillidx < d.hills.size(); ++hillidx)
			d.hills[hillidx]->needupdate = true;
	}

	// check if hills are still present
	for (uint hillidx = 0; hillidx < d.hills.size(); ++hillidx) {
		HillData * hd = d.hills[hillidx];
		Square & sq = state.grid[hd->pos.row][hd->pos.col];

		if (sq.isVisible && sq.hillPlayer != 0) {
			state.bug << "Lost hill at " << hd->pos << endl;
			delete hd;

			d.hills[hillidx] = d.hills.back();
			d.hills.pop_back();
			hillidx--;
			d.hilldestroyed = true;
		}
	}

	// check if new hills are known now
	for (uint statehillidx = 0; statehillidx < state.myHills.size(); ++statehillidx) {
		Location pos = state.myHills[statehillidx];

		for (uint hillidx = 0; hillidx < d.hills.size(); ++hillidx) {
			if (d.hills[hillidx]->pos == pos)
				goto known;
		}

		state.bug << "Found new hill at " << pos << endl;
		d.hills.push_back(new HillData(state, pos));

	known: ;
	}

	if (!d.hills.size())
		return;

	// update the distance map of one hill
	if (d.lastupdated >= d.hills.size())
		d.lastupdated = 0;
	uint hillidx = d.lastupdated;
	do {
		hillidx++;
		if (hillidx >= d.hills.size())
			hillidx = 0;

		if (d.hills[hillidx]->needupdate) {
			update_hill_distances(hillidx);
			d.lastupdated = hillidx;
			break;
		}
	} while (hillidx != d.lastupdated);
}

void HillDefense::run()
{
	state.bug << "Defense turn " << state.turn << endl;

	update_hills();

	//
	vector<HillData *> defense;

	for (uint hillidx = 0; hillidx < d.hills.size(); ++hillidx) {
		HillData * hd = d.hills[hillidx];
		if (bot.m_zoc.m_enemy[hd->pos] >= DefenseEnemyDistance)
			continue;

		hd->nrdefenders = MinDefenders;
		hd->candidates.clear();
		hd->defenders.clear();

		for (uint enemyantidx = 0; enemyantidx < state.enemyAnts.size(); ++enemyantidx) {
			if (state.eucliddist2(state.enemyAnts[enemyantidx], hd->pos) <= DefensePanicRadius2)
				hd->nrdefenders++;
		}

		for (uint antidx = 0; antidx < bot.m_ants.size(); ++antidx) {
			Ant & ant = bot.m_ants[antidx];
			hd->candidates.push_back(HillData::CandidateDefender(antidx, hd->dist[ant.where]));
		}
		sort(hd->candidates.begin(), hd->candidates.end());

		state.bug << "  defense at " << hd->pos << " requires " << hd->nrdefenders << " defenders" << endl;
		defense.push_back(hd);
	}

	//
	vector<bool> claimed;
	claimed.insert(claimed.begin(), bot.m_ants.size(), false);

	while (!defense.empty()) {
		int besthillidx = -1;
		int bestantidx = -1;
		uint bestdistance = numeric_limits<uint>::max();

		uint p = getprime();
		uint ofs = fastrng() % defense.size();
		for (uint prehillidx = 0; prehillidx < defense.size(); ++prehillidx) {
			uint hillidx = (p * prehillidx + ofs) % defense.size();
			HillData * hd = defense[hillidx];

			while (!hd->candidates.empty()) {
				uint antidx = hd->candidates.back().antidx;
				uint dist = hd->candidates.back().dist;

//				state.bug << "  hill at " << hd->pos << " candidate ant " << antidx
//					<< " at " << bot.m_ants[antidx].where << " dist " << dist << endl;

				if (claimed[antidx]) {
//					state.bug << "    already claimed" << endl;
					hd->candidates.pop_back();
					continue;
				}

				if (dist < bestdistance) {
//					state.bug << "   best so far" << endl;
					bestdistance = dist;
					bestantidx = antidx;
					besthillidx = hillidx;
				}
				break;
			}
		}

		if (besthillidx < 0)
			break;

		HillData * hd = defense[besthillidx];

		hd->nrdefenders--;
		if (hd->nrdefenders == 0) {
			defense[besthillidx] = defense.back();
			defense.pop_back();
		}

		claimed[bestantidx] = true;

		//
		uint hilldanger = bot.m_zoc.m_enemy[hd->pos];
		Ant & ant = bot.m_ants[bestantidx];

		state.bug << "  assign ant " << bestantidx << " at " << ant.where << " to defend " << hd->pos << endl;

		if (ant.assigneddirection)
			continue;

		uint equivdanger = bot.m_zoc.m_enemy[ant.where] + hd->dist[ant.where];
		if (equivdanger <= hilldanger + EquivalentDangerSlack) {
			state.bug << "    equivalent hilldanger " << hilldanger << " is low enough" << endl;
			continue;
		}

		int bestdirection = -1;
		uint bestequivdanger = equivdanger;
		const int * dirperm = getdirperm();
		for (int predir = 0; predir < TDIRECTIONS; ++predir) {
			int dir = dirperm[predir];
			Location n = state.getLocation(ant.where, dir);
			uint newequivdanger = bot.m_zoc.m_enemy[n] + hd->dist[n];
			if (newequivdanger < bestequivdanger) {
				bestequivdanger = newequivdanger;
				bestdirection = dir;
			}
		}

		if (bestdirection == -1) {
			uint bestdist = hd->dist[ant.where];
			for (int predir = 0; predir < TDIRECTIONS; ++predir) {
				int dir = dirperm[predir];
				Location n = state.getLocation(ant.where, dir);
				uint newdist = hd->dist[n];
				if (newdist < bestdist) {
					bestdist = newdist;
					bestdirection = dir;
				}
			}
		}

		//assert(bestdirection != -1);
		state.bug << "    defense move: " << cdir(bestdirection) << endl;

		ant.direction = bestdirection;
		ant.assigneddirection = true;
	}
}

uint HillDefense::getnrhills()
{
	return d.hills.size();
}

const Location & HillDefense::gethill(uint idx)
{
	//assert(idx < d.hills.size());
	return d.hills[idx]->pos;
}

bool HillDefense::hilldestroyed()
{
	return d.hilldestroyed;
}
