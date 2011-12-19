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

#include "scout.h"

#include <algorithm>
#include <list>
#include <queue>

#include "Bot.h"
#include "map.h"
#include "zoc.h"
#include <limits>

using namespace std;

static const uint ScoutFracNom = 1;
static const uint ScoutFracDenom = 50;
static const uint ScoutCountMargin = 10;

static const int ScoutMinAge = 50;

struct FoundInvisbleArea {
	vector<Location> where;
	uint oldest;
	bool used;
};

struct Scout::Data {
	uint nrscouts;
	Map<int> scoutmap;

	struct LocationHeight {
		LocationHeight(const Location & p, int h) : pos(p), height(h) {}

		bool operator<(const LocationHeight & other) const {return height < other.height;}

		Location pos;
		int height;
	};
	priority_queue<LocationHeight> updatequeue;

	Data() : nrscouts(0)
	{
	}
};

Scout::Scout(Bot & b) :
	d(*new Data),
	bot(b),
	state(b.state)
{
}

Scout::~Scout()
{
	delete &d;
}

void Scout::init()
{
	d.scoutmap.resize(state.rows, state.cols);
}

void Scout::update_nrscouts()
{
	uint upper = 1 + (bot.m_ants.size() * ScoutFracNom) / ScoutFracDenom;
	uint lower;
	if (bot.m_ants.size() < ScoutCountMargin)
		lower = 0;
	else
		lower = 1 + ((bot.m_ants.size() - ScoutCountMargin) * ScoutFracNom) / ScoutFracDenom;

	state.bug << "previous scouts: " << d.nrscouts << " lower " << lower << " upper " << upper;
	d.nrscouts = max(lower, min(upper, d.nrscouts));
	state.bug << " new: " << d.nrscouts << endl;
}

struct ShowScoutMap {
	ShowScoutMap(Scout & scout_) :
		scout(scout_),
		d(scout_.d),
		bot(scout_.bot),
		state(scout_.state)
	{}

	void print(ostream & out) const
	{
		Location cur;
		for (cur.row = 0; cur.row < state.rows; ++cur.row) {
			for (cur.col = 0; cur.col < state.cols; ++cur.col) {
				Square & sq = state.grid[cur.row][cur.col];
				if (sq.hillPlayer >= 0)
					out << (char)('A' + sq.hillPlayer);
				else if (sq.ant >= 0)
					out << (char)('a' + sq.ant);
				else if (sq.isFood)
					out << '*';
				else if (sq.isWater)
					out << '%';
				else {
					if (bot.m_zoc.m_enemy[cur] < bot.m_zoc.m_me[cur]) {
						if (sq.isVisible)
							out << '+';
						else
							out << '@';
					} else {
						uint increases = 0;
						for (int dir = 0; dir < TDIRECTIONS; ++dir) {
							Location n = state.getLocation(cur, dir);
							if (d.scoutmap[n] > d.scoutmap[cur])
								increases |= (1 << dir);
						}
						if (increases == 1)
							out << '^';
						else if (increases == 2)
							out << '>';
						else if (increases == 4)
							out << 'v';
						else if (increases == 8)
							out << '<';
						else if (increases == 0)
							out << 'o';
						else
							out << ' ';
					}
				}
			}
			out << endl;
		}
	}

	Scout & scout;
	Scout::Data & d;
	Bot & bot;
	State & state;
};

ostream & operator<<(ostream & out, const ShowScoutMap & sm)
{
	sm.print(out);
	return out;
}

void Scout::recompute_maps()
{
	d.scoutmap.fill(numeric_limits<int>::min());

	Location cur;
	for (cur.row = 0; cur.row < state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < state.cols; ++cur.col) {
			Square & sq = state.grid[cur.row][cur.col];
			if (sq.isWater || sq.isVisible)
				continue;
			if (state.turn - sq.lastseen < ScoutMinAge)
				continue;
			if (bot.m_zoc.m_enemy[cur] < bot.m_zoc.m_me[cur])
				continue;

			d.scoutmap[cur] = state.turn - sq.lastseen;

			d.updatequeue.push(Data::LocationHeight(cur, d.scoutmap[cur]));
		}
	}

	while (!d.updatequeue.empty()) {
		cur = d.updatequeue.top().pos;
		d.updatequeue.pop();

		int height = d.scoutmap[cur] - 1;
		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n = state.getLocation(cur, dir);
			Square & sq = state.grid[n.row][n.col];
			if (sq.isWater || bot.m_zoc.m_enemy[n] < bot.m_zoc.m_me[n])
				continue;
			if (height > d.scoutmap[n]) {
				d.scoutmap[n] = height;
				d.updatequeue.push(Data::LocationHeight(n, height));
			}
		}
	}

	//state.bug << ShowScoutMap(*this);
}

struct SortAntByHeight {
	SortAntByHeight(Bot & b, Map<int> & m) : bot(b), map(m) {}

	bool operator()(uint a, uint b) const {
		return map[bot.m_ants[a].where] > map[bot.m_ants[b].where];
	}

	Bot & bot;
	Map<int> & map;
};

void Scout::run()
{
	state.bug << "Scout turn " << state.turn << endl;

	update_nrscouts();
	if (!d.nrscouts)
		return;

	recompute_maps();

	vector<uint> ants;
	ants.reserve(bot.m_ants.size());
	for (uint idx = 0; idx < bot.m_ants.size(); ++idx)
		ants.push_back(idx);

	sort(ants.begin(), ants.end(), SortAntByHeight(bot, d.scoutmap));

	for (uint i = 0; i < d.nrscouts; ++i) {
		Ant & ant = bot.m_ants[ants[i]];
		if (ant.assigneddirection)
			continue;

		int my = d.scoutmap[ant.where];
		if (my < 0)
			break;

		const int * dirperm = getdirperm();
		for (int predir = 0; predir < TDIRECTIONS; ++predir) {
			int dir = dirperm[predir];
			Location n = state.getLocation(ant.where, dir);
			if (d.scoutmap[n] > my) {
				state.bug << "Scout " << i << "  ant " << ants[i] << " at " << ant.where << " move " << cdir(dir) << endl;
				ant.direction = dir;
				ant.assigneddirection = true;
				break;
			}
		}
	}
}
