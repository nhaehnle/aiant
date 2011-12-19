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

#include <cassert>
#include <map>

#include "Bot.h"
#include "astar.h"
#include "diffusion.h"
#include "foodseeker.h"
#include "hilldefense.h"
#include "offense.h"
#include "opportunisticattack.h"
#include "scout.h"
#include "symmetry.h"
#include "tactical_sm.h"
#include "zoc.h"

using namespace std;

static const uint KeepDirpermTurns = 64;
static const float SafetyMarginMs = 40.0;
static const float SafetyMult = 5.0;

void Ant::reset()
{
	assigneddirection = false;
	direction = -1;
	if (!dirperm || (fastrng() % KeepDirpermTurns) == 0)
		dirperm = getdirperm();
}


struct Bot::Data {
	struct Argument {
		std::string value;
		bool used;

		Argument(const std::string & s) : value(s), used(false) {}
	};
	std::map<std::string, Argument> args;

	double turntime;
	double lastchecktime;
	double checkavg[4];
	Map<bool> claimed;
};

//constructor
Bot::Bot(int argc, char * * argv) :
	d(*new Data),
	m_zoc(*new Zoc(state)),
	m_symmetry(*new SymmetryFinder(*this)),
	m_foodseeker(*new FoodSeeker(*this)),
	m_scout(*new Scout(*this)),
	m_hilldefense(*new HillDefense(*this)),
	m_opportunisticattack(*new OpportunisticAttack(*this)),
	m_diffusion(*new Diffusion(*this)),
	m_offense(*new Offense(*this)),
	m_tactical(*new TacticalSm(*this))
{
	for (int i = 1; i < argc; i += 2) {
		if (argv[i][0] != '-' || argv[i][1] != '-') {
			cerr << "Argument must start with '--'" << endl;
			abort();
		}

		if (i + 1 >= argc) {
			cerr << "Argument missing" << endl;
			abort();
		}

		string key = argv[i] + 2;
		string value = argv[i+1];

		state.bug.time << "Command line argument: " << key << " = " << value << endl;
		d.args.insert(make_pair(key, Data::Argument(value)));
	}

	d.checkavg[0] = d.checkavg[1] = d.checkavg[2] = d.checkavg[3] = 1.0;
}

Bot::~Bot()
{
	delete &m_tactical;
	delete &m_offense;
	delete &m_diffusion;
	delete &m_opportunisticattack;
	delete &m_hilldefense;
	delete &m_scout;
	delete &m_foodseeker;
	delete &m_symmetry;
	delete &m_zoc;
	delete &d;
}

float Bot::getargfloat(const std::string & key, float def)
{
	map<string, Data::Argument>::iterator it = d.args.find(key);
	if (it == d.args.end())
		return def;

	it->second.used = true;
	return strtof(it->second.value.c_str(), 0);
}

bool Bot::timeover()
{
	double t = state.timer.getTime();
	if
		((t + SafetyMarginMs + SafetyMult * d.checkavg[0] > d.turntime) ||
		 (t + SafetyMarginMs + SafetyMult * d.checkavg[1] > d.turntime) ||
		 (t + SafetyMarginMs + SafetyMult * d.checkavg[2] > d.turntime) ||
		 (t + SafetyMarginMs + SafetyMult * d.checkavg[3] > d.turntime))
		return true;

	if (d.lastchecktime >= 0.0) {
		double delta = t - d.lastchecktime;
		d.checkavg[0] = (d.checkavg[0] * 0.66) + (delta * 0.34);
		d.checkavg[1] = (d.checkavg[1] * 0.9) + (delta * 0.1);
		d.checkavg[2] = (d.checkavg[2] * 0.99) + (delta * 0.01);
		d.checkavg[3] = (d.checkavg[3] * 0.999) + (delta * 0.001);
	}

	d.lastchecktime = t;
	return false;
}

//plays a single game of Ants.
void Bot::playGame()
{
	//reads the game parameters and sets up
	cin >> state;
	state.setup();
	endTurn();

	d.claimed.resize(state.rows, state.cols);
	d.turntime = getargfloat("turntime", state.turntime);

	//
	m_zoc.init();
	m_symmetry.init();
	m_foodseeker.init();
	m_scout.init();
	m_hilldefense.init();
	m_opportunisticattack.init();
	m_diffusion.init();
	m_offense.init();
	m_tactical.init();

	//
	for (map<string, Data::Argument>::const_iterator it = d.args.begin(); it != d.args.end(); ++it) {
		if (!it->second.used) {
			cerr << "Argument " << it->first << " not used" << endl;
			abort();
		}
	}

	//continues making moves while the game is not over
	while(cin >> state)
	{
		state.updateVisionInformation();
		makeMoves();
		endTurn();
	}
}

uint Bot::myantidx_at(const Location & pos)
{
	//assert(state.grid[pos.row][pos.col].ant == 0);

	for (uint idx = 0; idx < m_ants.size(); ++idx) {
		if (m_ants[idx].where == pos)
			return idx;
	}

	abort();
}

bool Bot::try_rotate_move(uint antidx)
{
	Ant & ant = m_ants[antidx];
	int altdir = (ant.direction + 1 + (fastrng() & 2)) % TDIRECTIONS;
	Location moveto = state.getLocation(ant.where, altdir);

	state.bug << "  try altdir " << cdir(altdir) << " to " << moveto << endl;

	if (!d.claimed[moveto] && !state.grid[moveto.row][moveto.col].isWater) {
		Location next = state.getLocation(moveto, ant.direction);
		if (!state.grid[next.row][next.col].isWater) {
			ant.direction = altdir;
			state.bug << "  ok" << endl;
			return true;
		}
	}

	state.bug << "  nope" << endl;

	return false;
}

void Bot::make_moves()
{
	for (uint antidx = 0; antidx < m_ants.size(); ++antidx) {
		Ant & ant = m_ants[antidx];
		Location moveto;

		moveto = ant.where;
		if (ant.direction >= 0)
			moveto = state.getLocation(moveto, ant.direction);

		state.bug << "ant at " << ant.where << " initial attempt " << moveto << " dir " << cdir(ant.direction) << endl;

		if (d.claimed[moveto] && (ant.direction < 0 || !try_rotate_move(antidx))) {
			if (!d.claimed[ant.where]) {
				state.bug << "  current pos unclaimed" << endl;
				ant.direction = -1;
			} else {
				state.bug << "  rotate" << endl;
				const int * dirperm = getdirperm();
				for (int predir = 0; predir < TDIRECTIONS; ++predir) {
					int dir = dirperm[predir];
					moveto = state.getLocation(ant.where, dir);
					if (!d.claimed[moveto]) {
						ant.direction = dir;
						break;
					}
				}
			}
		}

		state.bug << "  move " << cdir(ant.direction) << endl;

		if (ant.direction >= 0) {
			moveto = state.getLocation(ant.where, ant.direction);
			d.claimed[moveto] = true;
			state.makeMove(ant.where, ant.direction);
		} else {
			d.claimed[ant.where] = true;
		}
	}
}

void Bot::update_ants()
{
	d.claimed.fill(false);

	for (uint antidx = 0; antidx < m_ants.size(); ++antidx) {
		Ant & ant = m_ants[antidx];
		if (ant.direction >= 0)
			ant.where = state.getLocation(ant.where, ant.direction);

		if (state.grid[ant.where.row][ant.where.col].ant != 0 || d.claimed[ant.where]) {
			state.bug << "ant at " << ant.where << " move dir " << cdir(ant.direction) << " disappeared" << endl;
			m_ants[antidx] = m_ants.back();
			m_ants.pop_back();
			antidx--;
		} else {
			state.bug << "ant at " << ant.where << " remembered" << endl;
			d.claimed[ant.where] = true;

			ant.reset();
		}
	}

	for (uint antidx = 0; antidx < state.myAnts.size(); ++antidx) {
		const Location & pos = state.myAnts[antidx];
		if (!d.claimed[pos]) {
			state.bug << "new ant at " << pos << endl;
			d.claimed[pos] = true;
			m_ants.push_back(Ant());
			Ant & ant = m_ants.back();
			ant.where = pos;
		}
	}
}

//makes the bots moves for the turn
void Bot::makeMoves()
{
	state.bug.time << "turn " << state.turn << ":" << endl;
	state.bug << state << endl;

	d.lastchecktime = -1.0;

	m_tactical.learn();

	m_zoc.update();
	//m_symmetry.run();

	//
	update_ants();

	d.claimed.fill(false);

	for (uint foodidx = 0; foodidx < state.food.size(); ++foodidx)
		d.claimed[state.food[foodidx]] = true;

	state.bug.time << "time after maintenance: " << state.timer.getTime() << "ms" << endl;

	//
	m_foodseeker.run();

	state.bug.time << "time after food: " << state.timer.getTime() << "ms" << endl;

	//
	m_hilldefense.run();

	//
	m_opportunisticattack.run();

	//
	m_scout.run();

	state.bug.time << "time before diffusion: " << state.timer.getTime() << "ms" << endl;

	//
	m_diffusion.run();

	//
	//m_offense.run();

	//
	for (uint antidx = 0; antidx < m_ants.size(); ++antidx) {
		Ant & ant = m_ants[antidx];

		state.bug << "ant " << antidx << " at " << ant.where << " dir " << cdir(ant.direction) << endl;

		if (ant.direction >= 0) {
			Location n = state.getLocation(ant.where, ant.direction);
			if (state.grid[n.row][n.col].isWater) {
				state.bug << "  attempt to walk into water" << endl;
				ant.direction = -1;
				ant.assigneddirection = false;
			}
		}

		if (!ant.assigneddirection) {
			// not looking for food, go towards enemy territory
			uint my = m_zoc.m_enemy[ant.where];
			const int * dirperm = ant.dirperm;
			for (int predir = 0; predir < TDIRECTIONS; predir++) {
				int dir = dirperm[predir];
				Location loc = state.getLocation(ant.where, dir);

				if (!state.grid[loc.row][loc.col].isWater) {
					if (m_zoc.m_enemy[loc] < my) {
						state.bug << "  zoc move " << cdir(dir) << endl;
						ant.direction = dir;
						ant.assigneddirection = true;
						break;
					}
				}
			}
		}
	}

	state.bug.time << "time before tactical: " << state.timer.getTime() << "ms" << endl;

	//
	m_tactical.run();

	//
	make_moves();

	state.bug.time << "time taken in turn " << state.turn << ": " << state.timer.getTime() << "ms, checkavg: "
		<< d.checkavg[0] << " " << d.checkavg[1] << " " << d.checkavg[2] << " " << d.checkavg[3] << endl << endl;
}

//finishes the turn
void Bot::endTurn()
{
	if(state.turn > 0)
		state.reset();
	state.turn++;

	cout << "go" << endl;
}
