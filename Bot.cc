#include <cassert>

#include "Bot.h"
#include "astar.h"
#include "hilldefense.h"
#include "scout.h"
#include "tactical.h"
#include "zoc.h"

using namespace std;

static const unsigned int TacticalClose = 2;
static const unsigned int TacticalMid = 3;
static const unsigned int TacticalFar = 5;


//constructor
Bot::Bot() :
	m_zoc(*new Zoc(state)),
	m_scout(*new Scout(*this)),
	m_hilldefense(*new HillDefense(*this))
{

}

Bot::~Bot()
{
	delete &m_hilldefense;
	delete &m_scout;
	delete &m_zoc;
}

//plays a single game of Ants.
void Bot::playGame()
{
	//reads the game parameters and sets up
	cin >> state;
	state.setup();
	endTurn();

	m_zoc.init();
	m_scout.init();
	m_hilldefense.init();

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
	assert(state.grid[pos.row][pos.col].ant == 0);

	for (uint idx = 0; idx < m_ants.size(); ++idx) {
		if (m_ants[idx].where == pos)
			return idx;
	}

	abort();
}

uint Bot::foodidx_at(const Location & pos)
{
	state.bug << "foodidx " << pos << endl;

	assert(state.grid[pos.row][pos.col].isFood);

	for (uint idx = 0; idx < m_foods.size(); ++idx) {
		if (m_foods[idx].where == pos)
			return idx;
	}

	state.bug << "nope " << pos << endl;

	abort();
}


void Bot::assign_food()
{
	AStar<LocationEvalZero, StepEvalDefault> astar(state, LocationEvalZero(), StepEvalDefault(state));
	uint unassigned = min(m_foods.size(), m_ants.size());

	while (unassigned > 0) {
		state.bug << unassigned << " unassigned" << endl;
		for (uint foodidx = 0; foodidx < m_foods.size(); ++foodidx) {
			if (!m_foods[foodidx].claimed) {
				astar.push(m_foods[foodidx].where);
			}
		}

		Location where;
		int32_t cost;
		while (astar.step(where, cost)) {
			if (cost > 3 * state.viewradius)
				return;

			if (state.grid[where.row][where.col].ant != 0)
				continue;

			uint antidx = myantidx_at(where);
			Ant & ant = m_ants[antidx];

			if (ant.direction >= 0)
				continue;

			Location foodwhere = astar.getsource(where);
			ant.direction = reversedir(astar.getlaststep(where));

			state.bug << "ant " << antidx << " at " << where
				<< " assigned to food at " << foodwhere << " (" << CDIRECTIONS[ant.direction] << "), distance " << cost << endl;

			m_foods[foodidx_at(foodwhere)].claimed = true;
			unassigned--;
			break;
		}

		astar.clear();
	}
}

bool Bot::try_rotate_move(uint antidx, const Map<bool> & claims)
{
	static int rotate = 1;
	if (rotate == 1)
		rotate = TDIRECTIONS - 1;
	else
		rotate = 1;

	Ant & ant = m_ants[antidx];
	int altdir = (ant.direction + rotate) % TDIRECTIONS;
	Location moveto = state.getLocation(ant.where, altdir);

	state.bug << "  try altdir " << cdir(altdir) << " to " << moveto << endl;

	if (!claims[moveto] && !state.grid[moveto.row][moveto.col].isWater) {
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
	Map<bool> claimed(state.rows, state.cols);

	for (uint foodidx = 0; foodidx < m_foods.size(); ++foodidx)
		claimed[m_foods[foodidx].where] = true;

	for (uint antidx = 0; antidx < m_ants.size(); ++antidx) {
		Ant & ant = m_ants[antidx];
		Location moveto;

		moveto = ant.where;
		if (ant.direction >= 0)
			moveto = state.getLocation(moveto, ant.direction);

		state.bug << "ant at " << ant.where << " initial attempt " << moveto << " dir " << cdir(ant.direction) << endl;

		if (claimed[moveto] && (ant.direction < 0 || !try_rotate_move(antidx, claimed))) {
			if (!claimed[ant.where]) {
				state.bug << "  current pos unclaimed" << endl;
				ant.direction = -1;
			} else {
				state.bug << "  rotate" << endl;
				static int rotate = 0;
				rotate = (rotate + 1) % TDIRECTIONS;
				for (int predir = 0; predir < TDIRECTIONS; ++predir) {
					int dir = (predir + rotate) % TDIRECTIONS;
					moveto = state.getLocation(ant.where, dir);
					if (!claimed[moveto]) {
						ant.direction = dir;
						break;
					}
				}
			}
		}

		state.bug << "  move " << cdir(ant.direction) << endl;

		if (ant.direction >= 0) {
			moveto = state.getLocation(ant.where, ant.direction);
			claimed[moveto] = true;
			state.makeMove(ant.where, ant.direction);
		} else {
			claimed[ant.where] = true;
		}
	}
}


//makes the bots moves for the turn
void Bot::makeMoves()
{
	state.bug << "turn " << state.turn << ":" << endl;
	state.bug << state << endl;

	m_zoc.update();

	//
	m_ants.clear();
	m_ants.reserve(state.myAnts.size());
	for (uint antidx = 0; antidx < state.myAnts.size(); ++antidx) {
		m_ants.push_back(Ant());
		Ant & ant = m_ants.back();
		ant.where = state.myAnts[antidx];
	}

	//
	m_foods.clear();
	m_foods.reserve(state.food.size());
	for (uint foodidx = 0; foodidx < state.food.size(); ++foodidx) {
		m_foods.push_back(Food());
		Food & food = m_foods.back();
		food.where = state.food[foodidx];
	}

	//
	assign_food();

	//
	m_hilldefense.run();

	//
	m_scout.run();

	//
	for (uint antidx = 0; antidx < state.myAnts.size(); ++antidx) {
		Ant & ant = m_ants[antidx];

		state.bug << "ant " << antidx << " at " << ant.where << " dir " << cdir(ant.direction) << endl;

		if (ant.direction < 0) {
			static int rotate = 0;
			rotate = (rotate + 1) % TDIRECTIONS;

			// not looking for food, go towards enemy territory
			uint my = m_zoc.m_enemy[ant.where];
			const int * dirperm = getdirperm();
			for (int predir = 0; predir < TDIRECTIONS; predir++) {
				int dir = dirperm[predir];
				Location loc = state.getLocation(ant.where, dir);

				if (!state.grid[loc.row][loc.col].isWater) {
					if (m_zoc.m_enemy[loc] < my) {
						state.bug << "  zoc move " << cdir(dir) << endl;
						ant.direction = dir;
						break;
					}
				}
			}
		}
	}

	//
	vector<uint> tactical_close;
	vector<uint> tactical_mid;
	vector<uint> tactical_far;

	for (uint antidx = 0; antidx < state.myAnts.size(); ++antidx) {
		Ant & ant = m_ants[antidx];
		uint enemydist = m_zoc.m_enemy[ant.where];
		if (enemydist <= TacticalFar) {
			if (enemydist <= TacticalClose)
				tactical_close.push_back(antidx);
			else if (enemydist <= TacticalMid)
				tactical_mid.push_back(antidx);
			else
				tactical_far.push_back(antidx);
		}
	}

	while (!tactical_close.empty()) {
		Ant & ant = m_ants[tactical_close.back()];
		tactical_close.pop_back();
		if (!ant.hastactical) {
			Tactical tactical(*this);
			tactical.make_moves(ant.where);
			ant.hastactical = true;
		}
	}

	while (!tactical_mid.empty()) {
		Ant & ant = m_ants[tactical_mid.back()];
		tactical_mid.pop_back();
		if (!ant.hastactical) {
			Tactical tactical(*this);
			tactical.make_moves(ant.where);
			ant.hastactical = true;
		}
	}

	while (!tactical_far.empty()) {
		Ant & ant = m_ants[tactical_far.back()];
		tactical_far.pop_back();
		if (!ant.hastactical) {
			Tactical tactical(*this);
			tactical.make_moves(ant.where);
			ant.hastactical = true;
		}
	}

	//
	make_moves();

	state.bug.time << "time taken in turn " << state.turn << ": " << state.timer.getTime() << "ms" << endl << endl;
}

//finishes the turn
void Bot::endTurn()
{
	if(state.turn > 0)
		state.reset();
	state.turn++;

	cout << "go" << endl;
}
