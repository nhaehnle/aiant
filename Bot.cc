#include <cassert>

#include "Bot.h"
#include "astar.h"

using namespace std;

//constructor
Bot::Bot()
{

}

//plays a single game of Ants.
void Bot::playGame()
{
	//reads the game parameters and sets up
	cin >> state;
	state.setup();
	endTurn();

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
				state.bug << "push " << m_foods[foodidx].where << endl;
			}
		}

		Location where;
		int32_t cost;
		while (astar.step(where, cost)) {
			state.bug << "step " << where << " " << cost << endl;

			if (cost > 3 * state.viewradius)
				return;

			if (state.grid[where.row][where.col].ant != 0)
				continue;

			uint antidx = myantidx_at(where);
			Ant & ant = m_ants[antidx];

			if (ant.hasfood)
				continue;

			ant.hasfood = true;
			ant.food.where = astar.getsource(where);
			ant.food.direction = reversedir(astar.getlaststep(where));
			ant.food.distance = cost;

			state.bug << "ant " << antidx << " at " << where
				<< " assigned to food at " << ant.food.where << " (" << CDIRECTIONS[ant.food.direction] << "), distance " << cost << endl;

			m_foods[foodidx_at(ant.food.where)].claimed = true;
			unassigned--;
			break;
		}

		astar.clear();
	}
}

void Bot::make_moves()
{
	for (uint antidx = 0; antidx < m_ants.size(); ++antidx) {
		Ant & ant = m_ants[antidx];
		if (ant.direction >= 0)
			state.makeMove(ant.where, ant.direction);
	}
}


//makes the bots moves for the turn
void Bot::makeMoves()
{
	state.bug << "turn " << state.turn << ":" << endl;
	state.bug << state << endl;

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
	for (uint antidx = 0; antidx < state.myAnts.size(); ++antidx) {
		Ant & ant = m_ants[antidx];
		ant.goal = ant.where;
		ant.direction = -1;

		state.bug << "ant " << antidx << " at " << ant.where << endl;

		if (ant.hasfood) {
			if (ant.food.distance > 1) {
				ant.goal = ant.food.where;
				ant.direction = ant.food.direction;
			}
		} else {
			for (int d = 0; d < TDIRECTIONS; d++) {
				Location loc = state.getLocation(ant.where, d);

				if (!state.grid[loc.row][loc.col].isWater) {
					ant.goal = loc;
					ant.direction = d;
					break;
				}
			}
		}

		state.bug << "  goal " << ant.goal << " dir " << ((ant.direction >= 0) ? CDIRECTIONS[ant.direction] : '-') << endl;
	}

	//
	make_moves();

	state.bug << "time taken: " << state.timer.getTime() << "ms" << endl << endl;
}

//finishes the turn
void Bot::endTurn()
{
	if(state.turn > 0)
		state.reset();
	state.turn++;

	cout << "go" << endl;
}
