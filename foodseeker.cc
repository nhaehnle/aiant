#include "foodseeker.h"

#include <algorithm>

#include "astar.h"
#include "Bot.h"

using namespace std;

struct Food {
	Location where;
	bool claimed;

	Food() : claimed(false) {}
};


struct FoodSeeker::Data {
	vector<Food> foods;
};

FoodSeeker::FoodSeeker(Bot & b) :
	bot(b),
	state(b.state),
	d(*new Data)
{
}

FoodSeeker::~FoodSeeker()
{
	delete &d;
}

void FoodSeeker::init()
{
}

uint FoodSeeker::foodidx_at(const Location & pos)
{
	state.bug << "foodidx " << pos << endl;

	assert(state.grid[pos.row][pos.col].isFood);

	for (uint idx = 0; idx < d.foods.size(); ++idx) {
		if (d.foods[idx].where == pos)
			return idx;
	}

	state.bug << "nope " << pos << endl;

	abort();
}


void FoodSeeker::assign_food()
{
	AStar<LocationEvalZero, StepEvalDefault> astar(state, LocationEvalZero(), StepEvalDefault(state));
	uint unassigned = min(d.foods.size(), bot.m_ants.size());

	while (unassigned > 0) {
		state.bug << unassigned << " unassigned" << endl;
		for (uint foodidx = 0; foodidx < d.foods.size(); ++foodidx) {
			if (!d.foods[foodidx].claimed) {
				astar.push(d.foods[foodidx].where);
			}
		}

		Location where;
		int32_t cost;
		while (astar.step(where, cost)) {
			if (cost > 3 * state.viewradius)
				return;

			if (state.grid[where.row][where.col].ant != 0)
				continue;

			uint antidx = bot.myantidx_at(where);
			Ant & ant = bot.m_ants[antidx];

			if (ant.assigneddirection)
				continue;

			Location foodwhere = astar.getsource(where);
			if (cost <= 1)
				ant.direction = -1;
			else
				ant.direction = reversedir(astar.getlaststep(where));
			ant.assigneddirection = true;

			state.bug << "ant " << antidx << " at " << where
				<< " assigned to food at " << foodwhere << " (" << CDIRECTIONS[ant.direction] << "), distance " << cost << endl;

			d.foods[foodidx_at(foodwhere)].claimed = true;
			unassigned--;
			break;
		}

		astar.clear();
	}
}

void FoodSeeker::run()
{
	d.foods.clear();
	d.foods.reserve(state.food.size());
	for (uint foodidx = 0; foodidx < state.food.size(); ++foodidx) {
		d.foods.push_back(Food());
		Food & food = d.foods.back();
		food.where = state.food[foodidx];
	}

	//
	assign_food();
}
