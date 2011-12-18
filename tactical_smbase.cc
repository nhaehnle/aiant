#include "tactical_smbase.h"

#include "Bot.h"

using namespace std;

const uint8_t Submap::Mine; // must be the two least significant bits
const uint8_t Submap::Enemy;

ostream & operator<<(ostream & out, const Submap & sm)
{
	Location cur;
	for (cur.row = 0; cur.row < Submap::Size; ++cur.row) {
		for (cur.col = 0; cur.col < Submap::Size; ++cur.col) {
			unsigned char field = sm[cur];
			if (field & Submap::Water) {
				out << '%';
			} else if (field & Submap::Food) {
				out << '*';
			} else if (field & Submap::Hill) {
				if (field & Submap::Mine)
					out << 'H';
				else
					out << 'h';
			} else if (field & Submap::Ant) {
				if (field & Submap::Mine)
					out << 'A';
				else
					out << 'a';
			} else {
				out << ' ';
			}
		}
		out << endl;
	}
	return out;
}

struct TacticalSmBase::BaseData {
	vector<Location> queue;
};

TacticalSmBase::TacticalSmBase(Bot & b) :
	bot(b),
	state(b.state),
	bd(*new BaseData)
{
	ValueKill = 3.0;
	ValueLoss = 0.25;
	ValueHill = 2.0;
	ValueEnemyDist = 0.9999;
	ValueFoodDist = 0.99;
	ValueEat = 1.334;

	MinAggressionAnts = 4;
	AggressionThresholdShift = -6.0;
	AggressionScale = 1.0;
}

TacticalSmBase::~TacticalSmBase()
{
	delete &bd;
}

void TacticalSmBase::init()
{
	ValueKill = bot.getargfloat("ValueKill", ValueKill);
	ValueLoss = bot.getargfloat("ValueLoss", ValueLoss);
	ValueHill = bot.getargfloat("ValueHill", ValueHill);
	ValueEat = bot.getargfloat("ValueEat", ValueEat);
	ValueEnemyDist = bot.getargfloat("ValueEnemyDist", ValueEnemyDist);
	ValueFoodDist = bot.getargfloat("ValueFoodDist", ValueFoodDist);

	AggressionThresholdShift = bot.getargfloat("AggressionThresholdShift", AggressionThresholdShift);
	AggressionScale = bot.getargfloat("AggressionScale", AggressionScale);

	state.bug.time << "Tactical parameters: ValueKill = " << ValueKill
		<< ", ValueLoss = " << ValueLoss
		<< ", ValueHill = " << ValueHill << endl
		<< "  ValueEnemyDist = " << ValueEnemyDist
		<< ", ValueFoodDist = " << ValueFoodDist
		<< ", ValueEat = " << ValueEat << endl
		<< "  AggressionThresholdShift = " << AggressionThresholdShift
		<< ", AggressionScale = " << AggressionScale << endl;

	fooddist.resize(state.rows, state.cols);
}

void TacticalSmBase::gensubmap_field(Submap & sm, const Location & local, const Location & global)
{
	Square & sq = state.grid[global.row][global.col];
	unsigned char & field = sm[local];
	if (sq.isFood) {
		field = Submap::Food;
	} else {
		field = 0;
		if (sq.ant >= 0) {
			field |= Submap::Ant;
			if (sq.ant > 0)
				field |= Submap::Enemy;
			else
				field |= Submap::Mine;
		}
		if (sq.isHill) {
			if (sq.hillPlayer > 0)
				field |= Submap::Hill | Submap::Enemy;
			else if (sq.hillPlayer == 0)
				field |= Submap::Hill | Submap::Mine;
		}
	}
}

/**
 * Generate a tactical map centered at \p center
 */
void TacticalSmBase::gensubmap(Submap & sm, const Location & center)
{
	Location local;
	for (local.row = 0; local.row < Submap::Size; ++local.row) {
		for (local.col = 0; local.col < Submap::Size; ++local.col) {
			Location global
				((center.row + local.row + state.rows - Submap::Radius) % state.rows,
				 (center.col + local.col + state.cols - Submap::Radius) % state.cols);

			if (state.grid[global.row][global.col].isWater) {
				sm[local] = Submap::Water;
			} else {
				gensubmap_field(sm, local, global);
			}
		}
	}
}

void TacticalSmBase::compute_fooddists()
{
	fooddist.fill(MaxFoodDist);
	bd.queue.clear();

	for (uint idx = 0; idx < state.food.size(); ++idx) {
		const Location & food = state.food[idx];
		fooddist[food] = 0;
		bd.queue.push_back(food);
	}

	uint queue_head = 0;
	while (queue_head < bd.queue.size()) {
		Location cur = bd.queue[queue_head++];
		uint ndist = fooddist[cur] + 1;

		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n = state.getLocation(cur, dir);
			if (state.grid[n.row][n.col].isWater)
				continue;

			if (fooddist[n] > ndist) {
				fooddist[n] = ndist;
				bd.queue.push_back(n);
			}
		}
	}
}
