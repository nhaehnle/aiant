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

TacticalSmBase::TacticalSmBase(Bot & b) :
	bot(b),
	state(b.state)
{
	ValueKill = 3.0;
	ValueLoss = 0.25;
	ValueHill = 2.0;
	ValueEnemyDist = 0.999;

	MinAggressionAnts = 4;
	AggressionThresholdShift = -6.0;
	AggressionScale = 1.0;
}

void TacticalSmBase::init()
{
	ValueKill = bot.getargfloat("ValueKill", ValueKill);
	ValueLoss = bot.getargfloat("ValueLoss", ValueLoss);
	ValueHill = bot.getargfloat("ValueHill", ValueHill);
	ValueEnemyDist = bot.getargfloat("ValueEnemyDist", ValueEnemyDist);

	AggressionThresholdShift = bot.getargfloat("AggressionThresholdShift", AggressionThresholdShift);
	AggressionScale = bot.getargfloat("AggressionScale", AggressionScale);

	state.bug.time << "Tactical parameters: ValueKill = " << ValueKill
		<< ", ValueLoss = " << ValueLoss
		<< ", ValueHill = " << ValueHill
		<< ", ValueEnemyDist = " << ValueEnemyDist << endl
		<< "  AggressionThresholdShift = " << AggressionThresholdShift
		<< ", AggressionScale = " << AggressionScale << endl;
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

