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
