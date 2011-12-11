#include "opportunisticattack.h"

#include "Bot.h"
#include "astar.h"

using namespace std;

static const int HillAttackMaxSteps = 10;

OpportunisticAttack::OpportunisticAttack(Bot & b) :
	bot(b),
	state(b.state)
{
}

OpportunisticAttack::~OpportunisticAttack()
{
}

void OpportunisticAttack::init()
{
}

void OpportunisticAttack::run()
{
	for (uint hillidx = 0; hillidx < state.enemyHills.size(); ++hillidx) {
		const Location & enemypos = state.enemyHills[hillidx];

		state.bug << "check opportunistic attack of hill at " << enemypos << endl;

		AStar<LocationEvalZero, StepEvalDefault> astar
			(state, LocationEvalZero(), StepEvalDefault(state));

		astar.push(enemypos);

		Location cur;
		int32_t cost;
		while (astar.step(cur, cost)) {
			if (cost > HillAttackMaxSteps)
				break;

			if (state.grid[cur.row][cur.col].ant == 0) {
				uint antidx = bot.myantidx_at(cur);
				Ant & ant = bot.m_ants[antidx];
				if (!ant.assigneddirection) {
					ant.direction = reversedir(astar.getlaststep(cur));
					ant.assigneddirection = true;
					state.bug << "  send ant " << antidx << " at " << ant.where
						<< " in direction " << cdir(ant.direction) << endl;
				}
			}
		}
	}
}
