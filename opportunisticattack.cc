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
