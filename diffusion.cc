#include "diffusion.h"

#include "Bot.h"
#include "map.h"
#include "zoc.h"

using namespace std;

/**
 * Tweakable parameters
 */
//@{
static float EnemyAntWeight = 3.0;
static float EnemyHillWeight = 10.0;
static float BorderWeight = 0.2;
static float DiffusionDecay = 0.95;
static float DiffusionFlow = 0.4;
//@}

struct Diffusion::Data {
	Map<float> map;
	Map<float> buddy;
};

Diffusion::Diffusion(Bot & b) :
	bot(b),
	state(b.state),
	d(*new Data)
{

}

Diffusion::~Diffusion()
{
	delete &d;
}

void Diffusion::init()
{
	EnemyAntWeight = bot.getargfloat("EnemyAntWeight", EnemyAntWeight);
	EnemyHillWeight = bot.getargfloat("EnemyHillWeight", EnemyHillWeight);
	BorderWeight = bot.getargfloat("BorderWeight", BorderWeight);
	DiffusionDecay = bot.getargfloat("DiffusionDecay", DiffusionDecay);
	DiffusionFlow = bot.getargfloat("DiffusionFlow", DiffusionFlow);

	state.bug.time << "Diffusion parameters: EnemyAntWeight = " << EnemyAntWeight
		<< ", EnemyHillWeight = " << EnemyHillWeight
		<< ", BorderWeight = " << BorderWeight
		<< ", DiffusionDecay = " << DiffusionDecay
		<< ", DiffusionFlow = " << DiffusionFlow << endl;

	d.map.resize(state.rows, state.cols);
	d.buddy.resize(state.rows, state.cols);
}

void Diffusion::diffuse()
{
	Location cur;
	for (cur.row = 0; cur.row < state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < state.cols; ++cur.col) {
			float value = d.map[cur] * DiffusionDecay;
			d.map[cur] = 0.0;

			d.buddy[cur] += value * (1.0 - DiffusionFlow);
			for (int dir = 0; dir < TDIRECTIONS; ++dir) {
				Location n = state.getLocation(cur, dir);
				if (state.grid[n.row][n.col].isWater)
					continue;

				d.buddy[n] += value * DiffusionFlow * 0.25;
			}
		}
	}

	for (cur.row = 0; cur.row < state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < state.cols; ++cur.col) {
			float value = d.buddy[cur] * DiffusionDecay;
			d.buddy[cur] = 0.0;

			float flow = value * DiffusionFlow * 0.25;
			for (int dir = 0; dir < TDIRECTIONS; ++dir) {
				Location n = state.getLocation(cur, dir);
				if (state.grid[n.row][n.col].isWater)
					continue;

				d.map[n] += flow;
				value -= flow;
			}
			d.map[cur] += value;
		}
	}
}

void Diffusion::run()
{
	//
	for (uint antidx = 0; antidx < bot.m_ants.size(); ++antidx) {
		Ant & ant = bot.m_ants[antidx];

		d.map[ant.where] += 1.0;
	}

	for (uint enemyidx = 0; enemyidx < state.enemyAnts.size(); ++enemyidx) {
		d.map[state.enemyAnts[enemyidx]] -= EnemyAntWeight;
	}

	for (uint enemyhill = 0; enemyhill < state.enemyHills.size(); ++enemyhill) {
		d.map[state.enemyHills[enemyhill]] -= EnemyHillWeight;
	}

	Location cur;
	for (cur.row = 0; cur.row < state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < state.cols; ++cur.col) {
			int delta = (int)bot.m_zoc.m_enemy[cur] - (int)bot.m_zoc.m_me[cur];

			if (1 >= delta && delta >= -1)
				d.map[cur] -= BorderWeight;

			if (isnan(d.map[cur])) {
				state.bug.time << "Diffusion NaN at " << cur << endl;
				d.map[cur] = 0.0;
			}
		}
	}

	//
	for (int iters = 0; iters < 2; ++iters)
		diffuse();

	//
	for (uint antidx = 0; antidx < bot.m_ants.size(); ++antidx) {
		Ant & ant = bot.m_ants[antidx];

		if (ant.assigneddirection)
			continue;

		float bestvalue = 0.0;
		int bestdir = -1;
		for (int predir = 0; predir < TDIRECTIONS; ++predir) {
			int dir = ant.dirperm[predir];
			Location n = state.getLocation(ant.where, dir);

			if (state.grid[n.row][n.col].isWater)
				continue;

			float value = d.map[n];
			if (value < bestvalue) {
				bestvalue = value;
				bestdir = dir;
			}
		}

		if (bestdir >= 0 && bestvalue < d.map[ant.where]) {
			state.bug << "diffusion move ant at " << ant.where << " to " << cdir(bestdir) << endl;
			ant.direction = bestdir;
			ant.assigneddirection = true;
		}
	}
}
