#include "tactical_sample.h"

#include <cassert>

#include "Bot.h"
#include "map.h"
#include "zoc.h"

using namespace std;

/**
 * Tunable parameters
 */
//@{
static const uint MaxFoodDist = 16;
static const uint MaxHillDist = 3;
static float ValueAnts[2] = { 6.0, 5.0 };
static float ValueHill = 2.0;
static float ValueEat = 1.0;
static float ValueEnemyDist = 0.0001;
static float ValueFoodDist = 0.001;
static float ValueHillDist = 0.5;
static float ValueInitialMove = 0.05;

static const uint MinAggressionAnts = 2;
static float AggressionThresholdShift = -10.0;
static float AggressionScale = 1.0;

static float LearnGamma = 1.1;
//@}

static const float ValueEpsilon = 0.00000001;

static const int AggressionKernelRoughRadius = 7;
static const uint AggressionKernelRadius2 = 61;

static const int AttackOffsets[20][2] = {
	{ -2, -1 },
	{ -2, 0 },
	{ -2, 1 },
	{ -1, -2 },
	{ -1, -1 },
	{ -1, 0 },
	{ -1, 1 },
	{ -1, 2 },
	{ 0, -2 },
	{ 0, -1 },
	{ 0, 1 },
	{ 0, 2 },
	{ 1, -2 },
	{ 1, -1 },
	{ 1, 0 },
	{ 1, 1 },
	{ 1, 2 },
	{ 2, -1 },
	{ 2, 0 },
	{ 2, 1 }
};
static const uint NrAttackOffsets = 20;

static const uint MINE = 0;
static const uint ENEMY = 1;

struct TacticalSample::Data {
	struct Ant {
		struct Attacker {
			uint32_t antidx;
			uint32_t weakness;
		};

		Location origpos;
		uint8_t owner : 1;
		uint8_t tactical : 1;
		uint8_t aggressive : 1;
		uint8_t haveinitialdirection : 1;
		int8_t initialdirection;

		Location pos;
		uint8_t move;
		bool killed;
		Attacker attackers[20];
		uint nrattackers;

		Ant(const Location & origpos_, uint8_t owner_) :
			origpos(origpos_),
			owner(owner_),
			tactical(0),
			aggressive(0), haveinitialdirection(0), initialdirection(-1),
			move(0), killed(false),
			nrattackers(0)
		{
		}
	};
	struct Field {
		uint32_t fooddist;
		int32_t antidx;
		uint8_t influence[2];
		uint8_t hilldist[2];
	};
	struct Sampler {
		uint32_t antidx;
		float w[5];
		uint cool : 1;
		uint unchangedrounds : 3;

		void reset() {
			cool = 0;
			unchangedrounds = 0;
		}
	};

	Bot & bot;
	uint nrsampleruns;
	vector<Ant> ants;
	Map<Field> map;
	vector<Sampler *> samplers;
	vector<Sampler *> coolsamplers;

	float value;

	vector<Sampler *> unused_samplers;

	vector<Location> locqueue;

	Data(Bot & b);
	~Data();

	void reset();
	Sampler * allocsampler();

	void mark_hills(uint who, const vector<Location> & hills);

	void ant_killed(uint32_t antidx, bool killed);

	void add_attacker(uint32_t antidx, uint32_t attackeridx);
	void weakness_increase(uint32_t antidx);
	void remove_attacker(uint32_t antidx, uint32_t attackeridx);
	void weakness_decrease(uint32_t antidx);

	void mark_ant(uint antidx, const Location & pos);
	void unmark_ant(uint antidx);

	Location add_location(const Location & p, int row, int col) {
		return Location
			((p.row + row + map.rows) % map.rows,
			 (p.col + col + map.cols) % map.cols);
	}

	float position_value(const Location & p, uint who);
};

TacticalSample::Data::Data(Bot & b) :
	bot(b)
{
}

TacticalSample::Data::~Data()
{
	reset();

	while (!unused_samplers.empty()) {
		delete unused_samplers.back();
		unused_samplers.pop_back();
	}
}

void TacticalSample::Data::reset()
{
	nrsampleruns = 0;

	unused_samplers.insert(unused_samplers.end(), samplers.begin(), samplers.end());
	samplers.clear();

	unused_samplers.insert(unused_samplers.end(), coolsamplers.begin(), coolsamplers.end());
	coolsamplers.clear();
}

TacticalSample::Data::Sampler * TacticalSample::Data::allocsampler()
{
	Sampler * s;
	if (unused_samplers.empty()) {
		s = new Sampler;
	} else {
		s = unused_samplers.back();
		unused_samplers.pop_back();
	}
	s->reset();
	return s;
}

float TacticalSample::Data::position_value(const Location & p, uint who)
{
	float v = ValueFoodDist * (MaxFoodDist - map[p].fooddist);

	if (map[p].fooddist == 1)
		v += ValueEat;

	if (who == MINE)
		v += -ValueEnemyDist * bot.m_zoc.m_enemy[p];
	else
		v += -ValueEnemyDist * bot.m_zoc.m_me[p];

	v += ValueHillDist * (MaxHillDist - map[p].hilldist[1 - who]);
	if (map[p].hilldist[1 - who] == 0)
		v += ValueHill;

	return v;
}

void TacticalSample::Data::ant_killed(uint32_t antidx, bool killed)
{
	Ant & ant = ants[antidx];

	assert(ant.killed != killed);
	ant.killed = killed;

	float v = ValueAnts[ant.owner];
	if (ant.aggressive)
		v = ValueAnts[ENEMY];

	if ((ant.owner == MINE) == killed)
		value -= v;
	else
		value += v;
}

void TacticalSample::Data::add_attacker(uint32_t antidx, uint32_t attackeridx)
{
	Ant & ant = ants[antidx];
	Ant & attacker = ants[attackeridx];

	assert(ant.owner != attacker.owner);
	assert(ant.nrattackers < 20);

	uint attackerweakness = map[attacker.pos].influence[ant.owner];

	ant.attackers[ant.nrattackers].antidx = attackeridx;
	ant.attackers[ant.nrattackers].weakness = attackerweakness;
	ant.nrattackers++;

	if (!ant.killed && attackerweakness <= map[ant.pos].influence[attacker.owner])
		ant_killed(antidx, true);
}

void TacticalSample::Data::remove_attacker(uint32_t antidx, uint32_t attackeridx)
{
	Ant & ant = ants[antidx];
	uint lowestweakness = numeric_limits<uint>::max();

	for (uint i = 0; i < ant.nrattackers; ++i) {
		if (ant.attackers[i].antidx == attackeridx) {
			ant.attackers[i] = ant.attackers[ant.nrattackers - 1];
			ant.nrattackers--;
			i--;
		} else {
			lowestweakness = min(lowestweakness, ant.attackers[i].weakness);
		}
	}

	if (ant.killed && lowestweakness > map[ant.pos].influence[1 - ant.owner])
		ant_killed(antidx, false);
}

void TacticalSample::Data::weakness_increase(uint32_t antidx)
{
	Ant & ant = ants[antidx];
	uint newweakness = map[ant.pos].influence[1 - ant.owner];
	(void)newweakness;

	for (uint i = 0; i < ant.nrattackers; ++i) {
		uint otheridx = ant.attackers[i].antidx;
		Ant & other = ants[otheridx];

		uint lowestweakness = numeric_limits<uint>::max();
		for (uint j = 0; j < other.nrattackers; ++j) {
			if (other.attackers[j].antidx == antidx) {
				other.attackers[j].weakness++;
				assert(other.attackers[j].weakness == newweakness);
			}
			lowestweakness = min(lowestweakness, other.attackers[j].weakness);
		}

		if (other.killed && lowestweakness > map[other.pos].influence[ant.owner])
			ant_killed(otheridx, false);
	}
}

void TacticalSample::Data::weakness_decrease(uint32_t antidx)
{
	Ant & ant = ants[antidx];
	uint newweakness = map[ant.pos].influence[1 - ant.owner];
	(void)newweakness;

	for (uint i = 0; i < ant.nrattackers; ++i) {
		uint otheridx = ant.attackers[i].antidx;
		Ant & other = ants[otheridx];

		uint lowestweakness = numeric_limits<uint>::max();
		for (uint j = 0; j < other.nrattackers; ++j) {
			if (other.attackers[j].antidx == antidx) {
				other.attackers[j].weakness--;
				assert(other.attackers[j].weakness == newweakness);
			}
			lowestweakness = min(lowestweakness, other.attackers[j].weakness);
		}

		if (!other.killed && lowestweakness <= map[other.pos].influence[ant.owner])
			ant_killed(otheridx, true);
	}
}

void TacticalSample::Data::mark_ant(uint antidx, const Location & pos)
{
	assert(map[pos].antidx < 0);

	Ant & ant = ants[antidx];

	assert(!ant.killed);

	ant.pos = pos;
	map[pos].antidx = antidx;

	for (uint i = 0; i < NrAttackOffsets; ++i) {
		Location n = add_location(pos, AttackOffsets[i][0], AttackOffsets[i][1]);
		map[n].influence[ant.owner]++;
		int otheridx = map[n].antidx;
		if (otheridx >= 0) {
			Ant & other = ants[otheridx];
			if (other.owner != ant.owner) {
				weakness_increase(otheridx);
				add_attacker(otheridx, antidx);
				add_attacker(antidx, otheridx);
			}
		}
	}

	// death effects have been added implicitly, now take into account other parts of value
	float v = position_value(pos, ant.owner);
	if (ant.owner == MINE)
		value += v;
	else
		value -= v;
}

void TacticalSample::Data::unmark_ant(uint antidx)
{
	Ant & ant = ants[antidx];

	assert(map[ant.pos].antidx == (int32_t)antidx);

	map[ant.pos].antidx = -1;

	for (uint i = 0; i < NrAttackOffsets; ++i) {
		Location n = add_location(ant.pos, AttackOffsets[i][0], AttackOffsets[i][1]);
		map[n].influence[ant.owner]--;
		int otheridx = map[n].antidx;
		if (otheridx >= 0) {
			Ant & other = ants[otheridx];
			if (other.owner != ant.owner) {
				remove_attacker(otheridx, antidx);
				remove_attacker(antidx, otheridx);
				weakness_decrease(otheridx);
			}
		}
	}

	assert(!ant.killed);

	float v = position_value(ant.pos, ant.owner);
	if (ant.owner == MINE)
		value -= v;
	else
		value += v;
}


TacticalSample::TacticalSample(Bot& b) :
	bot(b),
	state(b.state),
	d(*new Data(b))
{

}

TacticalSample::~TacticalSample()
{
	delete &d;
}

void TacticalSample::init()
{
	d.map.resize(state.rows, state.cols);

	ValueAnts[0] = bot.getargfloat("ValueAnts0", ValueAnts[0]);
	ValueAnts[1] = bot.getargfloat("ValueAnts1", ValueAnts[1]);
	ValueHill = bot.getargfloat("ValueHill", ValueHill);
	ValueEat = bot.getargfloat("ValueEat", ValueHill);
	ValueEnemyDist = bot.getargfloat("ValueEnemyDist", ValueEnemyDist);
	ValueFoodDist = bot.getargfloat("ValueFoodDist", ValueFoodDist);

	AggressionThresholdShift = bot.getargfloat("AggressionThresholdShift", AggressionThresholdShift);
	AggressionScale = bot.getargfloat("AggressionScale", AggressionScale);

	LearnGamma = bot.getargfloat("LearnGamma", LearnGamma);

	state.bug.time << "Tactical parameters: ValueAnts[0] = " << ValueAnts[0]
		<< ", ValueAnts[1] = " << ValueAnts[1]
		<< ", ValueHill = " << ValueHill
		<< ", ValueEat = " << ValueEat << endl
		<< "  ValueEnemyDist = " << ValueEnemyDist
		<< ", ValueFoodDist = " << ValueFoodDist << endl
		<< "  AggressionThresholdShift = " << AggressionThresholdShift
		<< ", AggressionScale = " << AggressionScale << endl
		<< "  LearnGamma = " << LearnGamma << endl;
}

void TacticalSample::do_sample_markonly(uint sampleridx)
{
	Data::Sampler & s = *d.samplers[sampleridx];
	Data::Ant & ant = d.ants[s.antidx];

	bool blocked[5];
	float totalweight = 0.0;
	for (int move = 0; move < 5; ++move) {
		if (s.w[move] <= 0.0) {
			blocked[move] = true;
			continue;
		}

		Location n;
		if (move == 0)
			n = ant.origpos;
		else
			n = state.getLocation(ant.origpos, move-1);

		if (d.map[n].antidx >= 0) {
			blocked[move] = true;
			continue;
		}

		blocked[move] = false;
		totalweight += s.w[move];
	}

	float r = fastrngd() * totalweight;
	int domove = 0;
	for (int move = 0; move < 5 && r >= 0.0; ++move) {
		if (!blocked[move]) {
			domove = move;
			r -= s.w[move];
		}
	}

	Location n;
	if (domove == 0)
		n = ant.origpos;
	else
		n = state.getLocation(ant.origpos, domove-1);

	ant.move = domove;
	d.mark_ant(s.antidx, n);
}

void TacticalSample::run_samplers()
{
	uint p = getprime();
	uint ofs = fastrng() % d.samplers.size();
	for (uint preidx = 0; preidx < d.samplers.size(); ++preidx) {
		uint idx = (p * preidx + ofs) % d.samplers.size();
		Data::Sampler & s = *d.samplers[idx];
		Data::Ant & ant = d.ants[s.antidx];

		state.bug << "  s " << idx << ", ant @ " << ant.origpos;
		d.nrsampleruns++;

		float bestvalue = 0.0;
		float values[5];
		int origmove = ant.move;

		d.value = 0.0;
		values[origmove] = 0.0;
		d.unmark_ant(s.antidx);

		for (int move = (origmove + 1) % 5; move != origmove; move = (move + 1) % 5) {
			if (s.w[move] < 0.0) {
				values[move] = numeric_limits<float>::min();
				continue;
			}

			Location n;
			if (move == 0)
				n = ant.origpos;
			else
				n = state.getLocation(ant.origpos, move - 1);

			if (d.map[n].antidx >= 0) {
				values[move] = numeric_limits<float>::min();
				continue;
			}

			d.mark_ant(s.antidx, n);
			if (ant.owner == MINE) {
				values[move] = d.value;
				if (ant.haveinitialdirection && (move - 1) == ant.initialdirection)
					values[move] += ValueInitialMove;
			} else {
				values[move] = -d.value;
			}
			d.unmark_ant(s.antidx);

			bestvalue = max(bestvalue, values[move]);
		}

		state.bug << ", new w =";

		for (int move = 0; move < 5; ++move) {
			if (values[move] >= bestvalue - ValueEpsilon)
				s.w[move] *= LearnGamma;
			state.bug << " " << s.w[move];
		}

		state.bug << endl;

		do_sample_markonly(idx);

		if (ant.move == origmove) {
			s.unchangedrounds++;
			if (s.unchangedrounds >= 4)
				s.cool = 1;
		} else {
			s.unchangedrounds = 0;
		}
	}
}

void TacticalSample::mark_hills(uint who, const vector<Location> & hills)
{
	d.locqueue.clear();
	for (uint idx = 0; idx < hills.size(); ++idx) {
		const Location & hill = hills[idx];
		d.map[hill].hilldist[who] = 0;
		d.locqueue.push_back(hill);
	}

	uint queue_head = 0;
	while (queue_head < d.locqueue.size()) {
		Location cur = d.locqueue[queue_head++];
		uint ndist = d.map[cur].hilldist[who] + 1;

		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n = state.getLocation(cur, dir);
			if (state.grid[n.row][n.col].isWater)
				continue;

			if (d.map[n].hilldist[who] > ndist) {
				d.map[n].hilldist[who] = ndist;
				d.locqueue.push_back(n);
			}
		}
	}
}

void TacticalSample::init_map()
{
	Data::Field initfield;
	initfield.fooddist = MaxFoodDist;
	initfield.influence[0] = initfield.influence[1] = 0;
	initfield.hilldist[0] = initfield.hilldist[1] = MaxHillDist;
	initfield.antidx = -1;
	d.map.fill(initfield);

	//
	d.locqueue.clear();
	for (uint idx = 0; idx < state.food.size(); ++idx) {
		const Location & food = state.food[idx];
		d.map[food].fooddist = 0;
		d.locqueue.push_back(food);
	}

	uint queue_head = 0;
	while (queue_head < d.locqueue.size()) {
		Location cur = d.locqueue[queue_head++];
		uint ndist = d.map[cur].fooddist + 1;

		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n = state.getLocation(cur, dir);
			if (state.grid[n.row][n.col].isWater)
				continue;
			if (ndist < d.map[n].fooddist) {
				d.map[n].fooddist = ndist;
				d.locqueue.push_back(n);
			}
		}
	}

	//
	mark_hills(MINE, state.myHills);
	mark_hills(ENEMY, state.enemyHills);
}

void TacticalSample::init_ants()
{
	d.ants.clear();

	for (uint idx = 0; idx < bot.m_ants.size(); ++idx) {
		Ant & ant = bot.m_ants[idx];
		d.ants.push_back(Data::Ant(ant.where, MINE));
		if (ant.assigneddirection) {
			d.ants.back().haveinitialdirection = 1;
			d.ants.back().initialdirection = ant.direction;
		}

		d.mark_ant(d.ants.size() - 1, ant.where);
	}

	for (uint idx = 0; idx < state.enemyAnts.size(); ++idx) {
		const Location & enemy = state.enemyAnts[idx];

		d.ants.push_back(Data::Ant(enemy, ENEMY));

		d.mark_ant(d.ants.size() - 1, enemy);
	}

	//
	for (uint idx = 0; idx < d.ants.size(); ++idx) {
		Data::Ant & ant = d.ants[idx];
		uint ants[2] = { 0, 0 };

		Location ofs;
		for (ofs.row = -AggressionKernelRoughRadius; ofs.row <= AggressionKernelRoughRadius; ++ofs.row) {
			for (ofs.col = -AggressionKernelRoughRadius; ofs.col <= AggressionKernelRoughRadius; ++ofs.col) {
				Location n = state.addLocations(ant.origpos, ofs);
				int32_t otheridx = d.map[n].antidx;

				if
					(otheridx >= 0 && (uint)otheridx != idx &&
					 state.eucliddist2(ant.origpos, n) <= AggressionKernelRadius2)
				{
					Data::Ant & other = d.ants[d.map[n].antidx];
					ants[other.owner != ant.owner]++;
				}
			}
		}

		if (ants[1] == 0)
			continue;

		ant.tactical = 1;
		if (ants[0] >= MinAggressionAnts) {
			double antratio = (double)ants[0] / (double)ants[1];
			double p = 1.0 / (1.0 + exp(-AggressionScale * (log(antratio) + AggressionThresholdShift)));

			if (fastrngd() <= p)
				ant.aggressive = 1;
		}

		Data::Sampler & s = *d.allocsampler();
		s.antidx = idx;
		d.samplers.push_back(&s);

		//
		for (int move = 0; move < 5; ++move) {
			Location n;
			if (move == 0)
				n = ant.origpos;
			else
				n = state.getLocation(ant.origpos, move-1);

			Square & sq = state.grid[n.row][n.col];
			if (sq.isWater || sq.isFood)
				s.w[move] = -1.0;
			else
				s.w[move] = 1.0;
		}
	}
}

void TacticalSample::run()
{
	d.reset();

	init_map();
	init_ants();

	if (d.samplers.empty())
		return;

	state.bug << "Tactical turn " << state.turn << " on " << d.samplers.size() << " samplers" << endl;

	//
	uint reheated = 0;
	uint rounds = 0;
	while (!bot.timeover()) {
		rounds++;
		run_samplers();

		for (uint idx = 0; idx < d.samplers.size(); ++idx) {
			if (d.samplers[idx]->cool) {
				d.coolsamplers.push_back(d.samplers[idx]);
				d.samplers[idx] = d.samplers.back();
				d.samplers.pop_back();
				idx--;
			}
		}

		for (uint idx = 0; idx < d.coolsamplers.size(); ++idx) {
			if (!d.coolsamplers[idx]->cool) {
				d.samplers.push_back(d.coolsamplers[idx]);
				d.coolsamplers[idx] = d.coolsamplers.back();
				d.coolsamplers.pop_back();
				idx--;
			}
		}

		if (d.samplers.empty()) {
			if (reheated >= 3)
				break;

			state.bug << "Reheating after round " << rounds << endl;
			swap(d.samplers, d.coolsamplers);
			for (uint idx = 0; idx < d.samplers.size(); ++idx) {
				Data::Sampler & s = *d.samplers[idx];
				s.cool = 0;
				s.unchangedrounds = 0;
			}
			reheated++;
		}
	}

	//
	for (uint idx = 0; idx < bot.m_ants.size(); ++idx) {
		Data::Ant & dant = d.ants[idx];
		Ant & bant = bot.m_ants[idx];

		assert(dant.origpos == bant.where);

		if (dant.tactical) {
			bant.direction = dant.move - 1;
			state.bug << "Tactical decision: ant at " << bant.where << " move to " << cdir(bant.direction) << endl;
		}
	}

	state.bug.time << "Total sample runs: " << d.nrsampleruns
		<< ", samplers: " << (d.samplers.size() + d.coolsamplers.size())
		<< ", rounds: " << rounds << endl;
}
