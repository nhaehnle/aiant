#include "tactical.h"

#include <cassert>
#include <cstring>
#include <cstdarg>
#include <cstdio>
#include <stdlib.h>

#include "Bot.h"
#include "State.h"
#include <limits>

using namespace std;

/// Tweakable parameters
//@{
static const int ValueField = 5;
static const int ValueFood = 5;

static const int FoodPotentialRadius = 4;
static const int FoodPotential[2 * FoodPotentialRadius + 1][2 * FoodPotentialRadius + 1] = {
	{ 0, 0, 0, 0, 1, 0, 0, 0, 0, },
	{ 0, 0, 0, 1, 2, 1, 0, 0, 0, },
	{ 0, 0, 1, 2, 3, 2, 1, 0, 0, },
	{ 0, 1, 2, 3, 4, 3, 2, 1, 0, },
	{ 1, 2, 3, 4, 5, 4, 3, 2, 1, },
	{ 0, 1, 2, 3, 4, 3, 2, 1, 0, },
	{ 0, 0, 1, 2, 3, 2, 1, 0, 0, },
	{ 0, 0, 0, 1, 2, 1, 0, 0, 0, },
	{ 0, 0, 0, 0, 1, 0, 0, 0, 0, },
};

static const int ValueKill = 1000;
static const int ValueLoss = -1200;
static const int ValueEat = 200;

static const uint MaxEval = 50;
//@}

// area in which enemies are counted as attackers
static const int AttackKernelRadius = 2;
static const int AttackKernelSize = 2 * AttackKernelRadius + 1;
static const bool AttackKernel[AttackKernelSize][AttackKernelSize] = {
	{ false, true,  true,  true,  false },
	{ true,  true,  true,  true,  true  },
	{ true,  true,  false, true,  true  },
	{ true,  true,  true,  true,  true  },
	{ false, true,  true,  true,  false },
};

struct Submap {
	static const int Radius = 8;
	static const int Size = 2 * Radius + 1;

	static const unsigned char AntCollision = 0x80;
	static const unsigned char AntMoved = 0x40;
	static const unsigned char Water = 0x20;
	static const unsigned char Food = 0x10;
	static const unsigned char Hill = 0x08;
	static const unsigned char Ant = 0x04;
	static const unsigned char Mine = 0x02; // must be the two least significant bits
	static const unsigned char Enemy = 0x01;

	unsigned char map[Size * Size];

	void fill(unsigned char c) {memset(&map, c, sizeof(map));}
	unsigned char & operator[](const Location & pos) {return map[pos.row * Size + pos.col];}
	const unsigned char & operator[](const Location & pos) const {return map[pos.row * Size + pos.col];}
	bool getneighbour(Location local, int direction, Location & out) const {
		out.row = local.row + DIRECTIONS[direction][0];
		out.col = local.col + DIRECTIONS[direction][1];
		return out.row >= 0 && out.row < Size && out.col >= 0 && out.col < Size;
	}
	bool getneighbouropt(Location local, int direction, Location & out) const {
		if (direction >= 0) {
			out.row = local.row + DIRECTIONS[direction][0];
			out.col = local.col + DIRECTIONS[direction][1];
			return out.row >= 0 && out.row < Size && out.col >= 0 && out.col < Size;
		} else {
			out = local;
			return true;
		}
	}

	void flipsides()
	{
		for (uint idx = 0; idx < Size * Size; ++idx) {
			unsigned char field = map[idx];

			if ((field ^ (field >> 1)) & 1)
				field ^= Mine | Enemy;

			map[idx] = field;
		}
	}

	uint manhattandist(const Location & a, const Location & b) const {
		return abs(a.row - b.row) + abs(a.col - b.col);
	}
};

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
					out << 'A';
				else
					out << 'B';
			} else if (field & Submap::Ant) {
				if (field & Submap::Mine)
					out << 'a';
				else
					out << 'b';
			} else {
				if (field & Submap::Mine)
					out << '.';
				else if (field & Submap::Enemy)
					out << ',';
				else
					out << ' ';
			}
		}
		out << endl;
	}
	return out;
}

Tactical::Tactical(Bot & bot_) :
	bot(bot_),
	state(bot_.state)
{
}

void Tactical::gensubmap_field(Submap & sm, const Location & local, const Location & global)
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
	if ((field & (Submap::Enemy | Submap::Mine)) == 0) {
		int diff = bot.m_zoc.m_enemy[global] - bot.m_zoc.m_me[global];
		if (diff > 0)
			field |= Submap::Mine;
		else
			field |= Submap::Enemy;
	}
}

/**
 * Generate a tactical map centered at \p center
 */
void Tactical::gensubmap(Submap & sm, const Location & center)
{
	sm.fill(Submap::Water);

	Location queue[Submap::Size * Submap::Size];
	uint queue_tail = 0;
	uint queue_head = 0;

	gensubmap_field(sm, Location(Submap::Radius, Submap::Radius), center);
	queue[queue_head++] = Location(Submap::Radius, Submap::Radius);

	while (queue_tail < queue_head) {
		Location toplocal = queue[queue_tail++];

		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location nlocal;
			if (!sm.getneighbour(toplocal, dir, nlocal))
				continue;

			Location nglobal
				((center.row + nlocal.row + state.rows - Submap::Radius) % state.rows,
				 (center.col + nlocal.col + state.cols - Submap::Radius) % state.cols);

			if (state.grid[nglobal.row][nglobal.col].isWater || sm[nlocal] != Submap::Water)
				continue;

			gensubmap_field(sm, nlocal, nglobal);
			queue[queue_head++] = nlocal;
		}
	}
}

/**
 * Evaluate the given submap according to some cheap evaluation function.
 * Positive values are better for me, but negative values are quite untypical.
 *
 * Currently only based on field ownership. Could (should) be extended to
 * consider formations etc.
 */
int Tactical::evaluate(const Submap & sm)
{
	Submap potential;
	potential.fill(0);

	Location cur;
	for (cur.row = 0; cur.row < Submap::Size; ++cur.row) {
		for (cur.col = 0; cur.col < Submap::Size; ++cur.col) {
			unsigned char field = sm[cur];
			if (field & Submap::Food) {
				Location sub;
				for (sub.row = 0; sub.row <= 2 * FoodPotentialRadius; ++sub.row) {
					for (sub.col = 0; sub.col <= 2 * FoodPotentialRadius; ++sub.col) {
						Location rel
							(cur.row + sub.row - FoodPotentialRadius,
							 cur.col + sub.col - FoodPotentialRadius);

						if
							(rel.row >= 0 && rel.row < Submap::Size &&
							 rel.col >= 0 && rel.col < Submap::Size)
						{
							potential[rel] += FoodPotential[sub.row][sub.col];
						}
					}
				}
			}
		}
	}

	int value = 0;
	for (uint idx = 0; idx < Submap::Size * Submap::Size; ++idx) {
		int fieldvalue = 0;
		unsigned char field = sm.map[idx];

		if (!(field & Submap::Water))
			fieldvalue += ValueField;
		if (field & Submap::Food)
			fieldvalue += ValueFood;
		if (field & Submap::Ant)
			fieldvalue += potential.map[idx];

		if (field & Submap::Enemy)
			value -= fieldvalue;
		else if (field & Submap::Mine)
			value += fieldvalue;
	}

	return value;
}

struct AntMove {
	Location pos;
	int direction;
};

struct PlayerMove {
	char name[32];
	vector<AntMove> antmoves;
	bool refined;
	int worstvalue;
	int worstopposingidx;

	PlayerMove(const char * fmt, ...)
	{
		va_list va;
		va_start(va, fmt);
		vsnprintf(name, sizeof(name), fmt, va);
		va_end(va);

		refined = false;
		worstvalue = numeric_limits<int>::max();
		worstopposingidx = -1;
	}
};

struct Scenarios {
	Submap basesm;
	vector<PlayerMove> mymoves;
	vector<PlayerMove> enemymoves;
	uint myevaluated;
	uint enemyevaluated;

	Scenarios() {
		myevaluated = 0;
		enemyevaluated = 0;
	}

	void flipsides()
	{
		basesm.flipsides();
		swap(mymoves, enemymoves);
		swap(myevaluated, enemyevaluated);
	}

	uint mybestmove()
	{
		uint bestidx = 0;
		int bestvalue = numeric_limits<int>::min();

		for (uint idx = 0; idx < mymoves.size(); ++idx) {
			int value = mymoves[idx].worstvalue;
			if (value > bestvalue) {
				bestvalue = value;
				bestidx = idx;
			}
		}

		return bestidx;
	}

	uint enemybestmove()
	{
		uint bestidx = 0;
		int bestvalue = numeric_limits<int>::min();

		for (uint idx = 0; idx < enemymoves.size(); ++idx) {
			int value = enemymoves[idx].worstvalue;
			if (value > bestvalue) {
				bestvalue = value;
				bestidx = idx;
			}
		}

		return bestidx;
	}
};

void Tactical::generate_moves(Scenarios & scn, const Location & offset, bool predefined)
{
	state.bug << "generate_moves" << endl;
	state.bug << scn.basesm;

	// collect ants and predefined tactical moves
	vector<AntMove> basemoves;

	Location local;
	for (local.row = 0; local.row < Submap::Size; ++local.row) {
		for (local.col = 0; local.col < Submap::Size; ++local.col) {
			if ((scn.basesm[local] & (Submap::Ant | Submap::Mine)) == (Submap::Ant | Submap::Mine)) {
				basemoves.push_back(AntMove());
				AntMove & am = basemoves.back();
				am.pos = local;
				am.direction = -1;

				if (predefined) {
					Location global
						((local.row + offset.row) % state.rows,
						 (local.col + offset.col) % state.cols);
					Ant & ant = bot.m_ants[bot.myantidx_at(global)];
					if (ant.hastactical)
						am.direction = ant.direction;
				}
			}
		}
	}

	if (basemoves.empty()) {
		state.bug << "No initial tactical moves generated" << endl;
		return;
	}

	// now setup some standard moves
	{
		scn.mymoves.push_back(PlayerMove("stay"));
		PlayerMove & pm = scn.mymoves.back();
		pm.antmoves = basemoves;
	}

	int centermostidx = 0;
	uint centermostdist = numeric_limits<uint>::max();
	int centermostidleidx = -1;
	uint centermostidledist = numeric_limits<uint>::max();

	for (uint idx = 0; idx < basemoves.size(); ++idx) {
		const AntMove & am = basemoves[idx];
		uint dist = scn.basesm.manhattandist(am.pos, Location(Submap::Radius, Submap::Radius));

		if (dist < centermostdist) {
			centermostdist = dist;
			centermostidx = 0;
		}
		if (am.direction < 0 && dist < centermostdist) {
			centermostidledist = dist;
			centermostidleidx = idx;
		}
	}

	uint antidx = (centermostidleidx >= 0) ? centermostidleidx : centermostidx;

	for (int dir = 0; dir < TDIRECTIONS; ++dir) {
		Location pos = basemoves[antidx].pos;
		Location goal;
		if (!scn.basesm.getneighbour(pos, dir, goal))
			continue;

		if (scn.basesm[goal] & Submap::Water)
			continue;

		scn.mymoves.push_back(PlayerMove("midmove %c", cdir(dir)));
		PlayerMove & pm = scn.mymoves.back();
		pm.antmoves = basemoves;
		pm.antmoves[antidx].direction = dir;
	}
}

struct Outcome {
	int myvalue;
	int enemyvalue;

	Outcome() : myvalue(0), enemyvalue(0) {}
};

/**
 * Look at the outcome of myidx move vs enemyidx move, and try to generate some better
 * alternative moves for myself.
 */
void Tactical::improve(Scenarios & scn, uint myidx, uint enemyidx)
{
	PlayerMove & mymove = scn.mymoves[myidx];
	//PlayerMove & enemymove = scn.enemymoves[myidx];

	if (mymove.refined)
		return;

	mymove.refined = true;
}

void Tactical::apply_moves(Submap & sm, const PlayerMove & me, const PlayerMove & enemy, Outcome & outcome)
{
	// place moves
	for (uint amidx = 0; amidx < me.antmoves.size(); ++amidx) {
		const AntMove & am = me.antmoves[amidx];
		Location newpos;
		if (!sm.getneighbouropt(am.pos, am.direction, newpos))
			continue;

		unsigned char & field = sm[newpos];
		if ((field & (Submap::AntCollision | Submap::AntMoved)) == 0) {
			field &= ~Submap::Enemy;
			field |= Submap::AntMoved | Submap::Mine;
		} else {
			if (!(field & Submap::AntCollision)) {
				field |= Submap::AntCollision;
				outcome.myvalue += ValueLoss;
			}
			outcome.myvalue += ValueLoss;
		}
	}

	for (uint amidx = 0; amidx < enemy.antmoves.size(); ++amidx) {
		const AntMove & am = enemy.antmoves[amidx];
		Location newpos;
		if (!sm.getneighbouropt(am.pos, am.direction, newpos))
			continue;

		unsigned char & field = sm[newpos];
		if ((field & (Submap::AntCollision | Submap::AntMoved)) == 0) {
			field &= ~Submap::Mine;
			field |= Submap::AntMoved | Submap::Enemy;
		} else {
			if (!(field & Submap::AntCollision)) {
				field |= Submap::AntCollision;
				outcome.enemyvalue += ValueLoss;
			}
			outcome.enemyvalue += ValueLoss;
		}
	}

	// update ant positions
	Location newpositions[Submap::Size * Submap::Size];
	uint nrnewpositions = 0;

	Location cur;
	for (cur.row = 0; cur.row < Submap::Size; ++cur.row) {
		for (cur.col = 0; cur.col < Submap::Size; ++cur.col) {
			unsigned char & field = sm[cur];
			if ((field & (Submap::AntMoved | Submap::AntCollision)) == Submap::AntMoved) {
				field |= Submap::Ant;

				newpositions[nrnewpositions++] = cur;
			} else {
				field &= ~Submap::Ant;
			}
			field &= ~(Submap::AntMoved | Submap::AntCollision);
		}
	}

	// compute number of nearby enemies
	Submap enemies;

	for (uint idx = 0; idx < nrnewpositions; ++idx) {
		Location cur = newpositions[idx];
		unsigned char ownership = sm[cur] & (Submap::Mine | Submap::Enemy);

		int count = 0;
		Location kernel;
		for (kernel.row = 0; kernel.row < AttackKernelSize; ++kernel.row) {
			for (kernel.col = 0; kernel.col < AttackKernelSize; ++kernel.col) {
				if (!AttackKernel[kernel.row][kernel.col])
					continue;

				Location other
					((cur.row + kernel.row + state.rows - AttackKernelRadius) % state.rows,
					 (cur.col + kernel.col + state.cols - AttackKernelRadius) % state.cols);
				if
					(other.row < 0 || other.row >= Submap::Size ||
					 other.col < 0 || other.col >= Submap::Size)
					continue;

				unsigned char otherfield = sm[other];
				if (otherfield & Submap::Ant) {
					if ((otherfield & (Submap::Mine | Submap::Enemy)) != ownership)
						count++;
				}
			}
		}

		enemies[cur] = count;
	}

	// compute killed ants
	Location killants[Submap::Size * Submap::Size];
	uint nrkillants = 0;

	for (uint idx = 0; idx < nrnewpositions; ++idx) {
		Location cur = newpositions[idx];
		unsigned char ownership = sm[cur] & (Submap::Mine | Submap::Enemy);

		Location kernel;
		for (kernel.row = 0; kernel.row < AttackKernelSize; ++kernel.row) {
			for (kernel.col = 0; kernel.col < AttackKernelSize; ++kernel.col) {
				if (!AttackKernel[kernel.row][kernel.col])
					continue;

				Location other
					((cur.row + kernel.row + state.rows - AttackKernelRadius) % state.rows,
					 (cur.col + kernel.col + state.cols - AttackKernelRadius) % state.cols);
				if
					(other.row < 0 || other.row >= Submap::Size ||
					 other.col < 0 || other.col >= Submap::Size)
					continue;

				unsigned char otherfield = sm[other];
				if (!(otherfield & Submap::Ant) || (otherfield & (Submap::Mine | Submap::Enemy)) == ownership)
					continue;

				if (enemies[cur] >= enemies[other]) {
					killants[nrkillants++] = cur;
					goto killed;
				}
			}
		}
	killed: ;
	}

	// reap ants
	for (uint idx = 0; idx < nrkillants; ++idx) {
		Location cur = killants[idx];
		unsigned char & field = sm[cur];

		if (field & Submap::Mine) {
			outcome.myvalue += ValueLoss;
			outcome.enemyvalue += ValueKill;
		} else {
			outcome.myvalue += ValueKill;
			outcome.enemyvalue += ValueLoss;
		}

		field &= ~(Submap::Ant | Submap::AntCollision | Submap::AntMoved);
	}

	// eat food
	for (uint idx = 0; idx < nrnewpositions; ++idx) {
		Location cur = killants[idx];
		unsigned char field = sm[cur];
		if (!(field & Submap::Ant))
			continue;

		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n;
			if (!sm.getneighbour(cur, dir, n))
				continue;

			unsigned char & nfield = sm[cur];
			if (nfield & Submap::Food) {
				nfield &= ~Submap::Food;

				if (field & Submap::Enemy)
					outcome.enemyvalue += ValueEat;
				else
					outcome.myvalue += ValueEat;
			}
		}
	}

	// clear field ownership of interior fields
	for (cur.row = 1; cur.row < Submap::Size - 1; ++cur.row) {
		for (cur.col = 1; cur.col < Submap::Size - 1; ++cur.col) {
			unsigned char & field = sm[cur];

			if ((field & (Submap::Ant | Submap::Hill)) == 0)
				field &= ~(Submap::Mine | Submap::Enemy);
		}
	}

	// rebuild field ownership
	Location queue[Submap::Size * Submap::Size];
	uint queue_head = 0;
	uint queue_tail = 0;

#define update_neighbour_ownership(cur) \
	do { \
		unsigned char ownership = sm[cur] & (Submap::Mine | Submap::Enemy); \
		if (!ownership) break; \
		Location pos(cur); \
		for (int dir = 0; dir < TDIRECTIONS; ++dir) { \
			Location n; \
			if (!sm.getneighbour(pos, dir, n)) continue; \
			unsigned char & nf = sm[n]; \
			if (nf & (Submap::Mine | Submap::Enemy)) continue; \
			nf |= sm[pos] & (Submap::Mine | Submap::Enemy); \
			queue[queue_tail++] = n; \
		} \
	} while(0)

	for (uint pos = 1; pos < Submap::Size - 1; ++pos) {
		cur.col = pos;
		cur.row = 0;
		update_neighbour_ownership(cur);

		cur.row = Submap::Size - 1;
		update_neighbour_ownership(cur);

		cur.row = pos;
		cur.col = 0;
		update_neighbour_ownership(cur);

		cur.col = Submap::Size;
		update_neighbour_ownership(cur);
	}

	for (uint idx = 0; idx < nrnewpositions; ++idx) {
		Location cur = newpositions[idx];

		if (sm[cur] & Submap::Ant)
			update_neighbour_ownership(cur);
	}

	while (queue_head < queue_tail) {
		Location cur = queue[queue_head++];
		update_neighbour_ownership(cur);
	}

#undef update_neighbour_ownership
}

bool Tactical::evaluate_new_moves(Scenarios & scn)
{
	if (scn.myevaluated >= scn.mymoves.size() && scn.enemyevaluated >= scn.enemymoves.size())
		return false;

	for (uint myidx = 0; myidx < scn.mymoves.size(); ++myidx) {
		PlayerMove & mymove = scn.mymoves[myidx];

		for
			(uint enemyidx = (myidx < scn.myevaluated) ? scn.enemyevaluated : 0;
			 enemyidx < scn.enemymoves.size(); ++enemyidx)
		{
			Outcome outcome;
			PlayerMove & enemymove = scn.enemymoves[enemyidx];
			Submap sm(scn.basesm);
			apply_moves(sm, mymove, enemymove, outcome);

			int value = evaluate(sm);
			outcome.myvalue += value;
			outcome.enemyvalue -= value;

			state.bug << " " << myidx << " vs. " << enemyidx
				<< "  value " << outcome.myvalue << " vs " << outcome.enemyvalue << endl;
			state.bug << sm;

			if (outcome.myvalue < mymove.worstvalue) {
				mymove.worstvalue = outcome.myvalue;
				mymove.worstopposingidx = enemyidx;
			}
			if (outcome.enemyvalue < enemymove.worstvalue) {
				enemymove.worstvalue = outcome.enemyvalue;
				enemymove.worstopposingidx = myidx;
			}
		}
	}

	scn.myevaluated = scn.mymoves.size();
	scn.enemyevaluated = scn.enemymoves.size();
	return true;
}

void Tactical::make_moves(const Location & center)
{
	Scenarios scn;

	gensubmap(scn.basesm, center);

	state.bug << "Tactical around " << center << endl;
	state.bug << scn.basesm;

	Location offset
		((center.row - Submap::Radius + state.rows) % state.rows,
		 (center.col - Submap::Radius + state.cols) % state.cols);

	scn.flipsides();
	generate_moves(scn, offset, false);
	if (scn.mymoves.empty())
		return;

	scn.flipsides();
	generate_moves(scn, offset, true);
	assert(!scn.mymoves.empty());

	for (;;) {
		if (!evaluate_new_moves(scn))
			break;

		if (scn.mymoves.size() + scn.enemymoves.size() >= MaxEval)
			break;

		uint mybestidx = scn.mybestmove();
		uint mybestenemyidx = scn.mymoves[mybestidx].worstopposingidx;

		uint enemybestidx = scn.enemybestmove();
		uint enemybestmyidx = scn.enemymoves[enemybestidx].worstopposingidx;

		improve(scn, mybestidx, mybestenemyidx);
		scn.flipsides();
		improve(scn, mybestenemyidx, mybestidx);
		improve(scn, enemybestidx, enemybestmyidx);
		scn.flipsides();
	}

	uint bestmoveidx = scn.mybestmove();
	const PlayerMove & move = scn.mymoves[bestmoveidx];

	state.bug << "Best move " << bestmoveidx << " name " << move.name << endl;

	for (uint idx = 0; idx < move.antmoves.size(); ++idx) {
		const AntMove & antmove = move.antmoves[idx];

		if (antmove.direction < 0)
			continue;

		Location global
			((antmove.pos.row + offset.row) % state.rows,
			 (antmove.pos.col + offset.col) % state.cols);
		uint antidx = bot.myantidx_at(global);
		Ant & ant = bot.m_ants[antidx];

		ant.hastactical = true;
		ant.direction = antmove.direction;
		ant.goal = state.getLocation(ant.where, antmove.direction);

		state.bug << "  Ant at " << ant.where
			<< " tactical move to " << ant.goal << " (" << cdir(ant.direction) << ")" << endl;
	}
}
