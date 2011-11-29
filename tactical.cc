#include "tactical.h"

#include <cassert>
#include <cstring>
#include <cstdarg>
#include <cstdio>
#include <stdlib.h>
#include <limits>
#include <algorithm>

#include "Bot.h"
#include "State.h"

using namespace std;

/// Tweakable parameters
//@{
static const int MaxFoodDistance = 16;
static const int FoodFactor = 10;
static const int MaxFoodSeekers = 2;
static const int MaxZocValue = 16;
static const int ZocFactor = 1;

static const int ValueKill = 5000;
static const int ValueLoss = -6000;
static const int ValueEat = 1000;

static const uint MaxEval = 50; ///< max number of moves to evaluate
static const uint MaxIdleDist2 = 10;

static const uint TacticalCoreRadius2 = 26;
static const uint TacticalMidRadius2 = 42;
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

static const int AttackRadius2 = 5;
static const int AttackNeighboursRadius2 = 10;
static const int NrAttackNeighboursUpperBound = 49;

template<typename T>
struct BaseSubmap {
	static const int Radius = 8;
	static const int Size = 2 * Radius + 1;

	T map[Size * Size];

	void fill(T c) {
		for (uint idx = 0; idx < Size * Size; ++idx)
			map[idx] = c;
	}
	T & operator[](const Location & pos) {return map[pos.row * Size + pos.col];}
	const T & operator[](const Location & pos) const {return map[pos.row * Size + pos.col];}
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
	uint manhattandist(const Location & a, const Location & b) const {
		return abs(a.row - b.row) + abs(a.col - b.col);
	}

	uint eucliddist2(const Location & a, const Location & b) const {
		int dr = a.row - b.row;
		int dc = a.col - b.col;
		return (dr * dr) + (dc * dc);
	}
};

struct Submap : BaseSubmap<uint8_t> {
	static const uint8_t Water = 0x20;
	static const uint8_t Food = 0x10;
	static const uint8_t Hill = 0x08;
	static const uint8_t Ant = 0x04;
	static const uint8_t Mine = 0x02; // must be the two least significant bits
	static const uint8_t Enemy = 0x01;

	void flipsides()
	{
		for (uint idx = 0; idx < Size * Size; ++idx) {
			uint8_t & field = map[idx];

			if ((field ^ (field >> 1)) & 1)
				field ^= Mine | Enemy;
		}
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
				out << ' ';
			}
		}
		out << endl;
	}
	return out;
}

static void add_enemy_kernel(Submap & sm, const Location & center, int delta)
{
	Location kernel;
	for (kernel.row = 0; kernel.row < AttackKernelSize; ++kernel.row) {
		for (kernel.col = 0; kernel.col < AttackKernelSize; ++kernel.col) {
			if (!AttackKernel[kernel.row][kernel.col])
				continue;

			Location other
				((center.row + kernel.row - AttackKernelRadius),
				 (center.col + kernel.col - AttackKernelRadius);
			if
				(other.row < 0 || other.row >= Submap::Size ||
				 other.col < 0 || other.col >= Submap::Size)
				continue;

			sm[other] += delta;
		}
	}
}

struct PlayerMove {
	struct AntMove {
		Location pos;
		int8_t direction;
		uint8_t collided : 1; // death by collision
		uint8_t collider : 1; // we ran into the other guy; oops

		uint8_t override : 1;
		uint8_t killed : 1;
		uint8_t sentoffense : 1;

		AntMove(const Location & p, int8_t dir) : pos(p), direction(dir), collided(0), collider(0), override(0) {}
	};

	vector<AntMove> antmoves;
	uint hash;
	uint nrcollided;

	static uint8_t AntsMask = 0xc0;
	static uint8_t AntsShift = 5;
	static uint8_t AttackMask = 0x1f;
	Submap map;

	int worstvalue;
	int worstopposingidx;

	PlayerMove() {
		hash = 0;
		worstvalue = numeric_limits<int>::max();
		worstopposingidx = -1;
	}

	void computehash() {
		hash = 0;
		for (uint idx = 0; idx < antmoves.size(); ++idx)
			hash = (hash * 5) + antmoves[idx].direction;
	}

	void ant_mark(uint idx) {
		// assume: direction and pos already set
		AntMove & am = antmoves[idx];

		if (am.pos.row < 0 || am.pos.row >= Submap::Size || am.pos.col < 0 || am.pos.col >= Submap::Size)
			return;

		uint prevnrants = map[am.pos] >> AntsShift;

		map[am.pos] += 1 << AntsShift;
		if (prevnrants == 0) {
			add_enemy_kernel(map, am.pos, 1);
			am.collided = 0;
			am.collider = 0;
		} else if (prevnrants == 1) {
			add_enemy_kernel(map, am.pos, -1);
			for (uint other = 0; other < antmoves.size(); ++other) {
				AntMove & otherm = antmoves[other];
				if (other == idx || other.pos != otherm.pos)
					continue;
				otherm.collided = 1;
			}
			am.collided = 1;
			am.collider = 1;
			nrcollided += 2;
		} else {
			am.collider = 1;
			am.collided = 1;
			nrcollided++;
		}
	}

	void ant_unmark(uint idx) {
		// assume: direction and pos still set
		AntMove & am = antmoves[idx];

		if (am.pos.row < 0 || am.pos.row >= Submap::Size || am.pos.col < 0 || am.pos.col >= Submap::Size)
			return;

		uint prevnrants = map[am.pos] >> AntsShift;
		assert(prevnrants >= 1);

		map[am.pos] -= 1 << AntsShift;
		if (prevnrants == 1) {
			add_enemy_kernel(map, am.pos, -1);
		} else if (prevnrants == 2) {
			for (uint other = 0; other < antmoves.size(); ++other) {
				AntMove & otherm = antmoves[other];
				if (other == idx || other.pos != otherm.pos)
					continue;

				add_enemy_kernel(map, am.pos, 1);
				otherm.collided = 0;
				otherm.collider = 0;
				nrcollided -= 2;
				break;
			}
		} else {
			if (!am.collider) {
				for (uint other = 0; other < antmoves.size(); ++other) {
					AntMove & otherm = antmoves[other];
					if (other == idx || other.pos != otherm.pos)
						continue;

					otherm.collider = 0;
					break;
				}
			}
			nrcollided--;
		}
	}
};


struct Tactical::Data {
	bool ismyperspective;
	Location offset;
	Submap basesm;
	vector<Location> myants;
	vector<Location> enemyants;
	BaseSubmap<unsigned int> foodpotential;

	vector<PlayerMove> mymoves;
	vector<PlayerMove> enemymoves;
	uint myevaluated;
	uint enemyevaluated;

	void flipsides()
	{
		ismyperspective = !ismyperspective;
		basesm.flipsides();
		swap(myants, enemyants);
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

Tactical::Tactical(Bot & bot_) :
	bot(bot_),
	state(bot_.state),
	d(*new Data)
{
}

Tactical::~Tactical()
{
	delete *d;
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
}

/**
 * Generate a tactical map centered at \p center
 */
void Tactical::gensubmap(Submap & sm, const Location & center)
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

void Tactical::compute_foodpotential(const Submap & sm, BaseSubmap<uint> & potential)
{
	potential.fill(0);

	Location food;
	for (food.row = 0; food.row < Submap::Size; ++food.row) {
		for (food.col = 0; food.col < Submap::Size; ++food.col) {
			if (!(sm[food] & Submap::Food))
				continue;

			Location queue[Submap::Size * Submap::Size];
			uint queue_head = 0;
			uint queue_tail = 0;
			Submap thispotential;

			thispotential.fill(0);
			thispotential[food] = MaxFoodDistance;
			queue[queue_tail++] = food;

			while (queue_head < queue_tail) {
				Location cur = queue[queue_head++];
				potential[cur] += thispotential[cur] * FoodFactor;

				if (thispotential[cur] > 1) {
					for (int dir = 0; dir < TDIRECTIONS; ++dir) {
						Location n;
						if (!sm.getneighbour(cur, dir, n))
							continue;
						if (sm[n] & Submap::Water || thispotential[n] > 0)
							continue;

						thispotential[n] = thispotential[cur] - 1;
						queue[queue_tail++] = n;
					}
				}
			}
		}
	}
}

void Tactical::evaluate_food(const Submap & sm, int & myvalue, int & enemyvalue)
{
	Location food;
	for (food.row = 0; food.row < Submap::Size; ++food.row) {
		for (food.col = 0; food.col < Submap::Size; ++food.col) {
			if (!(sm[food] & Submap::Food))
				continue;

			Location queue[Submap::Size * Submap::Size];
			uint queue_head = 0;
			uint queue_tail = 0;
			Submap thispotential;
			uint myfound = 0;
			uint enemyfound = 0;

			thispotential.fill(0);
			thispotential[food] = MaxFoodDistance;
			queue[queue_tail++] = food;

			while (queue_head < queue_tail) {
				Location cur = queue[queue_head++];

				if (sm[cur] & Submap::Ant) {
					if (sm[cur] & Submap::Enemy) {
						enemyfound++;
						enemyvalue += thispotential[cur] * FoodFactor / enemyfound;
					} else {
						myfound++;
						myvalue += thispotential[cur] * FoodFactor / myfound;
					}
					if (enemyfound >= MaxFoodSeekers && myfound >= MaxFoodSeekers)
						break;
				}

				if (thispotential[cur] > 1) {
					for (int dir = 0; dir < TDIRECTIONS; ++dir) {
						Location n;
						if (!sm.getneighbour(cur, dir, n))
							continue;
						if (sm[n] & Submap::Water || thispotential[n] > 0)
							continue;

						thispotential[n] = thispotential[cur] - 1;
						queue[queue_tail++] = n;
					}
				}
			}
		}
	}
}


/**
 * Evaluate the given submap according to some simple and hopefully cheap evaluation function.
 */
void Tactical::evaluate(const Submap & sm, int & myvalue, int & enemyvalue)
{
	evaluate_food(sm, myvalue, enemyvalue);

	Location local;
	for (local.row = 0; local.row < Submap::Size; ++local.row) {
		for (local.col = 0; local.col < Submap::Size; ++local.col) {
			if (sm[local] & Submap::Ant) {
				Location global
					((local.row + d.offset.row) % state.rows,
					 (local.col + d.offset.col) % state.cols);

				if (sm[local] & Submap::Enemy) {
					enemyvalue += (MaxZocValue - min(bot.m_zoc.m_me[global], MaxZocValue)) * ZocFactor;
				} else {
					myvalue += (MaxZocValue - min(bot.m_zoc.m_enemy[global], MaxZocValue)) * ZocFactor;
				}
			}
		}
	}
}

struct CloserToCenterCompare {
	bool operator()(const Location & a, const Location & b) const {
		uint a2 = Submap::eucliddist2(a, Location(Submap::Radius, Submap::Radius);
		uint b2 = Submap::eucliddist2(b, Location(Submap::Radius, Submap::Radius);
		return a2 < b2;
	}
};

void Tactical::generate_initmoves()
{
	d.mymoves.push_back(PlayerMove());
	d.enemymoves.push_back(PlayerMove());

	// collect all ants and sort by distance to center
	Location local;
	for (local.row = 0; local.row < Submap::Size; ++local.row) {
		for (local.col = 0; local.col < Submap::Size; ++local.col) {
			unsigned char field = d.basesm[local];

			if (!(field & Submap::Ant))
				continue;

			if (field & Submap::Enemy) {
				d.enemyants.push_back(local);
			} else {
				d.myants.push_back(local);
			}
		}
	}

	sort(d.myants.begin(), d.myants.end(), CloserToCenterCompare());
	sort(d.enemyants.begin(), d.enemyants.end(), CloserToCenterCompare());

	// Fill in default moves
	d.mymoves[0].antmoves.reserve(d.myants.size());
	for (uint idx = 0; idx < d.myants.size(); ++idx) {
		const Location & local = d.myants[idx];
		Location global
			((local.row + d.offset.row) % state.rows,
			 (local.col + d.offset.col) % state.cols);
		Ant & ant = bot.m_ants[bot.myantidx_at(global)];
		Location nextlocal;
		d.basesm.getneighbouropt(local, ant.direction, nextlocal);
		d.mymoves[0].antmoves.push_back(PlayerMove::AntMove(nextlocal, ant.direction));
		d.mymoves[0].ant_mark(idx);
	}

	d.enemymoves[0].antmoves.reserve(d.enemyants.size());
	for (uint idx = 0; idx < d.enemyants.size(); ++idx) {
		const Location & local = d.enemyants[idx];
		d.enemymoves[0].antmoves.push_back(PlayerMove::AntMove(local, -1));
		d.enemymoves[0].ant_mark(idx);
	}

	// Finalize the moves
	d.mymoves[0].computehash();
	d.enemymoves[0].computehash();
}


struct Outcome {
	int myvalue;
	int enemyvalue;

	Location newpositions[Submap::Size * Submap::Size];
	uint nrnewpositions;

	Outcome() : myvalue(0), enemyvalue(0) {}
};

struct NewMove {
	NewMove(Tactical & t_, uint myidxbase, const char * why) :
		t(t_), d(t_.d), state(t_.state), shouldcommit(false)
	{
		state.bug << "Preliminary insert new move " << d.mymoves.size() << " due to " << why
			<< " perspective " << (d.ismyperspective ? "mine" : "enemy") << endl;
		d.mymoves.push_back(PlayerMove());
		d.mymoves.back() = d.mymoves[myidxbase];

		PlayerMove & pm = d.mymoves.back();
		pm.worstvalue = numeric_limits<int>::max();
		pm.worstopposingidx = -1;
	}

	~NewMove()
	{
		if (!shouldcommit)
			d.mymoves.pop_back();
	}

	PlayerMove & move() {return d.mymoves.back();}

	void antmove(uint myantidx, int direction)
	{
		PlayerMove & pm = move();
		pm.ant_unmark(myantidx);

		PlayerMove::AntMove & am = pm.antmoves[myantidx];
		am.direction = direction;
		d.basesm.getneighbouropt(d.myants[myantidx], direction, am.pos);
		pm.ant_mark(myantidx);
	}

	void commit()
	{
		assert(!shouldcommit);

		PlayerMove & pm = move();
		pm.computehash();

		for (uint idx = 0; idx < d.mymoves.size() - 1; ++idx) {
			PlayerMove & other = d.mymoves[idx];
			if (other.hash != pm.hash)
				continue;

			for (uint antidx = 0; antidx < pm.antmoves.size(); ++antidx) {
				if (pm.antmoves[antidx].direction != other.antmoves[antidx].direction)
					goto notequal;
			}

			state.bug << "New move is equal to old move " << idx << endl;

			return;
		notequal: ;
		}

		shouldcommit = true;
	}

	Tactical & t;
	Tactical::Data & d;
	State & state;
	bool shouldcommit;
};

struct Improve {
	Tactical & t;
	Tactical::Data & d;
	State & state;
	uint myidx;
	uint enemyidx;

	struct Death {
		uint myantidx;
		uint enemyantidx;
	};

	vector<Death> deaths;

	Improve(Tactical & t_, uint myidx_, uint enemyidx_) :
		t(t_), d(t_.d), state(t_.state),
		myidx(myidx_), enemyidx(enemyidx_)
	{
	}

	void compute_deaths()
	{
		PlayerMove & mymove = themymove();
		PlayerMove & enemymove = d.enemymoves[enemyidx];

		for (uint myantidx = 0; myantidx < mymove.antmoves.size(); ++myantidx) {
			PlayerMove::AntMove & am = mymove.antmoves[myantidx];
			if (am.collided)
				continue;

			am.sentoffense = 0;
			am.killed = 0;

			uint nrenemies = enemymove.map[am.pos] & PlayerMove::AttackMask;
			if (nrenemies == 0)
				continue;

			for (uint enemyantidx = 0; enemyantidx < enemymove.antmoves.size(); ++enemyantidx) {
				const PlayerMove::AntMove & em = enemymove.antmoves[enemyantidx];
				if (em.collided)
					continue;

				Location kernel
					(em.pos.row - am.pos.row + AttackKernelRadius,
					 em.pos.col - am.pos.col + AttackKernelRadius);

				if ((mymove.map[em.pos] & PlayerMove::AttackMask) > nrenemies)
					continue;
				if
					(kernel.row < 0 || kernel.row >= AttackKernelSize ||
					 kernel.col < 0 || kernel.col >= AttackKernelSize ||
					 !AttackKernel[kernel.row][kernel.col])
					continue;

				am.killed = 1;

				// ignore deaths that are too far from the center, because we cannot assess them accurately
				if (d.basesm.eucliddist2(am.pos, Location(Submap::Radius, Submap::Radius)) > TacticalCoreRadius2)
					break;

				deaths.push_back(Death());
				deaths.back().myantidx = myantidx;
				deaths.back().enemyantidx = enemyantidx;
				break;
			}
		}

		for (uint enemyantidx = 0; enemyantidx < enemymove.antmoves.size(); ++enemyantidx) {
			PlayerMove::AntMove & am = enemymove.antmoves[enemyantidx];
			if (am.collided)
				continue;

			am.sentoffense = 0;
			am.killed = 0;

			uint nrenemies = mymove.map[am.pos] & PlayerMove::AttackMask;
			if (nrenemies == 0)
				continue;

			for (uint myantidx = 0; myantidx < mymove.antmoves.size(); ++myantidx) {
				const PlayerMove::AntMove & em = mymove.antmoves[myantidx];
				if (em.collided)
					continue;

				Location kernel
					(em.pos.row - am.pos.row + AttackKernelRadius,
					 em.pos.col - am.pos.col + AttackKernelRadius);

				if ((enemymove.map[em.pos] & PlayerMove::AttackMask) > nrenemies)
					continue;
				if
					(kernel.row < 0 || kernel.row >= AttackKernelSize ||
					 kernel.col < 0 || kernel.col >= AttackKernelSize ||
					 !AttackKernel[kernel.row][kernel.col])
					continue;

				am.killed = 1;
				break;
			}
		}
	}

	void avoid_collisions()
	{
		NewMove nm(t, myidx, "collision avoidance");
		const PlayerMove & mymove = themymove();
		const PlayerMove & enemymove = theenemymove();

		for (uint myantidx = 0; myantidx < mymove.antmoves.size(); ++myantidx) {
			const PlayerMove::AntMove & am = mymove.antmoves[myantidx];
			if (!am.collider)
				continue;

#define checkacceptable(label) \
	do { \
		Location to; \
		if (!d.basesm.getneighbouropt(from, dir, to)) break; \
		if (d.basesm[to] & Submap::Water) break; \
		if (nm.move().map[to] & PlayerMove::AntsMask) break; \
		altdir = dir; \
		if ((enemymove.map[to] & PlayerMove::AttackMask) == 0) \
			goto label; \
	} while(0)

			Location from = d.myants[myantidx];
			int olddir = am.direction;
			if (olddir == -1) {
				int altdir = -1;
				int * dirperm = getdirperm();
				for (uint predir = 0; predir < TDIRECTIONS; ++predir) {
					int dir = dirperm[predir];
					checkacceptable(found1);
				}
			found1:
				nm.antmove(myantidx, altdir);
			} else {
				int altdir = olddir;
				int dir = (olddir + 1 + (fastrng() & 2)) % TDIRECTIONS;

				checkacceptable(found2);

				dir = (dir + 2) % TDIRECTIONS;
				checkacceptable(found2);

				dir = (olddir + 2) % TDIRECTIONS;
				checkacceptable(found2);

				dir = -1;
				checkacceptable(found2);
			found2:
				nm.antmove(myantidx, altdir);
			}
		}

#undef checkacceptable

		nm.commit();
	}

	// if one of ours is killed, see if we can kill the enemy instead
	void do_defense()
	{
		for (uint idx = 0; idx < deaths.size(); ++idx) {
			uint enemyantidx = deaths[idx].enemyantidx;
			if (!theenemymove().antmoves[enemyantidx].sentoffense)
				send_to_kill(enemyantidx, "defense");
		}
	}

	void do_retreat()
	{
		for (uint idx = 0; idx < deaths.size(); ++idx) {
			const PlayerMove & mymove = themymove();
			const PlayerMove & enemymove = theenemymove();
			Death & death = deaths[idx];
			const Location & orig = d.myants[death.myantidx];
			const Location & other = enemymove.antmoves[death.enemyantidx].pos;

			int origdir = mymove.antmoves[death.myantidx].direction;
			int alt1dir = origdir;
			int alt2dir = origdir;
			int * optdirperm = getoptdirperm();
			for (int predir = 0; predir < TDIRECTIONS+1; ++predir) {
				int dir = optdirperm[predir];
				if (dir == origdir)
					continue;

				Location n;
				if (!d.basesm.getneighbouropt(orig, dir, n))
					continue;
				if (d.basesm[n] & Submap::Water)
					continue;

				if ((enemymove.map[n] & PlayerMove::AttackMask) == 0) {
					alt1dir = dir;
					if (alt2dir != origdir)
						break;
				} else {
					Location kernel
						(other.row - n.row + AttackKernelRadius,
						 other.col - n.col + AttackKernelRadius);
					if
						(kernel.row < 0 || kernel.row >= AttackKernelSize ||
						 kernel.col < 0 || kernel.col >= AttackKernelSize ||
						 !AttackKernel[kernel.row][kernel.col])
					{
						alt2dir = dir;
						if (alt1dir != origdir)
							break;
					}
				}
			}

			if (alt1dir != origdir) {
				NewMove nm(t, myidx, "retreat: unattacked field");
				if (nm.antmove(death.myantidx, alt1dir, NewMove::PushRetreat))
					nm.commit();
			}
			if (alt2dir != origdir && alt2dir != alt1dir) {
				NewMove nm(t, myidx, "retreat: get away from killer");
				if (nm.antmove(death.myantidx, alt2dir, NewMove::PushRetreat))
					nm.commit();
			}
		}
	}

	void send_to_kill(uint enemyantidx, const char * why)
	{
		const PlayerMove & mymove = themymove();
		const PlayerMove & enemymove = theenemymove();
		const PlayerMove::AntMove & enemyant = enemymove.antmoves[enemyantidx];
		Location enemypos = enemyant.pos;

		uint additionalidx[NrAttackNeighboursUpperBound];
		int additionaldir[NrAttackNeighboursUpperBound];
		uint additionalvalue[NrAttackNeighboursUpperBound];
		uint nradditional;

		uint nrattackers = mymove.map[enemypos] & PlayerMove::AttackMask;
		uint bestattacker = numeric_limits<uint>::max();

		uint p = getprime();
		uint ofs = fastrng() % d.myants.size();
		for (uint preidx = 0; preidx < d.myants.size(); ++preidx) {
			uint myantidx = (preidx * p + ofs) % d.myants.size();
			const Location & orig = d.myants[myantidx];
			const Location & current = mymove.antmoves[myantidx].pos;

			if (d.basesm.eucliddist2(current, enemypos) <= AttackRadius2) {
				// already attacking the target
				bestattacker = min(bestattacker, enemymove.map[current] & PlayerMove::AttackMask);
				continue;
			}

			if (d.basesm.eucliddist2(orig, enemypos) > AttackNeighboursRadius2)
				continue; // too far away

			int bestdirection = -1;
			uint bestvalue = numeric_limits<uint>::max();
			int * optdirperm = getoptdirperm();
			for (int predir = 0; predir < TDIRECTIONS + 1; ++predir) {
				int dir = optdirperm[predir];
				Location n;
				if (!d.basesm.getneighbouropt(orig, dir, n))
					continue;
				if (d.basesm.eucliddist2(n, enemypos) > AttackRadius2)
					continue;
				if (d.basesm[n] & Submap::Water)
					continue;

				uint value = enemymove.map[n] & PlayerMove::AttackMask;
				if (value < bestvalue) {
					bestvalue = value;
					bestdirection = dir;
				}
			}

			additionalidx[nradditional] = myantidx;
			additionaldir[nradditional] = bestdirection;
			additionalvalue[nradditional] = additionalvalue;
			nradditional++;
		}

		if (!nradditional)
			return;

		NewMove nm(t, myidx, why);
		for (uint idx = 0; idx < nradditional; ++idx) {
			if (nrattackers > bestattacker)
				break;

			uint myantidx = additionalidx[idx];
			const Location & orig = d.myants[myantidx];
			int dir = additionaldir[idx];
			Location n;
			if (!d.basesm.getneighbouropt(orig, dir, n)) {
				abort();
			}

			if (nm.move().map[n] & PlayerMove::AntsMask)
				continue;

			nm.antmove(myantidx, dir);
			bestattacker = min(bestattacker, additionalvalue[idx]);
		}

		if (nrattackers >= bestattacker) {
			theenemymove().antmoves[enemyantidx].sentoffense = 1;
			nm.commit();
		}
	}

	void do_offense()
	{
		for (uint enemyantidx = 0; enemyantidx < d.enemyants.size(); ++enemyantidx) {
			const PlayerMove & enemymove = theenemymove();
			const PlayerMove::AntMove & enemyant = enemymove.antmoves[enemyantidx];

			if (enemyant.collided || enemyant.killed || enemyant.sentoffense)
				continue;
			if (d.basesm.eucliddist2(enemyant.pos, Location(Submap::Radius, Submap::Radius)) > TacticalMidRadius2)
				continue;

			send_to_kill(enemyantidx, "offense");
		}
	}

	void do_improve()
	{
		t.state.bug << "improve " << myidx << " vs " << enemyidx
			<< " perspective: " << (d.ismyperspective ? "mine" : "enemy") << endl;

		//
		if (themymove().nrcollided) {
			avoid_collisions();
		}

		//
		compute_deaths();

		if (!deaths.empty()) {
			// if one of ours is killed, see if we can retreat
			do_retreat();
		}

		// try an unprovoked offense
		do_offense();
	}

	PlayerMove & themymove() const {return d.mymoves[myidx];}
	PlayerMove & theenemymove() const {return d.enemymoves[enemyidx];}
};

/**
 * Look at the outcome of myidx move vs enemyidx move, and try to generate some better
 * alternative moves for myself.
 */
void Tactical::improve(uint myidx, uint enemyidx)
{
	Improve imp(*this, myidx, enemyidx);
	imp.do_improve();
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
					(cur.row + kernel.row - AttackKernelRadius,
					 cur.col + kernel.col - AttackKernelRadius);
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
}

bool Tactical::evaluate_new_moves()
{
	if (d.myevaluated >= d.mymoves.size() && d.enemyevaluated >= d.enemymoves.size())
		return false;

	for (uint myidx = 0; myidx < d.mymoves.size(); ++myidx) {
		PlayerMove & mymove = d.mymoves[myidx];

		for
			(uint enemyidx = (myidx < d.myevaluated) ? d.enemyevaluated : 0;
			 enemyidx < d.enemymoves.size(); ++enemyidx)
		{
			PlayerMove & enemymove = d.enemymoves[enemyidx];
			int myvalue = mymove.nrcollided * ValueLoss;
			int enemyvalue = enemymove.nrcollided * ValueLoss;

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
	d.ismyperspective = true;
	gensubmap(d.basesm, center);
	d.offset = Location
		((center.row - Submap::Radius + state.rows) % state.rows,
		 (center.col - Submap::Radius + state.cols) % state.cols);

	d.myevaluated = 0;
	d.enemyevaluated = 0;

	state.bug << "Tactical around " << center << endl;
	state.bug << d.basesm;

	generate_initmoves();

	// iterative improvement of moves
	for (;;) {
		if (!evaluate_new_moves())
			break;

		if (d.mymoves.size() + d.enemymoves.size() >= MaxEval)
			break;

		uint mybestidx = d.mybestmove();
		uint mybestenemyidx = d.mymoves[mybestidx].worstopposingidx;

		uint enemybestidx = d.enemybestmove();
		uint enemybestmyidx = d.enemymoves[enemybestidx].worstopposingidx;

		improve(mybestidx, mybestenemyidx);
		d.flipsides();
		improve(mybestenemyidx, mybestidx);
		improve(enemybestidx, enemybestmyidx);
		d.flipsides();
	}

	// execute best move that was found
	uint bestmoveidx = d.mybestmove();
	const PlayerMove & move = d.mymoves[bestmoveidx];

	state.bug << "Best move " << bestmoveidx << endl;

	for (uint idx = 0; idx < move.antmoves.size(); ++idx) {
		int dir = move.antmoves[idx];
		Location antpos = d.myants[idx];
		Location global
			((antpos.row + d.offset.row) % state.rows,
			 (antpos.col + d.offset.col) % state.cols);
		uint antidx = bot.myantidx_at(global);
		Ant & ant = bot.m_ants[antidx];

		if (move.override[idx])
			ant.hastactical = true;
		ant.direction = -1;
		ant.goal = state.getLocation(ant.where, dir);

		state.bug << "  Ant at " << ant.where
			<< " tactical move to " << ant.goal << " (" << cdir(ant.direction) << ")" << endl;
	}
}
