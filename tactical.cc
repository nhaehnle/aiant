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
#include "zoc.h"

using namespace std;

/// Tweakable parameters
//@{
static const uint MaxEval = 1000; ///< max number of moves to evaluate

static const uint AttackableManhattanDistance = 5; ///< max distance for an ant to trigger theater creation
static const uint TheaterCoreInfinity = 2; ///< max infinity norm distance for an ant to "belong" to the theater
static const uint MaxAttackRadius2 = 34; ///< radius up to which we consider attacking an enemy (this is the maximum distance a core ant can possibly attack)

static const int MaxFoodDistance = 16;
static const int FoodFactor = 10;
static const int MaxZocValue = 16;
static const int ZocFactor = 1;
static const int HillValue = 300;

static const int ValueKill = 5000;
static const int ValueLoss = -6000;
static const int ValueEat = 1000;
//static const int ValueUselessDeath = -1500;
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

static const uint AttackRadius2 = 5;
static const uint AttackNeighboursRadius2 = 10;
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
	static bool inside(const Location & pos) {
		return pos.row >= 0 && pos.row < Size && pos.col >= 0 && pos.col < Size;
	}
	T & operator[](const Location & pos) {
		assert(inside(pos));
		return map[pos.row * Size + pos.col];
	}
	const T & operator[](const Location & pos) const {
		assert(inside(pos));
		return map[pos.row * Size + pos.col];
	}
	static bool getneighbour(Location local, int direction, Location & out) {
		out.row = local.row + DIRECTIONS[direction][0];
		out.col = local.col + DIRECTIONS[direction][1];
		return inside(out);
	}
	static bool getneighbouropt(Location local, int direction, Location & out) {
		if (direction >= 0) {
			out.row = local.row + DIRECTIONS[direction][0];
			out.col = local.col + DIRECTIONS[direction][1];
		} else {
			out = local;
		}
		return inside(out);
	}
	static uint manhattandist(const Location & a, const Location & b) {
		return abs(a.row - b.row) + abs(a.col - b.col);
	}

	static uint eucliddist2(const Location & a, const Location & b) {
		int dr = a.row - b.row;
		int dc = a.col - b.col;
		return (dr * dr) + (dc * dc);
	}

	static uint infinitydist(const Location & a, const Location & b) {
		int dr = a.row - b.row;
		int dc = a.col - b.col;
		return max(abs(dr), abs(dc));
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

static void add_enemy_kernel(Submap & sm, const Location & center, int delta)
{
	Location kernel;
	for (kernel.row = 0; kernel.row < AttackKernelSize; ++kernel.row) {
		for (kernel.col = 0; kernel.col < AttackKernelSize; ++kernel.col) {
			if (!AttackKernel[kernel.row][kernel.col])
				continue;

			Location other
				((center.row + kernel.row - AttackKernelRadius),
				 (center.col + kernel.col - AttackKernelRadius));
			if (!sm.inside(other))
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

		uint8_t killed : 1;
		uint8_t killer : 1;
		uint8_t sentoffense : 1;

		AntMove(const Location & p, int8_t dir) :
			pos(p), direction(dir),
			collided(0), collider(0), killed(0), sentoffense(0)
		{}
	};

	vector<AntMove> antmoves;
	uint hash;
	uint nrcollided;

	static const uint8_t AntsMask = 0xe0;
	static const uint8_t AntsShift = 5;
	static const uint8_t AttackMask = 0x1f;
	Submap map;

	int worstvalue;
	int worstopposingidx;

	PlayerMove() {
		reset();
	}

	void reset() {
		hash = 0;
		nrcollided = 0;
		worstvalue = numeric_limits<int>::max();
		worstopposingidx = -1;
		antmoves.clear();
	}

	void computehash() {
		hash = 0;
		for (uint idx = 0; idx < antmoves.size(); ++idx)
			hash = (hash * 5) + antmoves[idx].direction;
	}

	void ant_mark(uint idx) {
		// assume: direction and pos already set
		AntMove & am = antmoves[idx];

		if (!map.inside(am.pos))
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
				if (other == idx || am.pos != otherm.pos)
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

		if (!map.inside(am.pos))
			return;

		uint prevnrants = map[am.pos] >> AntsShift;
		assert(prevnrants >= 1);

		map[am.pos] -= 1 << AntsShift;
		if (prevnrants == 1) {
			add_enemy_kernel(map, am.pos, -1);
		} else if (prevnrants == 2) {
			for (uint other = 0; other < antmoves.size(); ++other) {
				AntMove & otherm = antmoves[other];
				if (other == idx || am.pos != otherm.pos)
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
					if (other == idx || am.pos != otherm.pos)
						continue;

					otherm.collider = 0;
					break;
				}
			}
			nrcollided--;
		}
	}

	void reset_killed()
	{
		for (uint idx = 0; idx < antmoves.size(); ++idx) {
			PlayerMove::AntMove & am = antmoves[idx];
			am.sentoffense = 0;
			am.killed = 0;
			am.killer = 0;
		}
	}
};

struct Tactical::Theater {
	bool needupdate;
	bool ismyperspective;
	Location offset;
	Submap basesm;
	vector<Location> myants;
	vector<Location> enemyants;
	BaseSubmap<unsigned int> foodpotential;

	vector<PlayerMove *> mymoves;
	vector<PlayerMove *> enemymoves;
	uint myevaluated;
	uint enemyevaluated;

	Theater(Data & d) {reset(d);}

	void reset(Data & d);
	void flipsides();
	uint mybestmove();
	uint enemybestmove();
};

struct Tactical::ShadowAnt {
	Location pos;
	bool hastactical;
};

struct Tactical::Data {
	vector<ShadowAnt> myshadowants;
	vector<Theater *> theaters;
	uint nrmoves;

	vector<PlayerMove *> unusedmoves;
	vector<Theater *> unusedtheaters;

	Data() {reset();}
	~Data() {
		reset();

		while (!unusedtheaters.empty()) {
			delete unusedtheaters.back();
			unusedtheaters.pop_back();
		}
		while (!unusedmoves.empty()) {
			delete unusedmoves.back();
			unusedmoves.pop_back();
		}
	}

	void reset() {
		nrmoves = 0;
		myshadowants.clear();

		while (!theaters.empty()) {
			Theater * theater = theaters.back();
			theaters.pop_back();
			theater->reset(*this);
			unusedtheaters.push_back(theater);
		}
	}

	PlayerMove * allocmove() {
		if (!unusedmoves.empty()) {
			PlayerMove * pm = unusedmoves.back();
			unusedmoves.pop_back();
			pm->reset();
			return pm;
		}
		return new PlayerMove;
	}

	Theater * alloctheater() {
		if (!unusedtheaters.empty()) {
			Theater * th = unusedtheaters.back();
			unusedtheaters.pop_back();
			return th;
		}
		return new Theater(*this);
	}
};

void Tactical::Theater::reset(Data & d)
{
	needupdate = true;
	ismyperspective = true;
	myants.clear();
	enemyants.clear();

	while (!mymoves.empty()) {
		d.unusedmoves.push_back(mymoves.back());
		mymoves.pop_back();
	}
	while (!enemymoves.empty()) {
		d.unusedmoves.push_back(enemymoves.back());
		enemymoves.pop_back();
	}
	myevaluated = 0;
	enemyevaluated = 0;
}

void Tactical::Theater::flipsides()
{
	ismyperspective = !ismyperspective;
	basesm.flipsides();
	swap(myants, enemyants);
	swap(mymoves, enemymoves);
	swap(myevaluated, enemyevaluated);
}

uint Tactical::Theater::mybestmove()
{
	uint bestidx = 0;
	int bestvalue = numeric_limits<int>::min();

	for (uint idx = 0; idx < mymoves.size(); ++idx) {
		int value = mymoves[idx]->worstvalue;
		if (value > bestvalue) {
			bestvalue = value;
			bestidx = idx;
		}
	}

	return bestidx;
}

uint Tactical::Theater::enemybestmove()
{
	uint bestidx = 0;
	int bestvalue = numeric_limits<int>::min();

	for (uint idx = 0; idx < enemymoves.size(); ++idx) {
		int value = enemymoves[idx]->worstvalue;
		if (value > bestvalue) {
			bestvalue = value;
			bestidx = idx;
		}
	}

	return bestidx;
}

struct Outcome {
	Tactical::Theater & th;
	BaseSubmap<char> map;

	Outcome(Tactical::Theater & th, int myidx, int enemyidx) : th(th) {
		Location local;
		for (local.row = 0; local.row < Submap::Size; ++local.row) {
			for (local.col = 0; local.col < Submap::Size; ++local.col) {
				unsigned char field = th.basesm[local];

				if (field & Submap::Water) {
					map[local] = '%';
				} else if (field & Submap::Food) {
					map[local] = '*';
				} else if (field & Submap::Hill) {
					if (field & Submap::Mine)
						map[local] = 'H';
					else
						map[local] = 'h';
				} else {
					uint8_t myfield = th.mymoves[myidx]->map[local];
					uint8_t enemyfield = th.enemymoves[enemyidx]->map[local];

					if (myfield & PlayerMove::AntsMask) {
						map[local] = 'Q' + (myfield >> PlayerMove::AntsShift);
					} else if (enemyfield & PlayerMove::AntsMask) {
						map[local] = 'q' + (enemyfield >> PlayerMove::AntsShift);
					} else {
						if (myfield & PlayerMove::AttackMask) {
							if (enemyfield & PlayerMove::AttackMask)
								map[local] = '+';
							else
								map[local] = '<';
						} else if (enemyfield & PlayerMove::AttackMask) {
							map[local] = '>';
						} else {
							map[local] = '.';
						}
					}
				}
			}
		}

		markmove(*th.mymoves[myidx], "ADC");
		markmove(*th.enemymoves[enemyidx], "adc");
	}

	void markmove(PlayerMove & move, const char * markers) {
		for (uint antidx = 0; antidx < move.antmoves.size(); ++antidx) {
			PlayerMove::AntMove & am = move.antmoves[antidx];
			if (!th.basesm.inside(am.pos))
				continue;

			if (am.collided)
				map[am.pos] = markers[2];
			else if (am.killed)
				map[am.pos] = markers[1];
			else
				map[am.pos] = markers[0];
		}
	}
};

ostream & operator<<(ostream & out, const Outcome & outcome) {
	Location local;
	for (local.row = 0; local.row < Submap::Size; ++local.row) {
		for (local.col = 0; local.col < Submap::Size; ++local.col) {
			out << outcome.map[local];
		}
		out << endl;
	}
	return out;
}

Tactical::Tactical(Bot & bot_) :
	bot(bot_),
	state(bot_.state),
	d(*new Data)
{
}

Tactical::~Tactical()
{
	delete &d;
}

void Tactical::init()
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

struct CloserToCenterCompare {
	bool operator()(const Location & a, const Location & b) const {
		uint a2 = Submap::eucliddist2(a, Location(Submap::Radius, Submap::Radius));
		uint b2 = Submap::eucliddist2(b, Location(Submap::Radius, Submap::Radius));
		return a2 < b2;
	}
};

struct NewMove {
	NewMove(Tactical & t_, Tactical::Theater & th_, uint myidxbase, const char * why) :
		t(t_), d(t_.d), th(th_), state(t_.state), shouldcommit(false)
	{
		state.bug << "Preliminary insert new move " << th.mymoves.size() << " due to " << why
			<< " perspective " << (th.ismyperspective ? "mine" : "enemy") << endl;
		th.mymoves.push_back(d.allocmove());
		*th.mymoves.back() = *th.mymoves[myidxbase];

		PlayerMove & pm = *th.mymoves.back();
		pm.worstvalue = numeric_limits<int>::max();
		pm.worstopposingidx = -1;
	}

	~NewMove()
	{
		if (!shouldcommit) {
			d.unusedmoves.push_back(th.mymoves.back());
			th.mymoves.pop_back();
			state.bug << "Abort new move " << th.mymoves.size() << endl;
		} else {
			state.bug << "Commit new move " << (th.mymoves.size() - 1) << endl;
			th.needupdate = true;
			d.nrmoves++;
		}
	}

	PlayerMove & move() {return *th.mymoves.back();}

	void antmove(uint myantidx, int direction)
	{
		PlayerMove & pm = move();
		pm.ant_unmark(myantidx);

		PlayerMove::AntMove & am = pm.antmoves[myantidx];
		am.direction = direction;
		th.basesm.getneighbouropt(th.myants[myantidx], direction, am.pos);

		pm.ant_mark(myantidx);

		state.bug << "antmove " << myantidx << " to " << cdir(direction) << " aka " << am.pos << endl;
	}

	void commit()
	{
		assert(!shouldcommit);

		PlayerMove & pm = move();
		pm.computehash();

		for (uint idx = 0; idx < th.mymoves.size() - 1; ++idx) {
			PlayerMove & other = *th.mymoves[idx];
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
	Tactical::Theater & th;
	State & state;
	bool shouldcommit;
};

struct Improve {
	Tactical & t;
	Tactical::Data & d;
	Tactical::Theater & th;
	State & state;
	uint myidx;
	uint enemyidx;

	struct Death {
		uint myantidx;
		uint enemyantidx;
	};

	vector<Death> deaths;

	Improve(Tactical & t_, Tactical::Theater & th_, uint myidx_, uint enemyidx_) :
		t(t_), d(t_.d), th(th_), state(t_.state),
		myidx(myidx_), enemyidx(enemyidx_)
	{
	}

	void compute_deaths()
	{
		PlayerMove & mymove = themymove();
		PlayerMove & enemymove = theenemymove();

		mymove.reset_killed();
		enemymove.reset_killed();

		for (uint myantidx = 0; myantidx < mymove.antmoves.size(); ++myantidx) {
			PlayerMove::AntMove & am = mymove.antmoves[myantidx];
			if (am.collided || !th.basesm.inside(am.pos))
				continue;

			uint nrenemies = enemymove.map[am.pos] & PlayerMove::AttackMask;
			if (nrenemies == 0)
				continue;

			for (uint enemyantidx = 0; enemyantidx < enemymove.antmoves.size(); ++enemyantidx) {
				PlayerMove::AntMove & em = enemymove.antmoves[enemyantidx];
				if (em.collided || !th.basesm.inside(em.pos))
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
				em.killer = 1;

				// ignore deaths that are too far from the center, because we cannot assess them accurately
				if (th.basesm.eucliddist2(am.pos, Location(Submap::Radius, Submap::Radius)) > MaxAttackRadius2)
					break;

				deaths.push_back(Death());
				deaths.back().myantidx = myantidx;
				deaths.back().enemyantidx = enemyantidx;
				break;
			}
		}

		for (uint enemyantidx = 0; enemyantidx < enemymove.antmoves.size(); ++enemyantidx) {
			PlayerMove::AntMove & am = enemymove.antmoves[enemyantidx];
			if (am.collided || !th.basesm.inside(am.pos))
				continue;

			uint nrenemies = mymove.map[am.pos] & PlayerMove::AttackMask;
			if (nrenemies == 0)
				continue;

			for (uint myantidx = 0; myantidx < mymove.antmoves.size(); ++myantidx) {
				PlayerMove::AntMove & em = mymove.antmoves[myantidx];
				if (em.collided || !th.basesm.inside(em.pos))
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
				em.killer = 1;
				break;
			}
		}
	}

	void avoid_collisions()
	{
		NewMove nm(t, th, myidx, "collision avoidance");
		const PlayerMove & mymove = themymove();
		const PlayerMove & enemymove = theenemymove();

		uint p = getprime();
		uint ofs = fastrng() % mymove.antmoves.size();
		for (uint preidx = 0; preidx < mymove.antmoves.size(); ++preidx) {
			uint myantidx = (p * preidx + ofs) % mymove.antmoves.size();
			const PlayerMove::AntMove & am = mymove.antmoves[myantidx];
			if (!am.collider)
				continue;

#define checkacceptable(label) \
	do { \
		Location to; \
		if (!th.basesm.getneighbouropt(from, dir, to)) break; \
		if (th.basesm[to] & Submap::Water) break; \
		if (nm.move().map[to] & PlayerMove::AntsMask) break; \
		altdir = dir; \
		if ((enemymove.map[to] & PlayerMove::AttackMask) == 0) \
			goto label; \
	} while(0)

			Location from = th.myants[myantidx];
			int olddir = am.direction;
			if (olddir == -1) {
				int altdir = -1;
				const int * dirperm = getdirperm();
				for (int predir = 0; predir < TDIRECTIONS; ++predir) {
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

	bool pushmove(NewMove & nm, uint myantidx, int direction)
	{
		state.bug << "pushmove " << myantidx << " at " << nm.move().antmoves[myantidx].pos << " to " << cdir(direction) << endl;

		// This has a bias towards retreating away from the enemy
		for (int iters = 0; iters < 5; ++iters) {
			PlayerMove::AntMove & am = nm.move().antmoves[myantidx];
			nm.antmove(myantidx, direction);

			if (!am.collided) {
				state.bug << "pushmove success in iter " << iters << endl;
				return true;
			}

			uint p = getprime();
			uint ofs = fastrng() % th.myants.size();
			for (uint preidx = 0; preidx < th.myants.size(); ++preidx) {
				uint otheridx = (p * preidx + ofs) % th.myants.size();
				PlayerMove::AntMove & otherm = nm.move().antmoves[otheridx];

				if (otheridx == myantidx || otherm.pos != am.pos)
					continue;

				int allfreedirection = -2;
				int collisiondirection = -2;
				int attackeddirection = -2;

#define check_candidate(dir) \
	do { \
		Location n; \
		if (!th.basesm.getneighbouropt(otherm.pos, (dir), n)) break; \
		if (th.basesm[n] & Submap::Water) break; \
		bool collision = nm.move().map[n] & PlayerMove::AntsMask; \
		bool attacked = theenemymove().map[n] & PlayerMove::AttackMask; \
		if (!collision && !attacked) { \
			allfreedirection = dir; \
			goto found; \
		} else if (!attacked) { \
			if (collisiondirection <= -2) collisiondirection = dir; \
		} else if (!collision) { \
			if (attackeddirection <= -2) attackeddirection = dir; \
		} \
	} while(0)

				int altdir = (direction + 1 + (fastrng() & 2)) % TDIRECTIONS;

				check_candidate(direction);

				if (altdir != otherm.direction)
					check_candidate(altdir);

				altdir = (altdir + 2) % TDIRECTIONS;
				if (altdir != otherm.direction)
					check_candidate(altdir);

				if (-1 != otherm.direction)
					check_candidate(-1);

			found:
				if (allfreedirection >= -1)
					direction = allfreedirection;
				else if (collisiondirection >= -1)
					direction = collisiondirection;
				else if (attackeddirection >= -1)
					direction = attackeddirection;

				myantidx = otheridx;
				break;
			}
		}

		state.bug << "pushmove failed" << endl;
		return false;
	}

	void compute_retreat_dir(const Location & orig, int origdir, const Location & attacker, int & noattackdir, int & awaydir)
	{
		const PlayerMove & mymove = themymove();
		const PlayerMove & enemymove = theenemymove();

		bool attackdirblocked = true;
		bool awaydirblocked = true;
		noattackdir = origdir;
		awaydir = origdir;

		const int * optdirperm = getoptdirperm();
		for (int predir = 0; predir < TDIRECTIONS+1; ++predir) {
			int dir = optdirperm[predir];
			if (dir == origdir)
				continue;

			Location n;
			if (!th.basesm.getneighbouropt(orig, dir, n))
				continue;
			if (th.basesm[n] & Submap::Water)
				continue;

			if ((enemymove.map[n] & PlayerMove::AttackMask) == 0) {
				if (attackdirblocked) {
					noattackdir = dir;
					if (!(mymove.map[n] & PlayerMove::AntsMask))
						attackdirblocked = false;
				}
			} else {
				Location kernel
					(attacker.row - n.row + AttackKernelRadius,
					 attacker.col - n.col + AttackKernelRadius);
				if
					(kernel.row < 0 || kernel.row >= AttackKernelSize ||
					 kernel.col < 0 || kernel.col >= AttackKernelSize ||
					 !AttackKernel[kernel.row][kernel.col])
				{
					if (awaydirblocked) {
						awaydir = dir;
						if (!(mymove.map[n] & PlayerMove::AntsMask))
							awaydirblocked = false;
					}
				}
			}
		}
	}

	void do_retreat()
	{
		for (uint idx = 0; idx < deaths.size(); ++idx) {
			const PlayerMove & mymove = themymove();
			const PlayerMove & enemymove = theenemymove();
			Death & death = deaths[idx];
			const Location & orig = th.myants[death.myantidx];
			const Location & other = enemymove.antmoves[death.enemyantidx].pos;

			int origdir = mymove.antmoves[death.myantidx].direction;
			int alt1dir;
			int alt2dir;

			compute_retreat_dir(orig, origdir, other, alt1dir, alt2dir);

			if (alt1dir != origdir) {
				NewMove nm(t, th, myidx, "retreat: unattacked field");
				if (pushmove(nm, death.myantidx, alt1dir))
					nm.commit();
			}
			if (alt2dir != origdir && alt2dir != alt1dir) {
				NewMove nm(t, th, myidx, "retreat: get away from killer");
				if (pushmove(nm, death.myantidx, alt2dir))
					nm.commit();
			}
		}

		if (deaths.size() >= 2) {
			NewMove nm(t, th, myidx, "retreat: all at once");
			const PlayerMove & mymove = themymove();
			const PlayerMove & enemymove = theenemymove();

			uint p = getprime();
			uint ofs = fastrng() % deaths.size();
			for (uint preidx = 0; preidx < deaths.size(); ++preidx) {
				uint idx = (p * preidx + ofs) % deaths.size();
				Death & death = deaths[idx];
				const Location & orig = th.myants[death.myantidx];
				const Location & other = enemymove.antmoves[death.enemyantidx].pos;

				int origdir = mymove.antmoves[death.myantidx].direction;
				int alt1dir;
				int alt2dir;

				compute_retreat_dir(orig, origdir, other, alt1dir, alt2dir);

				bool moved = false;
				if (alt1dir != origdir) {
					if (pushmove(nm, death.myantidx, alt1dir))
						moved = true;
				}
				if (!moved && alt2dir != origdir) {
					pushmove(nm, death.myantidx, alt2dir);
				}
			}

			nm.commit();
		}
	}

	void send_to_kill(uint enemyantidx, const char * why)
	{
		PlayerMove & mymove = themymove();
		PlayerMove & enemymove = theenemymove();
		PlayerMove::AntMove & enemyant = enemymove.antmoves[enemyantidx];
		Location enemypos = enemyant.pos;

		enemymove.antmoves[enemyantidx].sentoffense = 1;

		if (!th.basesm.inside(enemypos))
			return;
		if (th.myants.empty())
			return;

		state.bug << "check whether we can send ants to kill " << enemyantidx << " at " << enemypos << endl;

		uint additionalidx[NrAttackNeighboursUpperBound];
		int additionaldir[NrAttackNeighboursUpperBound];
		uint additionalvalue[NrAttackNeighboursUpperBound];
		uint nradditional = 0;

		uint nrattackers = mymove.map[enemypos] & PlayerMove::AttackMask;
		uint bestattacker = numeric_limits<uint>::max();

		uint p = getprime();
		uint ofs = fastrng() % th.myants.size();
		for (uint preidx = 0; preidx < th.myants.size(); ++preidx) {
			uint myantidx = (preidx * p + ofs) % th.myants.size();
			const Location & orig = th.myants[myantidx];
			const Location & current = mymove.antmoves[myantidx].pos;

			state.bug << "  consider ant " << myantidx << " at " << current << " (from " << orig << ")" << endl;

			if (th.basesm.eucliddist2(current, enemypos) <= AttackRadius2) {
				state.bug << "    already attacking" << endl;
				// already attacking the target
				bestattacker = min(bestattacker, uint(enemymove.map[current] & PlayerMove::AttackMask));
				continue;
			}

			if (th.basesm.eucliddist2(orig, enemypos) > AttackNeighboursRadius2)
				continue; // too far away

			state.bug << "  close enough" << endl;

			int bestdirection = -1;
			uint bestvalue = numeric_limits<uint>::max();
			const int * optdirperm = getoptdirperm();
			for (int predir = 0; predir < TDIRECTIONS + 1; ++predir) {
				int dir = optdirperm[predir];
				Location n;
				if (!th.basesm.getneighbouropt(orig, dir, n))
					continue;
				if (th.basesm.eucliddist2(n, enemypos) > AttackRadius2)
					continue;
				if (th.basesm[n] & Submap::Water)
					continue;

				state.bug << "  consider moving to " << n << " " << cdir(dir) << endl;

				uint value = enemymove.map[n] & PlayerMove::AttackMask;
				if (value < bestvalue) {
					bestvalue = value;
					bestdirection = dir;
				}
			}

			state.bug << "  best choice: " << cdir(bestdirection) << " attackers " << bestvalue << endl;

			additionalidx[nradditional] = myantidx;
			additionaldir[nradditional] = bestdirection;
			additionalvalue[nradditional] = bestvalue;
			nradditional++;
		}

		if (!nradditional)
			return;

		NewMove nm(t, th, myidx, why);
		for (uint idx = 0; idx < nradditional; ++idx) {
			state.bug << idx << " nrattackers vs. bestattacker " << nrattackers << " " << bestattacker << endl;
			if (nrattackers > bestattacker)
				break;

			uint myantidx = additionalidx[idx];
			const Location & orig = th.myants[myantidx];
			int dir = additionaldir[idx];
			Location n;
			if (!th.basesm.getneighbouropt(orig, dir, n)) {
				abort();
			}

			if (nm.move().map[n] & PlayerMove::AntsMask)
				continue;

			state.bug << "moving " << myantidx << " to " << n << endl;

			nm.antmove(myantidx, dir);
			bestattacker = min(bestattacker, additionalvalue[idx]);
			nrattackers++;
		}

		state.bug << "final nrattackers " << nrattackers << " vs " << bestattacker << endl;

		if (nrattackers >= bestattacker) {
			nm.commit();
		}
	}

	void do_offense()
	{
		for (uint enemyantidx = 0; enemyantidx < th.enemyants.size(); ++enemyantidx) {
			const PlayerMove & enemymove = theenemymove();
			const PlayerMove::AntMove & enemyant = enemymove.antmoves[enemyantidx];

			if (enemyant.collided || enemyant.sentoffense || (enemyant.killed && !enemyant.killer))
				continue;

			if (th.basesm.eucliddist2(enemyant.pos, Location(Submap::Radius, Submap::Radius)) > MaxAttackRadius2)
				continue;

			send_to_kill(enemyantidx, "offense");
		}
	}

	void do_improve()
	{
		state.bug << "improve " << myidx << " vs " << enemyidx
			<< " perspective: " << (th.ismyperspective ? "mine" : "enemy") << endl;

		//
		compute_deaths();

		state.bug << Outcome(th, myidx, enemyidx);

		//
		if (themymove().nrcollided) {
			state.bug << "avoid collisions" << endl;
			avoid_collisions();
		}

		if (!deaths.empty()) {
			// if one of ours is killed, see if we can retreat
			do_retreat();
		}

		// try an unprovoked offense
		do_offense();
	}

	PlayerMove & themymove() const {return *th.mymoves[myidx];}
	PlayerMove & theenemymove() const {return *th.enemymoves[enemyidx];}
};

/**
 * Look at the outcome of myidx move vs enemyidx move, and try to generate some better
 * alternative moves for myself.
 */
void Tactical::improve(uint theateridx, uint myidx, uint enemyidx)
{
	Improve imp(*this, *d.theaters[theateridx], myidx, enemyidx);
	imp.do_improve();
}

int Tactical::evaluate_ant_positions(Theater & th, PlayerMove & mymove, bool myperspective)
{
	int value = 0;
	uint8_t enemyflag = myperspective ? Submap::Enemy : Submap::Mine;

	for (uint antidx = 0; antidx < mymove.antmoves.size(); ++antidx) {
		PlayerMove::AntMove & am = mymove.antmoves[antidx];
		if (am.collided || am.killed)
			continue;

		Location global
			((am.pos.row + th.offset.row + state.rows) % state.rows,
			 (am.pos.col + th.offset.col + state.cols) % state.cols);

		if (!myperspective) {
			value += (MaxZocValue - min(bot.m_zoc.m_me[global], uint(MaxZocValue))) * ZocFactor;
		} else {
			value += (MaxZocValue - min(bot.m_zoc.m_enemy[global], uint(MaxZocValue))) * ZocFactor;
		}

		if (th.basesm.inside(am.pos)) {
			value += th.foodpotential[am.pos] * FoodFactor;

			if ((th.basesm[am.pos] & (Submap::Hill | enemyflag)) == (Submap::Hill | enemyflag))
				value += HillValue;
		}

		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n;
			if (!th.basesm.getneighbour(am.pos, dir, n))
				continue;

			if ((th.basesm[n] & (Submap::Hill | enemyflag)) == (Submap::Hill | enemyflag))
				value += HillValue;
		}
	}

	return value;
}


void Tactical::compute_deaths(Theater & th, PlayerMove & mymove, PlayerMove & enemymove, int & mydeaths, int & enemydeaths)
{
	for (uint myantidx = 0; myantidx < mymove.antmoves.size(); ++myantidx) {
		PlayerMove::AntMove & am = mymove.antmoves[myantidx];
		if (am.collided || !th.basesm.inside(am.pos))
			continue;

		am.killed = 0;

		uint nrenemies = enemymove.map[am.pos] & PlayerMove::AttackMask;
		if (nrenemies == 0)
			continue;

		for (uint enemyantidx = 0; enemyantidx < enemymove.antmoves.size(); ++enemyantidx) {
			const PlayerMove::AntMove & em = enemymove.antmoves[enemyantidx];
			if (em.collided || !th.basesm.inside(em.pos))
				continue;

			if (th.basesm.eucliddist2(em.pos, am.pos) > AttackRadius2)
				continue;
			if ((mymove.map[em.pos] & PlayerMove::AttackMask) > nrenemies)
				continue;

			am.killed = 1;
			mydeaths++;
			break;
		}
	}

	for (uint enemyantidx = 0; enemyantidx < enemymove.antmoves.size(); ++enemyantidx) {
		PlayerMove::AntMove & am = enemymove.antmoves[enemyantidx];
		if (am.collided || !th.basesm.inside(am.pos))
			continue;

		am.killed = 0;

		uint nrenemies = mymove.map[am.pos] & PlayerMove::AttackMask;
		if (nrenemies == 0)
			continue;

		for (uint myantidx = 0; myantidx < mymove.antmoves.size(); ++myantidx) {
			const PlayerMove::AntMove & em = mymove.antmoves[myantidx];
			if (em.collided || !th.basesm.inside(em.pos))
				continue;

			if (th.basesm.eucliddist2(am.pos, em.pos) > AttackRadius2)
				continue;
			if ((enemymove.map[em.pos] & PlayerMove::AttackMask) > nrenemies)
				continue;

			am.killed = 1;
			enemydeaths++;
			break;
		}
	}
}

void Tactical::evaluate_new_moves(uint theateridx)
{
	Theater & th = *d.theaters[theateridx];

	if (th.myevaluated >= th.mymoves.size() && th.enemyevaluated >= th.enemymoves.size())
		return;

	for (uint myidx = 0; myidx < th.mymoves.size(); ++myidx) {
		PlayerMove & mymove = *th.mymoves[myidx];

		for
			(uint enemyidx = (myidx < th.myevaluated) ? th.enemyevaluated : 0;
			 enemyidx < th.enemymoves.size(); ++enemyidx)
		{
			PlayerMove & enemymove = *th.enemymoves[enemyidx];
			int myvalue = mymove.nrcollided * ValueLoss;
			int enemyvalue = enemymove.nrcollided * ValueLoss;

			int mydeaths = 0;
			int enemydeaths = 0;
			compute_deaths(th, mymove, enemymove, mydeaths, enemydeaths);

			myvalue += mydeaths * ValueLoss + enemydeaths * ValueKill;
// 			if ((mydeaths || mymove.nrcollided) && !enemydeaths)
// 				myvalue += ValueUselessDeath;

			enemyvalue += mydeaths * ValueKill + enemydeaths * ValueLoss;
// 			if ((enemydeaths || enemymove.nrcollided) && !mydeaths)
// 				enemyvalue += ValueUselessDeath;

			int myposvalue = evaluate_ant_positions(th, mymove, true);
			int enemyposvalue = evaluate_ant_positions(th, enemymove, false);

			myvalue += myposvalue - enemyposvalue;
			enemyvalue += enemyposvalue - myposvalue;

			state.bug << " " << myidx << " vs. " << enemyidx
				<< "  value " << myvalue << " vs " << enemyvalue << endl;
			state.bug << Outcome(th, myidx, enemyidx);

			if (myvalue < mymove.worstvalue) {
				mymove.worstvalue = myvalue;
				mymove.worstopposingidx = enemyidx;
			}
			if (enemyvalue < enemymove.worstvalue) {
				enemymove.worstvalue = enemyvalue;
				enemymove.worstopposingidx = myidx;
			}
		}
	}

	th.myevaluated = th.mymoves.size();
	th.enemyevaluated = th.enemymoves.size();
}

void Tactical::generate_theater(const Location & center)
{
	uint theateridx = d.theaters.size();
	Theater & th = *d.alloctheater();
	d.theaters.push_back(&th);

	state.bug << "Tactical theater " << theateridx << " around " << center << endl;

	gensubmap(th.basesm, center);
	compute_foodpotential(th.basesm, th.foodpotential);
	th.offset = Location
		((center.row - Submap::Radius + state.rows) % state.rows,
		 (center.col - Submap::Radius + state.cols) % state.cols);

	// generate init moves
	th.mymoves.push_back(d.allocmove());
	th.enemymoves.push_back(d.allocmove());
	d.nrmoves += 2;

	th.mymoves[0]->map.fill(0);
	th.enemymoves[0]->map.fill(0);

	// collect all ants and sort by distance to center
	Location local;
	for (local.row = 0; local.row < Submap::Size; ++local.row) {
		for (local.col = 0; local.col < Submap::Size; ++local.col) {
			unsigned char field = th.basesm[local];

			if (!(field & Submap::Ant))
				continue;

			if (field & Submap::Enemy) {
				th.enemyants.push_back(local);
			} else {
				th.myants.push_back(local);

				if (th.basesm.infinitydist(local, Location(Submap::Radius, Submap::Radius)) <= TheaterCoreInfinity) {
					uint antidx = bot.myantidx_at(state.addLocations(th.offset, local));
					d.myshadowants[antidx].hastactical = true;
				}
			}
		}
	}

	sort(th.myants.begin(), th.myants.end(), CloserToCenterCompare());
	sort(th.enemyants.begin(), th.enemyants.end(), CloserToCenterCompare());

	// Fill in default moves
	th.mymoves[0]->antmoves.reserve(th.myants.size());
	for (uint idx = 0; idx < th.myants.size(); ++idx) {
		const Location & local = th.myants[idx];
		Location global
			((local.row + th.offset.row) % state.rows,
			 (local.col + th.offset.col) % state.cols);
		Ant & ant = bot.m_ants[bot.myantidx_at(global)];
		Location nextlocal;
		th.basesm.getneighbouropt(local, ant.direction, nextlocal);
		th.mymoves[0]->antmoves.push_back(PlayerMove::AntMove(nextlocal, ant.direction));
		th.mymoves[0]->ant_mark(idx);
	}

	th.enemymoves[0]->antmoves.reserve(th.enemyants.size());
	for (uint idx = 0; idx < th.enemyants.size(); ++idx) {
		const Location & local = th.enemyants[idx];
		th.enemymoves[0]->antmoves.push_back(PlayerMove::AntMove(local, -1));
		th.enemymoves[0]->ant_mark(idx);
	}

	// Finalize the moves
	th.mymoves[0]->computehash();
	th.enemymoves[0]->computehash();

	state.bug << Outcome(th, 0, 0);
}

void Tactical::run_theater(uint theateridx)
{
	state.bug << "Run tactical theater " << theateridx << " (turn " << state.turn << ")" << endl;

	Theater & th = *d.theaters[theateridx];
	th.needupdate = false;

	evaluate_new_moves(theateridx);

	uint mybestidx = th.mybestmove();
	uint mybestenemyidx = th.mymoves[mybestidx]->worstopposingidx;

	uint enemybestidx = th.enemybestmove();
	uint enemybestmyidx = th.enemymoves[enemybestidx]->worstopposingidx;

	improve(theateridx, mybestidx, mybestenemyidx);
	th.flipsides();
	improve(theateridx, mybestenemyidx, mybestidx);
	if (enemybestidx != mybestenemyidx || enemybestmyidx != mybestidx)
		improve(theateridx, enemybestidx, enemybestmyidx);
	th.flipsides();

	evaluate_new_moves(theateridx);

	state.bug << "-----------------------" << endl;
}

void Tactical::make_moves()
{
	for (uint theateridx = 0; theateridx < d.theaters.size(); ++theateridx) {
		Theater & th = *d.theaters[theateridx];
		uint bestmoveidx = th.mybestmove();
		const PlayerMove & move = *th.mymoves[bestmoveidx];

		state.bug << "Best move in theater " << theateridx << ": " << bestmoveidx << endl;

		for (uint idx = 0; idx < move.antmoves.size(); ++idx) {
			int dir = move.antmoves[idx].direction;
			Location antpos = th.myants[idx];
			Location global
				((antpos.row + th.offset.row) % state.rows,
				(antpos.col + th.offset.col) % state.cols);
			uint antidx = bot.myantidx_at(global);
			Ant & ant = bot.m_ants[antidx];

			ant.direction = dir;

//			state.bug << "  Ant at " << ant.where
//				<< " tactical move to " << state.getLocation(ant.where, dir) << " (" << cdir(ant.direction) << ")" << endl;
		}
	}
}

bool Tactical::timeover()
{
	if (d.nrmoves > MaxEval) {
		state.bug << "Max number of moves exceeded" << endl;
		return true;
	}

	return false;
}

void Tactical::run()
{
	d.reset();

	d.myshadowants.resize(bot.m_ants.size());
	for (uint antidx = 0; antidx < bot.m_ants.size(); ++antidx) {
		ShadowAnt & sa = d.myshadowants[antidx];
		sa.pos = bot.m_ants[antidx].where;
		sa.hastactical = false;
	}

	// generate theaters
	uint p = getprime();
	uint ofs = bot.m_ants.size();
	for (uint preidx = 0; preidx < bot.m_ants.size(); ++preidx) {
		uint antidx = (p * preidx + ofs) % bot.m_ants.size();
		Ant & ant = bot.m_ants[antidx];
		if (d.myshadowants[antidx].hastactical || bot.m_zoc.m_enemy[ant.where] > AttackableManhattanDistance)
			continue;

		generate_theater(ant.where);
	}

	if (d.theaters.empty())
		return;

	state.bug << "Tactical turn " << state.turn << " on " << d.theaters.size() << " theaters" << endl;

	//
	uint theateridx = d.theaters.size() - 1;
	while (!timeover()) {
		theateridx = (theateridx + 1) % d.theaters.size();
		uint starttheater = theateridx;
		do {
			if (d.theaters[theateridx]->needupdate)
				break;
			theateridx = (theateridx + 1) % d.theaters.size();
		} while (theateridx != starttheater);

		if (!d.theaters[theateridx]->needupdate)
			break;

		run_theater(theateridx);
	}

	state.bug << "Total number of generated moves: " << d.nrmoves << endl;

	//
	make_moves();
}
