#include "tactical_sms.h"

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
static const uint MaxEval = 100000; ///< max number of evaluations

static const uint AttackableManhattanDistance = 5; ///< max distance for an ant to trigger theater creation
static const uint TheaterCoreInfinity = 2; ///< max infinity norm distance for an ant to "belong" to the theater
static const uint MaxAttackRadius2 = 34; ///< radius up to which we consider attacking an enemy (this is the maximum distance a core ant can possibly attack)

static float LearnGamma = 1.5;
//@}

static const float EpsilonValue = 0.0000001;

static const uint AttackRadius2 = 5;
static const uint AttackNeighboursRadius2 = 10;
static const int NrAttackNeighboursUpperBound = 49;


struct TacticalSms::Theater {
	bool needupdate;
	bool ismyperspective;
	bool aggressive;
	Location offset;
	Submap basesm;
	vector<Location> myants;
	vector<Location> enemyants;

	vector<PlayerMove *> mymoves;
	vector<PlayerMove *> enemymoves;
	uint mymove;
	uint mymove_unchangedrounds;

	Theater(Data & d) {reset(d);}

	void reset(Data & d);
	void flipsides();
	bool is_duplicate_mymove(uint myidx);
};

struct TacticalSms::ShadowAnt {
	Location pos;
	bool hastactical;
};

struct TacticalSms::Data {
	vector<ShadowAnt> myshadowants;
	vector<Theater *> theaters;
	uint nrmoves;
	uint nrevals;

	vector<PlayerMove *> unusedmoves;
	vector<Theater *> unusedtheaters;

	vector<uint> tmp_candidates;
	vector<bool> tmp_mask;

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
		nrevals = 0;
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

void TacticalSms::Theater::reset(Data & d)
{
	needupdate = true;
	ismyperspective = true;
	aggressive = false;
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
	mymove = 0;
	mymove_unchangedrounds = 0;
}

void TacticalSms::Theater::flipsides()
{
	ismyperspective = !ismyperspective;
	basesm.flipsides();
	swap(myants, enemyants);
	swap(mymoves, enemymoves);
}

bool TacticalSms::Theater::is_duplicate_mymove(uint myidx)
{
	PlayerMove & pm = *mymoves[myidx];

	for (uint idx = 0; idx < mymoves.size(); ++idx) {
		if (idx == myidx)
			continue;

		PlayerMove & other = *mymoves[idx];
		if (other.hash != pm.hash)
			continue;

		for (uint antidx = 0; antidx < pm.antmoves.size(); ++antidx) {
			if (pm.antmoves[antidx].direction != other.antmoves[antidx].direction)
				goto notequal;
		}

		return true;
	notequal: ;
	}

	return false;
}


struct OutcomeSms {
	TacticalSms::Theater & th;
	BaseSubmap<char> map;

	OutcomeSms(TacticalSms::Theater & th, int myidx, int enemyidx) : th(th) {
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

ostream & operator<<(ostream & out, const OutcomeSms & outcome) {
	Location local;
	for (local.row = 0; local.row < Submap::Size; ++local.row) {
		for (local.col = 0; local.col < Submap::Size; ++local.col) {
			out << outcome.map[local];
		}
		out << endl;
	}
	return out;
}

TacticalSms::TacticalSms(Bot & bot_) :
	TacticalSmBase(bot_),
	d(*new Data)
{
}

TacticalSms::~TacticalSms()
{
	delete &d;
}

void TacticalSms::init()
{
	TacticalSmBase::init();

	state.bug.time << "Using TacticalSms" << endl;

	LearnGamma = bot.getargfloat("LearnGamma", LearnGamma);

	state.bug.time << "LearnGamma = " << LearnGamma << endl;
}

struct CloserToCenterCompare {
	bool operator()(const Location & a, const Location & b) const {
		uint a2 = Submap::eucliddist2(a, Location(Submap::Radius, Submap::Radius));
		uint b2 = Submap::eucliddist2(b, Location(Submap::Radius, Submap::Radius));
		return a2 < b2;
	}
};

struct NewMove {
	NewMove(TacticalSms & t_, TacticalSms::Theater & th_, uint myidxbase, const char * why) :
		t(t_), d(t_.d), th(th_), state(t_.state), shouldcommit(false)
	{
		state.bug << "Preliminary insert new move " << th.mymoves.size() << " due to " << why
			<< " perspective " << (th.ismyperspective ? "mine" : "enemy") << endl;
		th.mymoves.push_back(d.allocmove());
		th.mymoves.back()->init(*th.mymoves[myidxbase]);
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

		state.bug << "    antmove " << myantidx << " to " << cdir(direction) << " aka " << am.pos << endl;
	}

	void commit()
	{
		assert(!shouldcommit);

		PlayerMove & pm = move();
		pm.computehash();

		if (th.is_duplicate_mymove(th.mymoves.size() - 1)) {
			state.bug << "New move is equal to an old move" << endl;
			return;
		}

		shouldcommit = true;
	}

	TacticalSms & t;
	TacticalSms::Data & d;
	TacticalSms::Theater & th;
	State & state;
	bool shouldcommit;
};

struct Improve {
	TacticalSms & t;
	TacticalSms::Data & d;
	TacticalSms::Theater & th;
	State & state;
	uint myidx;
	uint enemyidx;

	struct Death {
		uint myantidx;
		uint enemyantidx;
	};

	vector<Death> deaths;

	Improve(TacticalSms & t_, TacticalSms::Theater & th_, uint myidx_, uint enemyidx_) :
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
				pushmove(nm, myantidx, altdir);
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
				pushmove(nm, myantidx, altdir);
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
		state.bug << "  pushmove " << myantidx << " at " << nm.move().antmoves[myantidx].pos << " to " << cdir(direction) << endl;

		// This has a bias towards retreating away from the enemy
		for (int iters = 0; iters < 5; ++iters) {
			PlayerMove::AntMove & am = nm.move().antmoves[myantidx];
			nm.antmove(myantidx, direction);

			if (!am.collided) {
				state.bug << "  pushmove success in iter " << iters << endl;
				return true;
			}

			uint p = getprime();
			uint ofs = fastrng() % th.myants.size();
			for (uint preidx = 0; preidx < th.myants.size(); ++preidx) {
				uint otheridx = (p * preidx + ofs) % th.myants.size();
				PlayerMove::AntMove & otherm = nm.move().antmoves[otheridx];

				if (otheridx == myantidx || otherm.pos != am.pos)
					continue;

				const Location & otherorig = th.myants[otheridx];
				int allfreedirection = -2;
				int collisiondirection = -2;
				int attackeddirection = -2;

#define check_candidate(dir) \
	do { \
		Location n; \
		if (!th.basesm.getneighbouropt(otherorig, (dir), n)) break; \
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

		state.bug << "  pushmove failed" << endl;
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

//			state.bug << "  consider ant " << myantidx << " at " << current << " (from " << orig << ")" << endl;

			if (th.basesm.eucliddist2(current, enemypos) <= AttackRadius2) {
//				state.bug << "    already attacking" << endl;
				// already attacking the target
				bestattacker = min(bestattacker, uint(enemymove.map[current] & PlayerMove::AttackMask));
				continue;
			}

			if (th.basesm.eucliddist2(orig, enemypos) > AttackNeighboursRadius2)
				continue; // too far away

//			state.bug << "  close enough" << endl;

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

//				state.bug << "  consider moving to " << n << " " << cdir(dir) << endl;

				uint value = enemymove.map[n] & PlayerMove::AttackMask;
				if (value < bestvalue) {
					bestvalue = value;
					bestdirection = dir;
				}
			}

//			state.bug << "  best choice: " << cdir(bestdirection) << " attackers " << bestvalue << endl;

			additionalidx[nradditional] = myantidx;
			additionaldir[nradditional] = bestdirection;
			additionalvalue[nradditional] = bestvalue;
			nradditional++;
		}

		if (!nradditional)
			return;

		NewMove nm(t, th, myidx, why);
		for (uint idx = 0; idx < nradditional; ++idx) {
//			state.bug << idx << " nrattackers vs. bestattacker " << nrattackers << " " << bestattacker << endl;
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

		state.bug << OutcomeSms(th, myidx, enemyidx);

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
void TacticalSms::improve(uint theateridx, uint myidx, uint enemyidx)
{
	d.theaters[theateridx]->mymoves[myidx]->outcomes[enemyidx].improved++;

	Improve imp(*this, *d.theaters[theateridx], myidx, enemyidx);
	imp.do_improve();
}

static float hillvalue(const Submap & sm, const Location & pos, bool mine, float ValueHill)
{
	if (sm[pos] & Submap::Hill) {
		if (mine != bool(sm[pos] & Submap::Mine))
			return ValueHill * ValueHill;
		return 1.0;
	}

	for (int dir = 0; dir < TDIRECTIONS; ++dir) {
		Location n;
		if (!sm.getneighbour(pos, dir, n))
			continue;

		if (sm[n] & Submap::Hill) {
			if (mine != bool(sm[pos] & Submap::Mine))
				return ValueHill;
			return 1.0;
		}
	}

	return 1.0;
}

void TacticalSms::evaluate_moves(Theater & th, PlayerMove & mymove, PlayerMove & enemymove, float & myvalue, float & enemyvalue)
{
	uint myenemydist = 0;

	d.nrevals++;

	if (mymove.nrcollided)
		myvalue *= pow(ValueLoss, mymove.nrcollided);
	if (enemymove.nrcollided)
		enemyvalue *= pow(ValueLoss, enemymove.nrcollided);

	for (uint myantidx = 0; myantidx < mymove.antmoves.size(); ++myantidx) {
		PlayerMove::AntMove & am = mymove.antmoves[myantidx];
		if (am.collided || !th.basesm.inside(am.pos))
			continue;

		am.killed = 0;

		uint nrenemies = enemymove.map[am.pos] & PlayerMove::AttackMask;
		if (nrenemies != 0) {
			for (uint enemyantidx = 0; enemyantidx < enemymove.antmoves.size(); ++enemyantidx) {
				const PlayerMove::AntMove & em = enemymove.antmoves[enemyantidx];
				if (em.collided || !th.basesm.inside(em.pos))
					continue;

				if (th.basesm.eucliddist2(em.pos, am.pos) > AttackRadius2)
					continue;
				if ((mymove.map[em.pos] & PlayerMove::AttackMask) > nrenemies)
					continue;

				am.killed = 1;
				break;
			}
		}

		if (am.killed) {
			if (th.ismyperspective && th.aggressive)
				myvalue /= ValueKill;
			else
				myvalue *= ValueLoss;
			enemyvalue *= ValueKill;
		} else {
			float hv = hillvalue(th.basesm, am.pos, th.ismyperspective, ValueHill);
			myvalue *= hv;
			enemyvalue *= 1.0 / hv;

			if (th.ismyperspective)
				myenemydist += bot.m_zoc.m_enemy[state.addLocations(am.pos, th.offset)];
		}
	}

	if (th.ismyperspective)
		myvalue *= pow(ValueEnemyDist, myenemydist);

	for (uint enemyantidx = 0; enemyantidx < enemymove.antmoves.size(); ++enemyantidx) {
		PlayerMove::AntMove & am = enemymove.antmoves[enemyantidx];
		if (am.collided || !th.basesm.inside(am.pos))
			continue;

		am.killed = 0;

		uint nrenemies = mymove.map[am.pos] & PlayerMove::AttackMask;
		if (nrenemies != 0) {
			for (uint myantidx = 0; myantidx < mymove.antmoves.size(); ++myantidx) {
				const PlayerMove::AntMove & em = mymove.antmoves[myantidx];
				if (em.collided || !th.basesm.inside(em.pos))
					continue;

				if (th.basesm.eucliddist2(am.pos, em.pos) > AttackRadius2)
					continue;
				if ((enemymove.map[em.pos] & PlayerMove::AttackMask) > nrenemies)
					continue;

				am.killed = 1;
				break;
			}
		}

		if (am.killed) {
			myvalue *= ValueKill;
			if (!th.ismyperspective && th.aggressive)
				enemyvalue /= ValueKill;
			else
				enemyvalue *= ValueLoss;
		} else {
			float hv = hillvalue(th.basesm, am.pos, !th.ismyperspective, ValueHill);
			enemyvalue *= hv;
			myvalue *= 1.0 / hv;

			if (!th.ismyperspective)
				myenemydist += bot.m_zoc.m_enemy[state.addLocations(am.pos, th.offset)];
		}
	}

	if (!th.ismyperspective)
		enemyvalue *= pow(ValueEnemyDist, myenemydist);
}

void TacticalSms::evaluate_pair(uint theateridx, uint myidx, uint enemyidx)
{
	Theater & th = *d.theaters[theateridx];
	PlayerMove & mymove = *th.mymoves[myidx];
	PlayerMove & enemymove = *th.enemymoves[enemyidx];

	while (mymove.outcomes.size() <= enemyidx) {
		mymove.outcomes.push_back(PlayerMove::VsOutcome(numeric_limits<float>::max(), -1));
	}
	while (enemymove.outcomes.size() <= myidx) {
		enemymove.outcomes.push_back(PlayerMove::VsOutcome(numeric_limits<float>::max(), -1));
	}

	if (mymove.outcomes[enemyidx].improved >= 0) {
		assert(enemymove.outcomes[myidx].improved >= 0);
		return;
	}
	assert(enemymove.outcomes[myidx].improved < 0);

	float myvalue = 1.0;
	float enemyvalue = 1.0;

	evaluate_moves(th, mymove, enemymove, myvalue, enemyvalue);

	state.bug << " " << theateridx << ": " << myidx << " vs. " << enemyidx
		<< "  value " << myvalue << " vs " << enemyvalue << endl;

	mymove.outcomes[enemyidx].value = myvalue;
	mymove.outcomes[enemyidx].improved = 0;
	if (myvalue < mymove.worstvalue)
		mymove.worstvalue = myvalue;

	enemymove.outcomes[myidx].value = enemyvalue;
	enemymove.outcomes[myidx].improved = 0;
	if (enemyvalue < enemymove.worstvalue)
		enemymove.worstvalue = enemyvalue;
}

void TacticalSms::update_weights(uint theateridx)
{
	Theater & th = *d.theaters[theateridx];

	state.bug << "update weights in theater " << theateridx << ", perspective: " << (th.ismyperspective ? "mine" : "enemy") << endl;

	d.tmp_candidates.clear();
	d.tmp_mask.clear();

	float totalweight = 0.0;
	float totalcounterweight = 0.0;
	for (uint idx = 0; idx < th.enemymoves.size(); ++idx) {
		totalweight += th.enemymoves[idx]->weight;
		totalcounterweight += th.enemymoves[idx]->counterweight;
	}

	float u = fastrngd() * totalweight * (1.0 / 4.0);
	float v = fastrngd() * totalcounterweight * (1.0 / 4.0);

	for (uint enemyidx = 0; enemyidx < th.enemymoves.size(); ++enemyidx) {
		bool select = false;

		u -= th.enemymoves[enemyidx]->weight;
		if (u <= 0.0) {
			select = true;
			while (u <= 0.0)
				u += totalweight * (1.0 / 4.0);
		}

		v -= th.enemymoves[enemyidx]->counterweight;
		if (v <= 0.0) {
			select = true;
			while (v <= 0.0)
				v += totalweight * (1.0 / 4.0);
		}

		if (select) {
			state.bug << "    counter candidate: " << enemyidx << endl;
			d.tmp_candidates.push_back(enemyidx);
			d.tmp_mask.push_back(false);
		}
	}

	//
	float bestvalue = numeric_limits<float>::min();

	for (uint myidx = 0; myidx < th.mymoves.size(); ++myidx) {
		PlayerMove & mymove = *th.mymoves[myidx];

		for (uint idx = 0; idx < d.tmp_candidates.size(); ++idx)
			evaluate_pair(theateridx, myidx, d.tmp_candidates[idx]);

		for (uint idx = 0; idx < d.tmp_candidates.size(); ++idx) {
			if (mymove.outcomes[d.tmp_candidates[idx]].value <= mymove.worstvalue + EpsilonValue)
				d.tmp_mask[idx] = true;
		}

		bestvalue = max(bestvalue, mymove.worstvalue);
	}

	//
	state.bug << "  best value: " << bestvalue << endl;

	for (uint myidx = 0; myidx < th.mymoves.size(); ++myidx) {
		PlayerMove & mymove = *th.mymoves[myidx];

		if (mymove.worstvalue >= bestvalue - EpsilonValue) {
			mymove.weight *= LearnGamma;
			state.bug << "    move " << myidx << ", new weight: " << mymove.weight << endl;
		}
	}

	for (uint idx = 0; idx < d.tmp_candidates.size(); ++idx) {
		if (!d.tmp_mask[idx]) {
			PlayerMove & enemymove = *th.enemymoves[d.tmp_candidates[idx]];
			enemymove.counterweight /= LearnGamma;
			state.bug << "    enemy move " << d.tmp_candidates[idx]
				<< ", new counter-weight: " << enemymove.counterweight << endl;
		}
	}
}

void TacticalSms::generate_theater(const Location & center)
{
	uint theateridx = d.theaters.size();
	Theater & th = *d.alloctheater();
	d.theaters.push_back(&th);

	state.bug << "Tactical theater " << theateridx << " around " << center << endl;

	gensubmap(th.basesm, center);
	th.offset = Location
		((center.row - Submap::Radius + state.rows) % state.rows,
		 (center.col - Submap::Radius + state.cols) % state.cols);

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

	// determine aggressiveness based on ratio of our vs. enemy ants
	if (th.myants.size() >= MinAggressionAnts) {
		double antratio = (double)th.myants.size() / (double)th.enemyants.size();
		double p = 1.0 / (1.0 + exp(-AggressionScale * (log(antratio) + AggressionThresholdShift)));

		if (fastrngd() <= p)
			th.aggressive = true;
	}

	// generate enemy init move
	th.enemymoves.push_back(d.allocmove());
	d.nrmoves++;
	th.enemymoves[0]->map.fill(0);

	th.enemymoves[0]->antmoves.reserve(th.enemyants.size());
	for (uint idx = 0; idx < th.enemyants.size(); ++idx) {
		const Location & local = th.enemyants[idx];
		th.enemymoves[0]->antmoves.push_back(PlayerMove::AntMove(local, -1));
		th.enemymoves[0]->ant_mark(idx);
	}

	th.enemymoves[0]->computehash();
}

void TacticalSms::pull_moves(uint theateridx)
{
	Theater & th = *d.theaters[theateridx];

	state.bug << "Theater " << theateridx << ": pulling  new moves " << th.mymoves.size() << endl;

	th.mymoves.push_back(d.allocmove());

	PlayerMove & pm = *th.mymoves.back();

	pm.map.fill(0);
	pm.antmoves.reserve(th.myants.size());
	for (uint idx = 0; idx < th.myants.size(); ++idx) {
		const Location & local = th.myants[idx];
		Location global
			((local.row + th.offset.row) % state.rows,
			 (local.col + th.offset.col) % state.cols);
		Ant & ant = bot.m_ants[bot.myantidx_at(global)];
		Location nextlocal;
		th.basesm.getneighbouropt(local, ant.direction, nextlocal);
		pm.antmoves.push_back(PlayerMove::AntMove(nextlocal, ant.direction));
		pm.ant_mark(idx);
	}

	pm.computehash();

	if (th.is_duplicate_mymove(th.mymoves.size() - 1)) {
		d.unusedmoves.push_back(&pm);
		th.mymoves.pop_back();
		state.bug << "  no new moves" << endl;
	} else {
		d.nrmoves++;
		th.needupdate = true;
		th.mymove = th.mymoves.size() - 1;
		th.mymove_unchangedrounds = 0;
	}
}

bool TacticalSms::get_improve_partner(const vector<PlayerMove *> & moves, uint myidx, uint & enemyidx)
{
	PlayerMove & mymove = *moves[myidx];

	d.tmp_candidates.clear();
	for (uint idx = 0; idx < mymove.outcomes.size(); ++idx) {
		if (!mymove.outcomes[idx].improved && mymove.outcomes[idx].value <= mymove.worstvalue + EpsilonValue)
			d.tmp_candidates.push_back(idx);
	}

	if (!d.tmp_candidates.empty()) {
		enemyidx = d.tmp_candidates[fastrng() % d.tmp_candidates.size()];
		return true;
	}

	return false;
}

void TacticalSms::run_theater(uint theateridx)
{
	state.bug << "Run tactical theater " << theateridx << " (turn " << state.turn << ")" << endl;

	Theater & th = *d.theaters[theateridx];

	uint init_mymoves = th.mymoves.size();
	uint init_enemymoves = th.enemymoves.size();

	// Phase 1: find improved moves and update weights
	uint mine_myidx = th.mymove;
	uint mine_enemyidx = 0;
	bool mine;

	for (uint idx = 0; idx < th.enemymoves.size(); ++idx)
		evaluate_pair(theateridx, mine_myidx, idx);
	mine = get_improve_partner(th.mymoves, mine_myidx, mine_enemyidx);

	//
	uint enemy_myidx = 0;
	uint enemy_enemyidx;
	bool enemy = false;

	float totalweight = 0.0;

	for (uint idx = 0; idx < th.enemymoves.size(); ++idx)
		totalweight += th.enemymoves[idx]->weight;

	for (uint tries = 0; tries < 3 && !enemy; ++tries) {
		float v = fastrngd() * totalweight;

		enemy_enemyidx = 0;
		do {
			v -= th.enemymoves[enemy_enemyidx++]->weight;
		} while (v >= 0.0 && enemy_enemyidx < th.enemymoves.size());
		enemy_enemyidx--;

		for (uint idx = 0; idx < th.mymoves.size(); ++idx)
			evaluate_pair(theateridx, idx, enemy_enemyidx);
		enemy = get_improve_partner(th.enemymoves, enemy_enemyidx, enemy_myidx);
	}

	//
	if (mine || enemy) {
		th.flipsides();
		if (mine)
			improve(theateridx, mine_enemyidx, mine_myidx);
		if (enemy && (!mine || enemy_enemyidx != mine_enemyidx || enemy_myidx != mine_myidx))
			improve(theateridx, enemy_enemyidx, enemy_myidx);
		update_weights(theateridx);
		th.flipsides();
		if (mine)
			improve(theateridx, mine_myidx, mine_enemyidx);
		if (enemy && (!mine || enemy_enemyidx != mine_enemyidx || enemy_myidx != mine_myidx))
			improve(theateridx, enemy_myidx, enemy_enemyidx);
		update_weights(theateridx);
	} else {
		th.flipsides();
		update_weights(theateridx);
		th.flipsides();
		update_weights(theateridx);
	}

	// Phase 2: push new sampled move
	if (init_mymoves != th.mymoves.size() || init_enemymoves != th.enemymoves.size())
		th.mymove_unchangedrounds = 0;

	totalweight = 0.0;
	for (uint idx = 0; idx < th.mymoves.size(); ++idx)
		totalweight += th.mymoves[idx]->weight;

	float v = fastrngd() * totalweight;
	uint myidx = 0;
	do {
		v -= th.mymoves[myidx++]->weight;
	} while (v >= 0.0 && myidx < th.mymoves.size());
	myidx--;

	state.bug << "Choosing " << myidx << endl;
	push_moves(theateridx, myidx);

	if (myidx == th.mymove) {
		th.mymove_unchangedrounds++;
		if (th.mymove_unchangedrounds >= 4)
			th.needupdate = false;
	} else {
		th.mymove = myidx;
		th.mymove_unchangedrounds = 0;
	}

	state.bug << "-----------------------" << endl;
}

void TacticalSms::push_moves(uint theateridx, uint myidx)
{
	Theater & th = *d.theaters[theateridx];
	const PlayerMove & move = *th.mymoves[myidx];

	for (uint idx = 0; idx < move.antmoves.size(); ++idx) {
		int dir = move.antmoves[idx].direction;
		Location antpos = th.myants[idx];
		Location global
			((antpos.row + th.offset.row) % state.rows,
			(antpos.col + th.offset.col) % state.cols);
		uint antidx = bot.myantidx_at(global);
		Ant & ant = bot.m_ants[antidx];

		ant.direction = dir;
	}
}

bool TacticalSms::timeover()
{
	return bot.timeover();
#if 0
	if (d.nrevals > MaxEval) {
		state.bug << "Max number of evals exceeded" << endl;
		return true;
	}

	return false;
#endif
}

void TacticalSms::run()
{
	uint nrtheaters = 0;

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
			pull_moves(theateridx);

			if (d.theaters[theateridx]->needupdate)
				break;
			theateridx = (theateridx + 1) % d.theaters.size();
		} while (theateridx != starttheater);

		if (!d.theaters[theateridx]->needupdate)
			break;

		run_theater(theateridx);
		nrtheaters++;
	}

	state.bug.time << "Total number of generated moves: " << d.nrmoves
		<< ", evaluated pairs: " << d.nrevals
		<< ", theater rounds: " << nrtheaters
		<< " in " << d.theaters.size() << " theaters" << endl;
}
