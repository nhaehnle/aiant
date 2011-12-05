#include "symmetry.h"

#include <map>
#include <set>

#include "Bot.h"

using namespace std;

static const uint ORIENTATIONS = 8;

Location orient_offset(const Location & ofs, uint orientation) {
	switch (orientation) {
	case 0: return ofs;
	case 1: return Location(-ofs.row, ofs.col);
	case 2: return Location(ofs.row, -ofs.col);
	case 3: return Location(-ofs.row, -ofs.col);
	case 4: return Location(ofs.col, ofs.row);
	case 5: return Location(-ofs.col, ofs.row);
	case 6: return Location(ofs.col, -ofs.row);
	case 7: return Location(-ofs.col, -ofs.row);
	}
	abort();
}

/**
 * Return the orientation f . g, where . is function composition.
 */
uint compose_orientations(uint f, uint g)
{
	uint out = g;
	if (f & 4) {
		out = ((out & 1) << 1) | ((out & 2) >> 1) | ((out & 4) ^ 4);
	}
	out ^= f & 3;
	return out;
}

/**
 * Represent the function (aka symmetry) that maps (0,0) at orientation 0
 * onto @ref origin with given @ref orientation.
 */
struct Symmetry {
	Location origin;
	uint orientation;

	Symmetry(const Location & origin_, uint orientation_) :
		origin(origin_),
		orientation(orientation_)
	{
	}

	bool isidentity() const {return origin.row == 0 && origin.col == 0 && orientation == 0;}

	bool operator==(const Symmetry & o) const {
		return origin == o.origin && orientation == o.orientation;
	}
	bool operator!=(const Symmetry & o) const {return !(*this == o);}

	Location apply(const State & state, const Location & loc) const {
		Location oriented = orient_offset(loc, orientation);
		uint shift = state.rows * state.cols;
		return Location
			((origin.row + oriented.row + shift) % state.rows,
			 (origin.col + oriented.col + shift) % state.cols);
	}
};

/**
 * return (f . g), where . is composition of functions.
 */
Symmetry compose_symmetries(const State & state, const Symmetry & f, const Symmetry & g) const
{
	Symmetry out(f.apply(state, g.origin), compose_orientations(f.orientation, g.orientation));
}

struct SymmetryFinder::Data {
	vector<Location> all_my_hills;
	vector<Location> all_enemy_hills;

	multimap<uint32_t, Location> fingerprint_multimap;
	Map<uint32_t> fingerprints;
	vector<Location> fingerprintstack;
};

SymmetryFinder::SymmetryFinder(Bot & b) :
	bot(b),
	state(b.state),
	d(*new Data),
	newwater(false)
{
}

SymmetryFinder::~SymmetryFinder()
{
	delete &d;
}

void SymmetryFinder::init()
{
	map.resize(state.rows, state.cols);
	d.fingerprints.resize(state.rows, state.cols);
}

void SymmetryFinder::update_hills()
{
	for (uint idx = 0; idx < state.myHills.size(); ++idx) {
		const Location & pos = state.myHills[idx];

		for (uint j = 0; j < d.all_my_hills.size(); ++j) {
			if (d.all_my_hills[j] == pos)
				goto found;
		}

		d.all_my_hills.push_back(pos);
	found: ;
	}

	for (uint idx = 0; idx < state.enemyHills.size(); ++idx) {
		const Location & pos = state.enemyHills[idx];

		for (uint j = 0; j < d.all_enemy_hills.size(); ++j) {
			if (d.all_enemy_hills[j] == pos)
				goto found;
		}

		d.all_enemy_hills.push_back(pos);
	found: ;
	}
}

void SymmetryFinder::compute_fingerprints()
{
	state.bug << "symmetry: compute new fingerprints" << endl;

	Location center;
	for (center.row = 0; center.row < state.rows; ++center.row) {
		for (center.col = 0; center.col < state.cols; ++center.col) {
			if (map[center] & MapFingerprinted)
				continue;
			if (!have_seen(center, 3))
				continue;

			uint32_t fingerprint = 0;
			const uint32_t prime = 312101;
			uint nrwater = 0;
			bool centerwater = state.grid[center.row][center.col].isWater;

			for (uint d = 1; d <= 3; ++d) {
				uint count1 = 0;
				uint count2 = 0;
				static const uint orientmap[4] = { 0, 2, 4, 5 };
				for (uint orient = 0; orient < 4; ++orient) {
					Location n(state.addLocations(center, orient_offset(Location(d, d), orient)));
					count1 += uint(state.grid[n.row][n.col].isWater != centerwater);

					n = state.addLocations(center, orient_offset(Location(0, d), orientmap[orient]));
					count2 += uint(state.grid[n.row][n.col].isWater != centerwater);
				}

				fingerprint = (fingerprint * prime) + count1;
				fingerprint = (fingerprint * prime) + count2;

				for (uint k = 1; k < d; ++k) {
					Location n[8];
					uint counts[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
					for (uint orient = 0; orient < ORIENTATIONS; ++orient)
						n[orient] = state.addLocations(center, orient_offset(Location(d, k), orient));

					for (uint orient = 0; orient < ORIENTATIONS; ++orient) {
						bool nwater = state.grid[n[orient].row][n[orient].col].isWater;
						counts[0] += uint(nwater != centerwater);

						for (uint rel = 1; rel < ORIENTATIONS; ++rel) {
							uint comp = compose_orientations(orient, rel);
							counts[rel] += uint(nwater != state.grid[n[comp].row][n[comp].col].isWater);
						}
					}

					for (uint orient = 0; orient < ORIENTATIONS; ++orient)
						fingerprint = (fingerprint * prime) + counts[orient];
				}
			}

			state.bug << "  fingerprint: " << center << " : " << fingerprint << endl;

			map[center] |= MapFingerprinted;
			d.fingerprints[center] = fingerprint;
			if (fingerprint)
				d.fingerprint_multimap.insert(std::make_pair(fingerprint, center));
			d.fingerprintstack.push_back(center);
		}
	}
}

void SymmetryFinder::find_fingerprint_symmetries(const Location & center)
{
	if (do_find_fingerprint_symmetries(center))
		return;

	const int * dirperm = getdirperm();
	for (int predir = 0; predir < TDIRECTIONS; ++predir) {
		int dir = dirperm[predir];
		if
			(do_find_fingerprint_symmetries
			 (state.addLocations(center, Location(3 * DIRECTIONS[dir][0], 3 * DIRECTIONS[dir][1]))))
			return;
	}
}

void SymmetryFinder::run()
{
	check_destroyed_hills();

	if (state.newsquare) {
		update_hills();

		compute_fingerprints();

		while (!d.fingerprintstack.empty()) {
			find_fingerprint_symmetries(d.fingerprintstack.back());
			d.fingerprintstack.pop_back();
		}
	}
}
