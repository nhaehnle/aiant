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

#include "symmetry.h"

#include <map>
#include <set>

#include "Bot.h"

using namespace std;

static const uint ORIENTATIONS = 8;
static const int ORIENTATIONDIRPERMS[8][4] = {
	{ 0, 1, 2, 3 }, // ORIENTATIONDIRPERMS[orientation][dir] = orientation . dir
	{ 2, 1, 0, 3 },
	{ 0, 3, 2, 1 },
	{ 2, 3, 0, 1 },
	{ 3, 2, 1, 0 },
	{ 1, 2, 3, 0 },
	{ 3, 0, 1, 2 },
	{ 1, 0, 3, 2 }
};

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

	/**
	 * Symmetry that maps @p from to @p to with given @p orientation_.
	 *
	 * Must satisfy @code to = origin + orient_offset(from, orientation) @endcode
	 */
	Symmetry(const State & state, const Location & from, const Location & to, uint orientation_) :
		orientation(orientation_)
	{
		Location oriented = orient_offset(from, orientation);
		uint shift = state.rows * state.cols;
		origin.row = (to.row - oriented.row + shift) % state.rows;
		origin.col = (to.col - oriented.col + shift) % state.cols;
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

ostream & operator<<(ostream & out, const Symmetry & s)
{
	return out << "[" << s.origin << ", " << s.orientation << "]";
}

/**
 * return (f . g), where . is composition of functions.
 */
Symmetry compose_symmetries(const State & state, const Symmetry & f, const Symmetry & g)
{
	return Symmetry(f.apply(state, g.origin), compose_orientations(f.orientation, g.orientation));
}

struct SymmetryCompare {
	bool operator()(const Symmetry & a, const Symmetry & b) const {
		if (a.origin.row < b.origin.row)
			return true;
		if (a.origin.row == b.origin.row) {
			if (a.origin.col < b.origin.col)
				return true;
			if (a.origin.col == b.origin.col)
				return a.orientation < b.orientation;
		}
		return false;
	}
};

typedef set<Symmetry, SymmetryCompare> SymmetrySet;

struct SymmetryFinder::Data {
	typedef multimap<uint32_t, Location> FingerprintMap;

	vector<Location> all_my_hills; // all hills that have ever actually been seen
	vector<Location> all_enemy_hills; // all hills that have ever actually been seen
	SymmetrySet symmetries;
	SymmetrySet rejected_symmetries;

	FingerprintMap fingerprint_multimap;
	Map<uint32_t> fingerprints;
	vector<Location> fingerprintstack;
};

SymmetryFinder::SymmetryFinder(Bot & b) :
	bot(b),
	state(b.state),
	d(*new Data),
	newwater(false),
	enemyhillschanged(false)
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

ostream & operator<<(ostream & out, const SymmetryFinder & s)
{
	Location cur;
	for (cur.row = 0; cur.row < s.state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < s.state.cols; ++cur.col) {
			Square & sq = s.state.grid[cur.row][cur.col];

			if (s.map[cur] & SymmetryFinder::MapWater) {
				if (sq.lastseen > 0)
					out << '%';
				else
					out << '@';
			} else {
				if (sq.hillPlayer >= 0) {
					out << (char)('A' + sq.hillPlayer);
				} else if (sq.ant >= 0) {
					out << (char)('a' + sq.ant);
				} else if (s.map[cur] & SymmetryFinder::MapEnemyHill) {
					out << 'H';
				} else {
					if (sq.isVisible)
						out << ' ';
					else if (sq.lastseen > 0)
						out << '.';
					else if (s.map[cur] & SymmetryFinder::MapKnown)
						out << ',';
					else
						out << '?';
				}
			}
		}
		out << endl;
	}
	return out;
}

void SymmetryFinder::add_possible_enemy_hill(const Location & pos)
{
	if (map[pos] & MapKnownNoHill)
		return;

	for (uint idx = 0; idx < d.all_my_hills.size(); ++idx) {
		if (d.all_my_hills[idx] == pos)
			return;
	}

	for (uint idx = 0; idx < d.all_enemy_hills.size(); ++idx) {
		if (d.all_enemy_hills[idx] == pos)
			return;
	}

	enemy_hills.push_back(pos);
	enemyhillschanged = true;
	map[pos] |= MapEnemyHill;
}

void SymmetryFinder::broadcast_hill(const Location & pos)
{
	for (set<Symmetry>::const_iterator it = d.symmetries.begin(); it != d.symmetries.end(); ++it) {
		Location img = it->apply(state, pos);

		add_possible_enemy_hill(img);
	}
}

void SymmetryFinder::recheck_symmetries()
{
	Location cur;
	for (cur.row = 0; cur.row < state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < state.cols; ++cur.col) {
			map[cur] &= ~MapKnown;
		}
	}

	SymmetrySet oldsymmetries;
	swap(d.symmetries, oldsymmetries);

	for (SymmetrySet::const_iterator it = oldsymmetries.begin(); it != oldsymmetries.end(); ++it) {
		state.bug << "recheck symmetry " << *it << endl;

		if (check_symmetry(*it)) {
			state.bug << "  good" << endl;
			d.symmetries.insert(*it);
		} else {
			state.bug << "  bad" << endl;
			d.rejected_symmetries.insert(*it);
		}
	}
}

void SymmetryFinder::update_map()
{
	// broadcast newly known fields
	Location cur;
	for (cur.row = 0; cur.row < state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < state.cols; ++cur.col) {
			Square & sq = state.grid[cur.row][cur.col];
			if (sq.lastseen == 0)
				continue;

			if (map[cur] & MapKnown) {
				if (sq.isWater != bool(map[cur] & MapWater)) {
					state.bug << "  OUCH! incorrect symmetry detected at " << cur << endl;
					recheck_symmetries();
				}
			} else {
				map[cur] |= MapKnown;
				if (sq.isWater) {
					map[cur] |= MapWater;
					newwater = true;
				} else {
					map[cur] &= ~MapWater;
				}

				for (set<Symmetry>::const_iterator it = d.symmetries.begin(); it != d.symmetries.end(); ++it) {
					Location img = it->apply(state, cur);
					if (map[img] & MapKnown) {
						if ((map[img] & MapWater) != (map[cur] & MapWater)) {
							state.bug << "  OUCH! incorrect symmetry detected at " << img << " from " << cur << endl;
							recheck_symmetries();
							break;
						}
					} else {
						map[img] = (map[img] & ~MapWater) | MapKnown | (map[cur] & MapWater);
					}
				}
			}
		}
	}

	// broadcast newly found hills
	bool newhills = false;
	for (uint idx = 0; idx < state.myHills.size(); ++idx) {
		const Location & pos = state.myHills[idx];

		for (uint j = 0; j < d.all_my_hills.size(); ++j) {
			if (d.all_my_hills[j] == pos)
				goto foundmine;
		}

		d.all_my_hills.push_back(pos);
		newhills = true;
	foundmine: ;
	}

	// broadcast enemy hills
	for (uint idx = 0; idx < state.enemyHills.size(); ++idx) {
		const Location & pos = state.enemyHills[idx];

		for (uint j = 0; j < d.all_enemy_hills.size(); ++j) {
			if (d.all_enemy_hills[j] == pos)
				goto foundenemy;
		}

		d.all_enemy_hills.push_back(pos);
		enemy_hills.push_back(pos);
		enemyhillschanged = true;
		map[pos] |= MapEnemyHill;
		newhills = true;
	foundenemy: ;
	}

	// do this afterwards, in case a symmetry maps own hills onto each other
	if (newhills) {
		for (uint idx = 0; idx < state.myHills.size(); ++idx)
			broadcast_hill(state.myHills[idx]);
		for (uint idx = 0; idx < state.enemyHills.size(); ++idx)
			broadcast_hill(state.enemyHills[idx]);
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
			bool centerwater = state.grid[center.row][center.col].isWater;

			for (uint dist = 1; dist <= 3; ++dist) {
				uint count1 = 0;
				uint count2 = 0;
				static const uint orientmap[4] = { 0, 2, 4, 5 };
				for (uint orient = 0; orient < 4; ++orient) {
					Location n(state.addLocations(center, orient_offset(Location(dist, dist), orient)));
					count1 += uint(state.grid[n.row][n.col].isWater != centerwater);

					n = state.addLocations(center, orient_offset(Location(0, dist), orientmap[orient]));
					count2 += uint(state.grid[n.row][n.col].isWater != centerwater);
				}

				fingerprint = (fingerprint * prime) + count1;
				fingerprint = (fingerprint * prime) + count2;

				for (uint k = 1; k < dist; ++k) {
					Location n[8];
					uint counts[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
					for (uint orient = 0; orient < ORIENTATIONS; ++orient)
						n[orient] = state.addLocations(center, orient_offset(Location(dist, k), orient));

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

			map[center] |= MapFingerprinted;
			d.fingerprints[center] = fingerprint;
			if (fingerprint)
				d.fingerprint_multimap.insert(std::make_pair(fingerprint, center));
			d.fingerprintstack.push_back(center);
		}
	}
}

void SymmetryFinder::add_symmetry(const Symmetry & s)
{
	if (d.symmetries.count(s))
		return;

	state.bug << "Adding symmetry " << s << endl;
	d.symmetries.insert(s);

	Location cur;
	for (cur.row = 0; cur.row < state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < state.cols; ++cur.col) {
			if (!(map[cur] & MapKnown))
				continue;

			Location img = s.apply(state, cur);
			if (!(map[img] & MapKnown)) {
				map[img] |= MapKnown | (map[cur] & MapWater);
			}
		}
	}

	for (uint idx = 0; idx < d.all_my_hills.size(); ++idx) {
		Location img = s.apply(state, d.all_my_hills[idx]);
		add_possible_enemy_hill(img);
	}
	for (uint idx = 0; idx < d.all_enemy_hills.size(); ++idx) {
		Location img = s.apply(state, d.all_enemy_hills[idx]);
		add_possible_enemy_hill(img);
	}
}

bool SymmetryFinder::check_symmetry(const Symmetry & s)
{
	Location cur;
	for (cur.row = 0; cur.row < state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < state.cols; ++cur.col) {
			Square & sq = state.grid[cur.row][cur.col];
			if (sq.lastseen == 0)
				continue;

			Location imgpos = s.apply(state, cur);
			Square & img = state.grid[imgpos.row][imgpos.col];

			if (img.lastseen == 0)
				continue;

			if (sq.isWater != img.isWater)
				return false;
		}
	}
	return true;
}

void SymmetryFinder::add_candidate_symmetry(const Location & from, const Location & to, uint orientation)
{
	Symmetry s(state, from, to, orientation);

	if (d.symmetries.count(s) || d.rejected_symmetries.count(s))
		return;

	state.bug << "Check candidate symmetry from " << from << " to " << to << " orientation " << orientation << ": " << s << endl;

	Symmetry f(s);
	uint power = 1;

	while (!f.isidentity()) {
		state.bug << "  check " << f << endl;

		if (!check_symmetry(f)) {
			state.bug << "  not valid" << endl;
			d.rejected_symmetries.insert(s);
			if (power > 1)
				d.rejected_symmetries.insert(f);
			return;
		}

		power++;
		if (power > 20) {
			state.bug << "  suspicious: symmetry has order > 20, skipping" << endl;
			d.rejected_symmetries.insert(s);
			return;
		}
		f = compose_symmetries(state, f, s);
	}

	state.bug << "  checks successful, adding symmetries based on " << s << endl;

	f = s;
	while (!f.isidentity()) {
		add_symmetry(f);
		f = compose_symmetries(state, f, s);
	}
}

bool SymmetryFinder::do_find_fingerprint_symmetries(const Location & center)
{
	uint32_t fingerprint = d.fingerprints[center];
	uint32_t fingerprints[4];

	if (!map[center] & MapFingerprinted)
		return false;
	if (fingerprint == 0)
		return false;

	for (int dir = 0; dir < TDIRECTIONS; ++dir) {
		Location n(state.addLocations(center, Location(3 * DIRECTIONS[dir][0], 3 * DIRECTIONS[dir][1])));
		if (!(map[n] & MapFingerprinted))
			return false;

		fingerprints[dir] = d.fingerprints[n];
	}

	for
		(Data::FingerprintMap::const_iterator it = d.fingerprint_multimap.lower_bound(fingerprint);
		 it != d.fingerprint_multimap.end() && it->first == fingerprint;
		 ++it)
	{
		if (it->second == center)
			continue;

		uint otherprints[4];

		for (int dir = 0; dir < TDIRECTIONS; ++dir) {
			Location n(state.addLocations(it->second, Location(3 * DIRECTIONS[dir][0], 3 * DIRECTIONS[dir][1])));
			if (!(map[n] & MapFingerprinted))
				goto notfingerprinted;
			otherprints[dir] = d.fingerprints[n];
		}

		for (uint orientation = 0; orientation < ORIENTATIONS; ++orientation) {
			for (int dir = 0; dir < TDIRECTIONS; ++dir) {
				if (otherprints[ORIENTATIONDIRPERMS[orientation][dir]] != fingerprints[dir])
					goto skip;
			}

			add_candidate_symmetry(center, it->second, orientation);
		skip: ;
		}

	notfingerprinted: ;
	}

	return true;
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

void SymmetryFinder::check_destroyed_hills()
{
	Location cur;
	for (cur.row = 0; cur.row < state.rows; ++cur.row) {
		for (cur.col = 0; cur.col < state.cols; ++cur.col) {
			Square & sq = state.grid[cur.row][cur.col];
			if (sq.isVisible && sq.hillPlayer < 0)
				map[cur] = (map[cur] & ~MapEnemyHill) | MapKnownNoHill;
		}
	}

	for (uint idx = 0; idx < enemy_hills.size(); ++idx) {
		if (map[enemy_hills[idx]] & MapKnownNoHill) {
			enemy_hills[idx] = enemy_hills.back();
			enemy_hills.pop_back();
			enemyhillschanged = true;
			idx--;
		}
	}
}

bool SymmetryFinder::have_seen(const Location & center, uint rel) const
{
	Location ofs;
	for (ofs.row = -(int)rel; ofs.row <= (int)rel; ++ofs.row) {
		for (ofs.col = -(int)rel; ofs.col <= (int)rel; ++ofs.col) {
			Location where = state.addLocations(center, ofs);
			if (state.grid[where.row][where.col].lastseen == 0)
				return false;
		}
	}
	return true;
}

void SymmetryFinder::run()
{
	newwater = false;
	enemyhillschanged = false;

	state.bug << "Symmetry turn " << state.turn << endl;

	check_destroyed_hills();

	if (state.newsquare) {
		update_map();

		compute_fingerprints();

		while (!d.fingerprintstack.empty()) {
			find_fingerprint_symmetries(d.fingerprintstack.back());
			d.fingerprintstack.pop_back();
		}
	}

	//state.bug << *this;
}
