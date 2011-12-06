#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <cstdlib>

#include "map.h"

struct Bot;
struct State;
struct Symmetry;

struct SymmetryFinder {
	static const uint8_t MapWater = 0x01;
	static const uint8_t MapEnemyHill = 0x02;
	static const uint8_t MapKnown = 0x20;
	static const uint8_t MapKnownNoHill = 0x40;
	static const uint8_t MapFingerprinted = 0x80;

	SymmetryFinder(Bot & b);
	~SymmetryFinder();

	void init();
	void run();

	void update_map();
	void compute_fingerprints();
	bool do_find_fingerprint_symmetries(const Location & center);
	void find_fingerprint_symmetries(const Location & center);
	void add_candidate_symmetry(const Location & from, const Location & to, uint orientation);
	bool check_symmetry(const Symmetry & s);
	void recheck_symmetries();
	void add_symmetry(const Symmetry & s);
	void check_destroyed_hills();
	void broadcast_hill(const Location & pos);
	void add_possible_enemy_hill(const Location & pos);

	bool have_seen(const Location & center, uint rel) const;

	struct Data;

	Bot & bot;
	State & state;
	Data & d;

	Map<uint8_t> map;
	bool newwater;
	bool enemyhillschanged;
	std::vector<Location> enemy_hills;
};

#endif // SYMMETRY_H
