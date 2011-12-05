#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <cstdlib>

#include "map.h"

struct Bot;
struct State;

struct SymmetryFinder {
	static const uint8_t MapWater = 0x01;
	static const uint8_t MapKnownNoHill = 0x40;
	static const uint8_t MapFingerprinted = 0x80;

	SymmetryFinder(Bot & b);
	~SymmetryFinder();

	void init();
	void run();

	void update_hills();
	void compute_fingerprints();
	void find_fingerprint_symmetries(const Location & center);

	bool have_seen(const Location & center, uint rel) const;

	struct Data;

	Bot & bot;
	State & state;
	Data & d;

	Map<uint8_t> map;
	bool newwater;
	std::vector<Location> enemy_hills;
};

#endif // SYMMETRY_H
