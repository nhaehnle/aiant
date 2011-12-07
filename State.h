#ifndef STATE_H_
#define STATE_H_

#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <stdint.h>

#include "Timer.h"
#include "Bug.h"
#include "Square.h"
#include "Location.h"

/*
    constants
*/
const int TDIRECTIONS = 4;
const char CDIRECTIONS[4] = {'N', 'E', 'S', 'W'};
const int DIRECTIONS[4][2] = { {-1, 0}, {0, 1}, {1, 0}, {0, -1} };      //{N, E, S, W}

inline int reversedir(int dir) {return (dir + 2) % 4;}
inline char cdir(int dir) {return (dir >= 0) ? CDIRECTIONS[dir] : '-';}

const int DIRPERMUTATIONS[24][4] = {
	{ 0, 1, 2, 3 },
	{ 0, 1, 3, 2 },
	{ 0, 2, 1, 3 },
	{ 0, 2, 3, 1 },
	{ 0, 3, 1, 2 },
	{ 0, 3, 2, 1 },
	{ 1, 0, 2, 3 },
	{ 1, 0, 3, 2 },
	{ 1, 2, 0, 3 },
	{ 1, 2, 3, 0 },
	{ 1, 3, 0, 2 },
	{ 1, 3, 2, 0 },
	{ 2, 0, 1, 3 },
	{ 2, 0, 3, 1 },
	{ 2, 1, 0, 3 },
	{ 2, 1, 3, 0 },
	{ 2, 3, 0, 1 },
	{ 2, 3, 1, 0 },
	{ 3, 0, 1, 2 },
	{ 3, 0, 2, 1 },
	{ 3, 1, 0, 2 },
	{ 3, 1, 2, 0 },
	{ 3, 2, 0, 1 },
	{ 3, 2, 1, 0 },
};

const int OPTDIRPERMUTATIONS[120][5] = {
	{ -1, 0, 1, 2, 3 },
	{ -1, 0, 1, 3, 2 },
	{ -1, 0, 2, 1, 3 },
	{ -1, 0, 2, 3, 1 },
	{ -1, 0, 3, 1, 2 },
	{ -1, 0, 3, 2, 1 },
	{ -1, 1, 0, 2, 3 },
	{ -1, 1, 0, 3, 2 },
	{ -1, 1, 2, 0, 3 },
	{ -1, 1, 2, 3, 0 },
	{ -1, 1, 3, 0, 2 },
	{ -1, 1, 3, 2, 0 },
	{ -1, 2, 0, 1, 3 },
	{ -1, 2, 0, 3, 1 },
	{ -1, 2, 1, 0, 3 },
	{ -1, 2, 1, 3, 0 },
	{ -1, 2, 3, 0, 1 },
	{ -1, 2, 3, 1, 0 },
	{ -1, 3, 0, 1, 2 },
	{ -1, 3, 0, 2, 1 },
	{ -1, 3, 1, 0, 2 },
	{ -1, 3, 1, 2, 0 },
	{ -1, 3, 2, 0, 1 },
	{ -1, 3, 2, 1, 0 },
	{ 0, -1, 1, 2, 3 },
	{ 0, -1, 1, 3, 2 },
	{ 0, -1, 2, 1, 3 },
	{ 0, -1, 2, 3, 1 },
	{ 0, -1, 3, 1, 2 },
	{ 0, -1, 3, 2, 1 },
	{ 0, 1, -1, 2, 3 },
	{ 0, 1, -1, 3, 2 },
	{ 0, 1, 2, -1, 3 },
	{ 0, 1, 2, 3, -1 },
	{ 0, 1, 3, -1, 2 },
	{ 0, 1, 3, 2, -1 },
	{ 0, 2, -1, 1, 3 },
	{ 0, 2, -1, 3, 1 },
	{ 0, 2, 1, -1, 3 },
	{ 0, 2, 1, 3, -1 },
	{ 0, 2, 3, -1, 1 },
	{ 0, 2, 3, 1, -1 },
	{ 0, 3, -1, 1, 2 },
	{ 0, 3, -1, 2, 1 },
	{ 0, 3, 1, -1, 2 },
	{ 0, 3, 1, 2, -1 },
	{ 0, 3, 2, -1, 1 },
	{ 0, 3, 2, 1, -1 },
	{ 1, -1, 0, 2, 3 },
	{ 1, -1, 0, 3, 2 },
	{ 1, -1, 2, 0, 3 },
	{ 1, -1, 2, 3, 0 },
	{ 1, -1, 3, 0, 2 },
	{ 1, -1, 3, 2, 0 },
	{ 1, 0, -1, 2, 3 },
	{ 1, 0, -1, 3, 2 },
	{ 1, 0, 2, -1, 3 },
	{ 1, 0, 2, 3, -1 },
	{ 1, 0, 3, -1, 2 },
	{ 1, 0, 3, 2, -1 },
	{ 1, 2, -1, 0, 3 },
	{ 1, 2, -1, 3, 0 },
	{ 1, 2, 0, -1, 3 },
	{ 1, 2, 0, 3, -1 },
	{ 1, 2, 3, -1, 0 },
	{ 1, 2, 3, 0, -1 },
	{ 1, 3, -1, 0, 2 },
	{ 1, 3, -1, 2, 0 },
	{ 1, 3, 0, -1, 2 },
	{ 1, 3, 0, 2, -1 },
	{ 1, 3, 2, -1, 0 },
	{ 1, 3, 2, 0, -1 },
	{ 2, -1, 0, 1, 3 },
	{ 2, -1, 0, 3, 1 },
	{ 2, -1, 1, 0, 3 },
	{ 2, -1, 1, 3, 0 },
	{ 2, -1, 3, 0, 1 },
	{ 2, -1, 3, 1, 0 },
	{ 2, 0, -1, 1, 3 },
	{ 2, 0, -1, 3, 1 },
	{ 2, 0, 1, -1, 3 },
	{ 2, 0, 1, 3, -1 },
	{ 2, 0, 3, -1, 1 },
	{ 2, 0, 3, 1, -1 },
	{ 2, 1, -1, 0, 3 },
	{ 2, 1, -1, 3, 0 },
	{ 2, 1, 0, -1, 3 },
	{ 2, 1, 0, 3, -1 },
	{ 2, 1, 3, -1, 0 },
	{ 2, 1, 3, 0, -1 },
	{ 2, 3, -1, 0, 1 },
	{ 2, 3, -1, 1, 0 },
	{ 2, 3, 0, -1, 1 },
	{ 2, 3, 0, 1, -1 },
	{ 2, 3, 1, -1, 0 },
	{ 2, 3, 1, 0, -1 },
	{ 3, -1, 0, 1, 2 },
	{ 3, -1, 0, 2, 1 },
	{ 3, -1, 1, 0, 2 },
	{ 3, -1, 1, 2, 0 },
	{ 3, -1, 2, 0, 1 },
	{ 3, -1, 2, 1, 0 },
	{ 3, 0, -1, 1, 2 },
	{ 3, 0, -1, 2, 1 },
	{ 3, 0, 1, -1, 2 },
	{ 3, 0, 1, 2, -1 },
	{ 3, 0, 2, -1, 1 },
	{ 3, 0, 2, 1, -1 },
	{ 3, 1, -1, 0, 2 },
	{ 3, 1, -1, 2, 0 },
	{ 3, 1, 0, -1, 2 },
	{ 3, 1, 0, 2, -1 },
	{ 3, 1, 2, -1, 0 },
	{ 3, 1, 2, 0, -1 },
	{ 3, 2, -1, 0, 1 },
	{ 3, 2, -1, 1, 0 },
	{ 3, 2, 0, -1, 1 },
	{ 3, 2, 0, 1, -1 },
	{ 3, 2, 1, -1, 0 },
	{ 3, 2, 1, 0, -1 },
};

extern uint32_t rngstate;

inline uint32_t fastrng() {
	rngstate = rngstate * 1664525 + 1013904223; // numerical recipes
	return rngstate;
}
inline double fastrngd() {
	return (double)fastrng() / (double)std::numeric_limits<uint32_t>::max();
}

inline const int * getdirperm() {
	return DIRPERMUTATIONS[fastrng() % 24];
}

inline const int * getoptdirperm() {
	return OPTDIRPERMUTATIONS[fastrng() % 120];
}

// some large prime, to obtain pseudo-random permutations
const int PRIMES[] = {
	312007,
	312023,
	312029,
	312031,
	312043,
	312047,
	312071,
	312073,
	312083,
	312089,
	312101,
	312107,
	312121,
	312161,
	312197,
	312199,
	312203,
	312209,
	312211,
	312217,
	312229,
	312233,
	312241,
	312251,
	312253,
	312269,
	312281,
	312283,
	312289,
	312311,
	312313,
	312331,
	312343,
	312349,
	312353,
	312371,
	312383,
	312397,
	312401,
	312407,
	312413,
};

inline int getprime() {
	return PRIMES[fastrng() % (sizeof(PRIMES) / sizeof(PRIMES[0]))];
}


/*
    struct to store current state information
*/
struct State
{
	/*
	Variables
	*/
	int rows, cols,
	turn, turns,
	noPlayers;
	double attackradius, spawnradius, viewradius;
	double loadtime, turntime;
	std::vector<double> scores;
	bool gameover;
	int64_t seed;
	bool newwater; // have we seen new water this turn?
	bool newsquare; // have we seen a square for the first time this turn?

	std::vector<std::vector<Square> > grid;
	std::vector<Location> myAnts, enemyAnts, myHills, enemyHills, food;

	Timer timer;
	Bug bug;

	/*
	Functions
	*/
	State();
	~State();

	void setup();
	void reset();

	void makeMove(const Location &loc, int direction);

	uint manhattanDistance(const Location & loc1, const Location & loc2) const;
	double distance(const Location &loc1, const Location &loc2) const;
	Location getLocation(const Location &startLoc, int direction) const;
	uint eucliddist2(const Location & loc1, const Location & loc2) const;
	Location addLocations(const Location & a, const Location & b) const {
		return Location((a.row + b.row + rows) % rows, (a.col + b.col + cols) % cols);
	}

	void updateVisionInformation();
};

std::ostream& operator<<(std::ostream &os, const State &state);
std::istream& operator>>(std::istream &is, State &state);

#endif //STATE_H_
