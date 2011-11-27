#include "tactical.h"

#include <cstring>

#include "Bot.h"
#include "State.h"

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
//@}

struct Submap {
	static const int Radius = 8;
	static const int Size = 2 * Radius + 1;

	static const unsigned char Water = 0x20;
	static const unsigned char Food = 0x10;
	static const unsigned char Hill = 0x08;
	static const unsigned char Ant = 0x04;
	static const unsigned char Mine = 0x02;
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
	for (cur.row = 0; cur.row < Submap::Size; ++cur.row) {
		for (cur.col = 0; cur.col < Submap::Size; ++cur.col) {
			int fieldvalue = 0;
			unsigned char field = sm[cur];

			if (!(field & Submap::Water))
				fieldvalue += ValueField;
			if (field & Submap::Food)
				fieldvalue += ValueFood;
			if (field & Submap::Ant)
				fieldvalue += potential[cur];

			if (field & Submap::Enemy)
				value -= fieldvalue;
			else if (field & Submap::Mine)
				value += fieldvalue;
		}
	}

	return value;
}

void Tactical::make_moves(const Location & center)
{
	Submap base_sm;
	gensubmap(base_sm, center);
	int base_value = evaluate(base_sm);

	state.bug << "Tactical around " << center <<", base value " << base_value << endl;
	state.bug << base_sm;
}
