#ifndef BOT_H_
#define BOT_H_

#include "State.h"
#include "map.h"

struct HillDefense;
struct Scout;
struct Zoc;

struct Food {
	Location where;
	bool claimed;

	Food() : claimed(false) {}
};

struct PointOfInterest {
	Location where;
	uint distance; ///< Manhattan distance
	int direction;
};

struct Ant {
	Location where;

	bool hastactical;

	int direction;

	Ant() : hastactical(false), direction(-1) {}
};

/*
    This struct represents your bot in the game of Ants
*/
struct Bot
{
	State state;

	Bot();
	~Bot();

	void playGame();    //plays a single game of Ants

	void makeMoves();   //makes moves for a single turn
	void endTurn();     //indicates to the engine that it has made its moves

	uint myantidx_at(const Location & pos);
	uint foodidx_at(const Location & pos);

	void assign_food();
	bool try_rotate_move(uint antidx, const Map<bool> & claims);
	void make_moves();

	Zoc & m_zoc;
	Scout & m_scout;
	HillDefense & m_hilldefense;
	std::vector<Ant> m_ants;
	std::vector<Food> m_foods;
};

#endif //BOT_H_
