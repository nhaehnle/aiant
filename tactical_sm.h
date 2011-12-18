#ifndef TACTICAL_SM_H
#define TACTICAL_SM_H

#include <stdlib.h>
#include <vector>

#include "tactical_smbase.h"

struct Bot;
struct Location;
struct PlayerMove;
template<typename T>
struct BaseSubmap;
struct Submap;
struct State;

struct TacticalSm : TacticalSmBase {
	struct Data;
	struct ShadowAnt;
	struct Theater;

	TacticalSm(Bot & bot_);
	~TacticalSm();

	void init();
	void run();
	void learn();

	bool timeover();
	void generate_theater(const Location & center);
	void run_theater(uint theateridx);
	void push_moves(uint theateridx, uint myidx);
	void pull_moves(uint theateridx);
	int pull_enemy_moves(uint theateridx);

	void evaluate_moves(Theater & th, PlayerMove & pm, PlayerMove & enemymove, float & myvalue, float & enemyvalue);
	void evaluate_new_moves(uint theateridx);
	bool get_improve_pair(const std::vector<PlayerMove *> & moves, uint & myidx, uint & enemyidx);
	void improve(uint theateridx, uint myidx, uint enemyidx);

	void choose_strategy(uint theateridx);
	uint choose_max_min_move(uint theateridx);
	uint choose_max_avg_move(uint theateridx);

	Data & d;
};

#endif // TACTICAL_SM_H
