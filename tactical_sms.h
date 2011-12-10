#ifndef TACTICAL_SMS_H
#define TACTICAL_SMS_H

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

struct TacticalSms : TacticalSmBase {
	struct Data;
	struct ShadowAnt;
	struct Theater;

	TacticalSms(Bot & bot_);
	~TacticalSms();

	void init();
	void run();

	bool timeover();
	void generate_theater(const Location & center);
	void run_theater(uint theateridx);
	void push_moves(uint theateridx, uint myidx);
	void pull_moves(uint theateridx);

	void evaluate_moves(Theater & th, PlayerMove & pm, PlayerMove & enemymove, float & myvalue, float & enemyvalue);
	void evaluate_new_moves(uint theateridx);
	bool get_improve_pair(const std::vector<PlayerMove *> & moves, uint & myidx, uint & enemyidx);
	void improve(uint theateridx, uint myidx, uint enemyidx);

	Data & d;
};

#endif // TACTICAL_SMS_H
