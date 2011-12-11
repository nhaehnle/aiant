#ifndef TACTICAL_SAMPLE_H
#define TACTICAL_SAMPLE_H

#include <cstdlib>
#include <vector>

#include "Location.h"
#include "module.h"

struct Bot;
struct State;

struct TacticalSample : Module {
	struct Data;

	TacticalSample(Bot & b);
	virtual ~TacticalSample();

	virtual void init();
	virtual void run();
	void run_samplers();
	void do_sample_markonly(uint sampleridx);

	void init_map();
	void init_ants();
	void mark_hills(uint who, const std::vector<Location> & hills);

	Bot & bot;
	State & state;
	Data & d;
};

#endif // TACTICAL_SAMPLE_H
