#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "module.h"

struct Bot;
struct State;

struct Diffusion : Module {
	Diffusion(Bot & b);
	virtual ~Diffusion();

	virtual void init();
	virtual void run();

	void diffuse();

	struct Data;

	Bot & bot;
	State & state;
	Data & d;
};

#endif // DIFFUSION_H
