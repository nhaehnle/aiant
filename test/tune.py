#!/usr/bin/env python
#
# Tune bot's parameters using metaheuristics
# Call as ./tune.py tunesettings.json

import json
import math
import os
import random
import subprocess
import sys

import helpers
import trueskill

DEFAULT_SIGMASCHEDULE = [1.0 / 3.0, 1.0 / 3.0] + [1.0/(n*n) for n in range(2,10)]
MINGAMES = 15

def weightedsample(population, count, key):
	assert count <= len(population)

	weights = [key(p) for p in population]
	totalweight = sum(weights)

	result = []
	while len(result) < count:
		f = random.uniform(0, totalweight)
		for idx in range(len(weights)):
			if weights[idx] > 0.0:
				f -= weights[idx]
				if f < 0.0:
					break

		assert weights[idx] >= 0.0
		result.append(population[idx])
		totalweight -= weights[idx]
		weights[idx] = -1.0
		while weights[-1] < 0.0:
			weights.pop()

	return result

class Bot(object):
	def __init__(self):
		self.skill = (trueskill.INITIAL_MU, trueskill.INITIAL_SIGMA)
		self.nrgames = 0

class StaticBot(Bot):
	def __init__(self, cmd):
		super(StaticBot, self).__init__()
		self.cmd = cmd

	def getcmd(self):
		return self.cmd

	def getname(self):
		return self.cmd[self.cmd.rfind('/')+1:]

class TuneBot(Bot):
	def __init__(self, group, parameters):
		super(TuneBot, self).__init__()
		self.group = group
		self.parameters = parameters.copy()
		mykeys = set(parameters.keys())
		normkeys = set(group.settings["parameters"].keys())
		diff = mykeys.symmetric_difference(normkeys)
		if len(diff) > 0:
			print "TuneBot: symmetric difference in parameters:", diff
			raise ValueError("bad keys")

	def getcmd(self):
		return self.group.settings["bot"] + " " + \
			" ".join(["--{0} {1}".format(key, value) for key, value in self.parameters.iteritems()])

	def getname(self):
		bot = self.group.settings["bot"]
		return bot[bot.rfind('/')+1:]

class MutationSchedule(object):
	def __init__(self, keys):
		self.keys = keys
		self.sigmaschedule = DEFAULT_SIGMASCHEDULE
		self.round = 0

	def nextmutation(self):
		sigma = self.sigmaschedule[self.round % len(self.sigmaschedule)]
		self.round += 1
		return [(key, sigma) for key in self.keys]

class TuneGroup(object):
	def __init__(self, settings):
		self.settings = settings
		self.schedule = MutationSchedule(settings["parameters"].keys())

		print "Initializing TuneGroup({0})".format(self.getname())
		try:
			with open(settings["population"], 'r') as filp:
				pop = json.load(filp)
				self.population = [TuneBot(self, ind) for ind in pop]
			print "  populated from file", settings["population"]
		except:
			print "  starting with a single seed individual"
			self.population = [
				TuneBot(self, dict([(key, settings["parameters"][key]["seed"]) for key in settings["parameters"].iterkeys()]))
			]
		self.generate_mutations()

	def getname(self):
		bot = self.settings["bot"]
		return bot[bot.rfind('/')+1:]

	def mutate(self, key, sigma, oldv):
		param = self.settings["parameters"][key]
		kind = param["kind"]
		if kind == "log":
			oldlog = math.log(oldv)
			s = abs(oldlog)
			if "minsigma" in param:
				s = max(s, param["minsigma"])
			else:
				s += 0.00001
			newlog = random.gauss(oldlog, s * sigma)
			return math.exp(newlog)
		elif kind == "linear":
			s = abs(oldv)
			if "minsigma" in param:
				s = max(s, param["minsigma"])
			else:
				s += 0.00001
			return random.gauss(oldv, s * sigma)
		elif kind == "lognormal":
			oldlog = math.log(oldv)
			newlog = random.gauss(oldlog, param["sigma"] * sigma)
			return math.exp(newlog)
		else:
			raise ValueError("bad key kind '{0}' for key '{1}'".format(kind, key))

	def generate_mutations(self):
		initialcount = len(self.population)
		keysigmas = self.schedule.nextmutation()

		idx = 0
		while len(self.population) < self.settings["popsize"]:
			newbot = TuneBot(self, self.population[idx].parameters)
			for key,sigma in keysigmas:
				newbot.parameters[key] = self.mutate(key, sigma, newbot.parameters[key])
			self.population.append(newbot)
			idx = (idx + 1) % initialcount

		self.save_population()

		print "After mutation", keysigmas
		self.print_population()

	def save_population(self):
		with open(self.settings["population"], 'w') as filp:
			json.dump([bot.parameters for bot in self.population], filp, indent=2)

	def print_population(self):
		for key in self.settings["parameters"].iterkeys():
			print "{0:>10}".format(key),
		print "    MU  SIGMA"
		for bot in self.population:
			for key in self.settings["parameters"].iterkeys():
				align = max(10, len(key))
				print "{value:>{align}.{width}}".format(value=bot.parameters[key], align=align, width=align-3),
			print "{mu:>6} {sigma:>6}".format(mu=bot.skill[0], sigma=bot.skill[1])

	def run(self):
		for bot in self.population:
			if bot.nrgames < MINGAMES:
				return bot

		print "Stabilization point reached, population is:"
		self.print_population()

		minskill = min(bot.skill[0] for bot in self.population)
		self.population = weightedsample(self.population, self.settings["popsize"] / 3, key=lambda bot: bot.skill[0] - minskill)

		print "After population culling:"
		self.print_population()

		self.generate_mutations()

		return self.population[-1]

def getfirstgameid():
	gameid = 0
	for name in os.listdir("tune_logs"):
		id = int(name.split('.')[0])
		gameid = max(gameid, id+1)
	return gameid

def main():
	global MINGAMES

	trueskill.SetParameters(beta=25.0/3.0, draw_probability=0.3)

	with open(sys.argv[1], 'r') as filp:
		settings = json.load(filp)

	if "mingames" in settings:
		MINGAMES = settings["mingames"]

	maps = helpers.listmaps(settings["ants"] + "maps")
	gameid = getfirstgameid()

	staticbots = [StaticBot(cmd) for cmd in settings["static"]]
	tunegroups = [TuneGroup(group) for group in settings["tune"]]

	while True:
		seedbot = tunegroups[0].run()
		otherbots = staticbots + [bot for group in tunegroups for bot in group.population if bot != seedbot]

		maprec = random.choice(maps)
		print "Game {gameid}: {maprec[players]} player map {maprec[path]}".format(**locals())

		mapbots = [seedbot] + random.sample(otherbots, maprec["players"] - 1)

		botcmds = [bot.getcmd() for bot in mapbots]
		for idx in range(len(botcmds)):
			cmd = botcmds[idx]
			skill = mapbots[idx].skill
			print "  ", "({mu:3}, {sigma:3}):".format(mu=skill[0], sigma=skill[1]), cmd[cmd.rfind('/')+1:]

		args = [
			settings["ants"] + "playgame.py",
			"--turns", "1000",
			"--turntime", "500",
			"--log_dir", "tune_logs",
			"--nolaunch", "--serial",
			"-g", str(gameid),
			"--map", maprec["path"],
		] + botcmds
		subprocess.check_call(args)

		with open('tune_logs/{0}.replay'.format(gameid), 'r') as filp:
			replay = json.load(filp)
			if "error" in replay:
				print "An error occurred"
				print replay["error"]
				return

			rank = replay["rank"]
			print "   Ranks:", ', '.join([str(r) for r in rank])
			for idx in range(len(rank)):
				mapbots[idx].rank = rank[idx]
				mapbots[idx].nrgames += 1

		trueskill.AdjustPlayers(mapbots)
		print "  ", ', '.join("{name} ({mu:3}, {sigma:3})".format(name=bot.getname(), mu=bot.skill[0], sigma=bot.skill[1]) for bot in mapbots)

		gameid += 1

if __name__ == "__main__":
	main()
