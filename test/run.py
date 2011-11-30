#!/usr/bin/env python

import argparse
import json
import os
import random
import stat
import subprocess


def listmaps(dirname):
	maps = []
	for name in os.listdir(dirname):
		path = dirname + '/' + name
		st = os.stat(path)
		if stat.S_ISDIR(st.st_mode):
			maps += listmaps(path)
		elif name.endswith(".map"):
			rec = {"path": path}
			with open(path, 'r') as filp:
				for line in filp:
					parts = line.split()
					if not parts:
						continue
					if parts[0] == "players":
						rec["players"] = int(parts[1])
			if not "players" in rec:
				print "Map {0} has no players info".format(path)
				continue

			print "Found {rec[players]} player map {rec[path]}".format(**locals())
			maps.append(rec)
	return maps

def getfirstgameid():
	gameid = 0
	for name in os.listdir("game_logs"):
		id = int(name.split('.')[0])
		gameid = max(gameid, id+1)
	return gameid

def main():
	with open('bots.json', 'r') as filp:
		settings = json.load(filp)

	maps = listmaps(settings["ants"] + "maps")
	allbots = [ (name, name[name.rfind('/')+1:]) for name in settings["bots"] ]
	gameid = getfirstgameid()

	while True:
		maprec = random.choice(maps)
		print "Game {gameid}: {maprec[players]} player map {maprec[path]}".format(**locals())

		mapbots = [random.choice(allbots) for i in range(maprec["players"])]
		for i in range(len(mapbots)):
			print "  {i}: {bot}".format(i=i+1, bot=mapbots[i][1])

		args = [
			settings["ants"] + "playgame.py",
			"--turns", "1000",
			"--log_dir", "game_logs",
			"--log_input", "--nolaunch", "--serial",
			"-g", str(gameid),
			"--map", maprec["path"],
		] + [bot[0] for bot in mapbots]
		subprocess.check_call(args)

		with open('game_logs/{0}.replay'.format(gameid), 'r') as filp:
			replay = json.load(filp)
			if "error" in replay:
				print "An error occurred"
				print replay["error"]
				return

			score = replay["score"]
			status = replay["status"]

			print "Outcome after {0} turns".format(max(replay["playerturns"]))
			for i in range(len(mapbots)):
				print "  {i}: {name}  score: {score}  status: {status}".format(
					i=i+1,
					name=replay["playernames"][i],
					score=score[i],
					status=status[i])

		gameid += 1

if __name__ == "__main__":
	main()
