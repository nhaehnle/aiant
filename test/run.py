#!/usr/bin/env python
#
# Needs settings in file bots.json, which should look this this:
# {
#	"bots": [
#		"python /home/prefect/gitdata/aichallenge/ants/dist/sample_bots/python/GreedyBot.py",
#		"python /home/prefect/gitdata/aichallenge/ants/dist/sample_bots/python/HoldBot.py",
#		"python /home/prefect/gitdata/aichallenge/ants/dist/sample_bots/python/HunterBot.py",
#		"python /home/prefect/gitdata/aichallenge/ants/dist/sample_bots/python/LeftyBot.py",
#		"/home/prefect/data/dev/aichallenge/v2/MyBot-v2",
#		"/home/prefect/data/dev/aichallenge/v4/MyBot-v4"
#	],
#	"ants": "/home/prefect/gitdata/aichallenge/ants/"
# }
# "ants" is the path to the game distribution (contains the playgame.py file)

import json
import os
import random
import subprocess

import helpers

def getfirstgameid():
	gameid = 0
	for name in os.listdir("game_logs"):
		id = int(name.split('.')[0])
		gameid = max(gameid, id+1)
	return gameid

def main():
	with open('bots.json', 'r') as filp:
		settings = json.load(filp)

	maps = [maprec for maprec in helpers.listmaps(settings["ants"] + "maps") if maprec["players"] <= len(settings["bots"])]
	allbots = [ (name, name[name.rfind('/')+1:]) for name in settings["bots"] ]
	gameid = getfirstgameid()

	while True:
		maprec = random.choice(maps)
		print "Game {gameid}: {maprec[players]} player map {maprec[path]}".format(**locals())

		mapbots = random.sample(allbots, maprec["players"])

		for i in range(len(mapbots)):
			print "  {i}: {bot}".format(i=i+1, bot=mapbots[i][1])

		args = [
			settings["ants"] + "playgame.py",
			"--turns", "1000",
			"--turntime", "500",
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
