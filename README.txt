This project is a submission to the 2011 Google AI Challenge.
Visit http://aichallenge.org/ for more information on the challenge.

To see the bot in action, get the tools package either from the contest
homepage or from the repository at https://github.com/aichallenge/aichallenge/

You can compile the submission using the Makefile, or via CMake using the
provided CMakeLists.txt.

The repository contains a number of tags following the development history of
the submission. These are my notes of the progress corresponding to the tags:

v2: search for food, expand zoc
v4: + tactical
v4-1: + scouting
v4-2: + hill defense
v4-3: + opportunistic attack
v4-5: + organized offense
v4-6: + slightly reduce number of offense ants, use symmetry finder for target selection, and some cleanups
v4-8: crash fix
v4-9: + various tweaks to tactical
v4x10: disable offense
v5 = v4x10
v5.01: tactical rework
v5.02: enable tactical parameterization
v5.03: add tactical aggressive mode
v5.04: add diffusion based movement
v5.05: limit computations via timer
v5.06: add TacticalSms (submap sampling)
v5.07: TacticalSample: a1k0n-style sampling
v5.08: use TacticalSm, with re-added valuation of food
v5.09: use TacticalSm, with max-avg instead of max-min (and improved duplicate move detection)
v6 = v5.08
v5.10: use TacticalSm, with learning which move selection to use
v5.11: higher aggressiveness
v5.12: add some more MoveSels
v5.13: disable MtMt movesels
v5.14 = v5.10
v5.15: remove old tactical code
v5.16: disable assertions
v10 = v5.16

My final submission (version 10 on the contest homepage) is equal to the tag v5.16.

A blog entry describing some of the thoughts behind the bot can be found here:
http://nhaehnle.blogspot.com/2011/12/ai-challenge-look-back.html
