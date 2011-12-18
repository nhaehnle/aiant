
set (CMAKE_CXX_FLAGS_DEBUG "-g -O3 -DDEBUG -Wall" CACHE STRING "Set" FORCE)
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wall" CACHE STRING "Set" FORCE)

add_executable(MyBot
	MyBot.cc
	Bot.cc
	State.cc
	tactical_sm.cc
	tactical_smbase.cc
	diffusion.cc
	zoc.cc
	scout.cc
	hilldefense.cc
	foodseeker.cc
	opportunisticattack.cc
	offense.cc
	symmetry.cc
)
