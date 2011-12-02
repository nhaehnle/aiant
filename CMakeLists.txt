
set (CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG -Wall" CACHE STRING "Set" FORCE)
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wall" CACHE STRING "Set" FORCE)

add_executable(MyBot MyBot.cc Bot.cc State.cc zoc.cc tactical.cc scout.cc hilldefense.cc foodseeker.cc)
