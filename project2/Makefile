CC=g++
CFLAGS=-std=c++11 
SOURCES=molecule.cpp 2.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=2

default:
	g++ -I /Users/ruheedcunha/Desktop/playground/eigen -std=c++11 -c molecule.cpp
	g++ -I /Users/ruheedcunha/Desktop/playground/eigen -std=c++11 -c 2.cpp
	g++ -I /Users/ruheedcunha/Desktop/playground/eigen -std=c++11 molecule.o 2.o -o 2
