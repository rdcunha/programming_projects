CC=g++
CFLAGS=-std=c++11 
SOURCES=hf.cpp 5.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=5

default:
	g++ -I /Users/ruheedcunha/Desktop/playground/eigen -std=c++11 -c hf.cpp
	g++ -I /Users/ruheedcunha/Desktop/playground/eigen -std=c++11 -c cc.cpp
	g++ -I /Users/ruheedcunha/Desktop/playground/eigen -std=c++11 -c 5.cpp
	g++ -I /Users/ruheedcunha/Desktop/playground/eigen -std=c++11 hf.o cc.o 5.o -o 5
