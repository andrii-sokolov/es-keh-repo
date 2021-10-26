CC = g++ 
CFLAGS = -Wall -O2
TARGET = eskeh
SOURCE = src/main.cc 
LDFLAGS = -lboost_system -lboost_iostreams -lboost_filesystem

all: $(SOURCE)
	$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) $(LDFLAGS)   

execute: $(TARGET)
	cp Data/Experimental_Data/*.csv ./;
	./$(TARGET);
	
clear: 
	rm -f *.csv;
	rm -f $(TARGET);
	rm -f plotfiles.py;
	rm -f *.txt;

plot: 
	cp pysrc/plotfiles.py ./
	python plotfiles.py 
