CXX = g++-8.1
CXXFLAGS = -std=c++17 -O2 -Wall
OBJS = Bitmap.o main.o
TARGET = main

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(CXXFLAGS)

main.o: main.cpp Bitmap.o
	$(CXX) -c main.cpp $(CXXFLAGS)

Bitmap.o: Bitmap.cpp Bitmap.h bmp.h 
	$(CXX) -c Bitmap.cpp $(CXXFLAGS)

clean:
	rm $(OBJS) $(TARGET)