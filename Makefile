TARGET=Bi_cutter.exe
all: Bi_cutter.exe

$(TARGET) : bi_cutter.C
	#g++ -c `root-config --cflags --libs` -lHistPainter
	#g++ -c `root-config --cflags --libs`  sndisplay.cc
	g++ -c `root-config --cflags --libs` bi_cutter.C
	#g++ -o $(TARGET) sndisplay.o test.o myDictionary.o `root-config --cflags --libs`
	g++ -o $(TARGET) bi_cutter.o myDictionary.o `root-config --cflags --libs`
clean :
	rm $(TARGET)

