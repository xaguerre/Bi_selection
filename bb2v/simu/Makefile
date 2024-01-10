TARGET=Bb_cutter.exe
all: Bb_cutter.exe

$(TARGET) : bb_cutter_simu.C
	g++ -c `root-config --cflags --libs` bb_cutter_simu.C
	g++ -o $(TARGET) bb_cutter_simu.o myDictionary.o `root-config --cflags --libs`
clean :
	rm $(TARGET)
