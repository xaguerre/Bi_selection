TARGET=Bi_cutter.exe
all: Bi_cutter.exe

$(TARGET) : bi_cutter.C
	g++ -c `root-config --cflags --libs` bi_cutter.C
	g++ -o $(TARGET) bi_cutter.o myDictionary.o `root-config --cflags --libs`
clean :
	rm $(TARGET)
