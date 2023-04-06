TARGET=Simu_Bi_cutter.exe
all: Simu_Bi_cutter.exe

$(TARGET) : simu_bi_cutter.C
	g++ -c `root-config --cflags --libs` simu_bi_cutter.C
	g++ -o $(TARGET) simu_bi_cutter.o myDictionary.o `root-config --cflags --libs`
clean :
	rm $(TARGET)
