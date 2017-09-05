CXXFLAGS=--std=gnu++0x -O2 -fno-math-errno -march=native -g -Wall
LDLIBS=-lavcodec -lavformat -lavutil -lswresample

all: synth decode sync cleaner

%.o: %.cpp
	$(CXX) -MMD -MP $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

OBJS=audioreader.o decode.o synth.o synth_main.o interpolate.o sync.o level.o filter.o

DEPS=$(OBJS:.o=.d)
-include $(DEPS)

decode: interpolate.o audioreader.o decode.o level.o filter.o
	$(CXX) -o $@ $^ $(LDLIBS) $(LDFLAGS)

synth: synth.o synth_main.o
	$(CXX) -o $@ $^ $(LDFLAGS)

sync: interpolate.o sync.o
	$(CXX) -o $@ $^ $(LDFLAGS)

cleaner: cleaner.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) synth decode sync cleaner $(OBJS) $(DEPS)
