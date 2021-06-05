LDFLAGS += -lm -std=c++17
CXXFLAGS += -O2

all: stattest

clean:
	$(RM) stattest
