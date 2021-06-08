LDFLAGS += -lm -std=c++17
CXXFLAGS += -O2 -Wall -Wextra -Werror

all: stattest

clean:
	$(RM) stattest
