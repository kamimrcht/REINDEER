SRC:=update_reindeer_info.cpp
BIN:=$(SRC:.cpp=)

all: $(BIN)

$(BIN): $(SRC)

# do test
test:
	cd test && $(MAKE) all

.PHONY: test
