CXX := g++
CXX_FLAGS := -ggdb -Wall -Wno-c++11-extensions

BIN := bin
SRC := src
INCLUDE := include
LIB := lib
LIBRARIES :=
EXECUTABLE := main

SOURCES := $(wildcard $(SRC)/*.cpp)
OBJECTS := $(patsubst $(SRC)/%.cpp,$(BIN)/%.o,$(SOURCES))

# Rule to build the executable from the object files
$(BIN)/$(EXECUTABLE): $(OBJECTS)
	@echo "Linking..."
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $(OBJECTS) -o $@ $(LIBRARIES)

# Rule to compile each .cpp file into a .o file and generate a .d dependency file
$(BIN)/%.o: $(SRC)/%.cpp
	@mkdir -p $(BIN)
	@echo "Compiling $<..."
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -c $< -o $@
	@$(CXX) -MM -MT $@ -I$(INCLUDE) $< > $(BIN)/$*.d

# Include dependency files (.d) for each object file, if they exist
-include $(OBJECTS:.o=.d)

# Phony target to clean the build directory
clean:
	@echo "Clearing..."
	-rm -rf $(BIN)/*

# Phony target to run the program
run: $(BIN)/$(EXECUTABLE)
	@clear
	@echo "Executing..."
	./$(BIN)/$(EXECUTABLE)

.PHONY: all clean run
