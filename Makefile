
SRC 	:= ./src
OBJ 	:= ./obj
INC 	:= ./include
UTIL	:= $(SRC)/behavior_code
BIN		:= ./bin

SOURCE 		:= $(wildcard $(SRC)/*.cc)
OBJS 		:= $(patsubst %.cc, $(OBJ)/%.o, $(notdir $(SOURCE)))
MAIN_SOURCE	:= $(wildcard $(UTIL)/*.cc)
MAIN_OBJS	:= $(patsubst %.cc, $(OBJ)/%.o, $(notdir $(MAIN_SOURCE)))

CC = g++
CFLAGS = -g -Wall

LIBS	:= -lntl -lgmp -lm -pthread
INCLUDE	:= -I $(INC)
TARGET	:= $(BIN)/main

all: main #$(TARGET)

$(MAIN_OBJS):$(OBJ)/%.o: $(UTIL)/%.cc
	$(CC) -c $< -o $@ $(INCLUDE)

$(OBJS):$(OBJ)/%.o: $(SRC)/%.cc
	$(CC) -c $< -o $@ $(INCLUDE)


$(TARGET): $(MAIN_OBJS) $(OBJS) 
	$(CC) -o $@ $(MAIN_OBJS) $(OBJS) $(LIBS)

main: $(TARGET)
	$(BIN)/main 0 8191 22 16 200

barrett: $(OBJS) $(MAIN_OBJS) 
	$(CC) -o $(BIN)/barrett.exe $(OBJS) $(MAIN_OBJS) $(LIBS)
	$(BIN)/barrett.exe

AE: $(OBJS) $(MAIN_OBJS) 
	$(CC) -o $(BIN)/AE.exe $(OBJS) $(MAIN_OBJS) $(LIBS)
	$(BIN)/AE.exe

.PHONY: clean


clean:
	rm -rf $(TARGET)
	rm -rf $(OBJ)/*.o
	rm -rf $(BIN)/*.exe
	rm -rf $(BIN)/main
	rm -rf *.txt
	rm -rf test_input/*.txt
	rm -rf ROM_Data/*.txt
	rm -rf my_print_data/*.txt
	rm -rf my_print_data/v_code_gen/*.txt
	rm -rf SPMB_tw/*.txt
	rm -rf ROM_Data/ROM0/*.txt
	rm -rf ROM_Data/ROM1/*.txt
	rm -rf ROM_Data/ROM2/*.txt
	rm -rf NWC_PrintData/*/*.txt
	