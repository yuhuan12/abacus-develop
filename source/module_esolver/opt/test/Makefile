
#==========================
# Compiler information 
#==========================
CPLUSPULS = g++
# CPLUSPULS = icpc
OBJ_DIR = opt_obj
BIN_DIR = testbin
NP=8
#==========================
# Objects
#==========================
VPATH=../../../module_base\
:../

OPT_OBJS_0=global_function.o\
tool_quit.o\
global_file.o\
memory.o\
global_variable.o\
timer.o\
opt_DCsrch.o\
opt_CG.o\
test_tools.o\


FILE=CG_test.exe\
TN_test.exe

#==========================
# Options
#==========================

OPT_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${OPT_OBJS_0})
TESTFILE=$(patsubst %.exe, ${BIN_DIR}/%.exe, ${FILE})

OPTS = -Ofast -std=c++11 -m64 -Wall -pedantic -g

LIBS = -lgtest -lpthread -lgtest_main

#==========================
# MAKING OPTIONS
#==========================

opt :
	@make init
	@make -j ${NP} ${OPT_OBJS}
	@make -j ${NP} ${TESTFILE}

init :
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@ if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	@ if [ ! -d $(OBJ_DIR)/README ]; then echo "This directory contains all of the .o files" > $(OBJ_DIR)/README; fi

${BIN_DIR}/%.exe: %.cpp ${OPT_OBJS}
	${CPLUSPULS} ${OPTS} $< ${OPT_OBJS} ${LIBS} -o $@

${OBJ_DIR}/%.o:%.cpp
	${CPLUSPULS} ${OPTS} -c $< -o $@

.PHONY: clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@ if [ -d $(BIN_DIR) ]; then rm -rf $(BIN_DIR); fi