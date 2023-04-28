CC := gcc
TARGET := arrangement

SRC := arrangement.c bitset.c search.c labelling.c
OBJ := $(SRC:.c=.o)
DEP := $(OBJ:%.o=%.d)

UNAME := $(shell uname)
LOCAL := 1
DEBUG := 0
FLTO := 1

COMMONFLAGS := #-Wall -Werror -Wextra
CFLAGS := $(COMMONFLAGS) -std=gnu99  -gdwarf-3 -fopt-info-vec-optimized -ftree-vectorize -mavx512f
LDFLAGS := $(COMMONFLAGS) -lm

ifeq ($(PARALLEL),1)
	CFLAGS += -DPARALLEL -fopencilk
	LDFLAGS += -fopencilk
endif

ifeq ($(UNAME),Darwin)
	CFLAGS += -DMACPORT -D__MACH__
endif

# for Cilksan and Cilkscale
ifeq ($(CILKSAN),1)
	CFLAGS += -fsanitize=cilk
	LDFLAGS += -fsanitize=cilk
else ifeq ($(CILKSCALE),1)
	CFLAGS += -fcilktool=cilkscale
	LDFLAGS += -fcilktool=cilkscale
endif

ifeq ($(DEBUG),1)
	CFLAGS += -O0 -DDEBUG -fno-inline
else
	CFLAGS += -O3 -DNDEBUG
endif

ifeq ($(FLTO),1)
	ifneq ($(UNAME),Darwin)
		CFLAGS += -flto
		LDFLAGS += -flto
	endif
endif

# DO NOT MODIFY THE ARCHITECTURES
ifeq ($(LOCAL),0)
	# You will be graded on haswell.
	CFLAGS += -march=haswell
else
	CFLAGS += -march=x86-64
endif

# First target is default
default: print_local_warning $(TARGET) ;

ifeq ($(LOCAL),0)
print_local_warning: ;
else
print_local_warning:
	@printf "\033[01;33mBE ADVISED: You have selected to build for your native architecture. This might be different than Haswell, which the AWS grading machines use.\033[00m\n"
	@printf "\033[01;33mUnless you have a good reason, unset the LOCAL flag or use LOCAL=0 before getting performance numbers on AWS.\033[00m\n"
endif

# Link the leiserchess binary
$(TARGET): $(OBJ)
	$(CC) $(PARAMS) $(LDFLAGS) $^ -o $@

%.o: %.c Makefile
	$(CC) $(PARAMS) $(CFLAGS) -MMD -c $< -o $@

-include $(DEP)

.PHONY: default format clean print_local_warning
clean:
	rm -f *.o *.d* *~ $(TARGET)

format:
	clang-format -i --style=file *.c *.h