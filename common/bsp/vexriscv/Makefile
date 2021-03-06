PLATFORM ?= murax
RISCV_ARCH ?= rv32im
RISCV_ABI ?= ilp32
RISCV_CMODEL ?= medany
RISCV_ARCHFLAGS += -march=$(RISCV_ARCH)
RISCV_ARCHFLAGS += -mabi=$(RISCV_ABI)
RISCV_ARCHFLAGS += -mcmodel=$(RISCV_CMODEL)

CFLAGS += -Os -Wall -Wextra
CFLAGS += $(RISCV_ARCHFLAGS)
CFLAGS += -fstrict-volatile-bitfields
CFLAGS += --specs=nosys.specs

ifeq ($(PLATFORM),murax)
CFLAGS += -DVEXRISCV_VOLATILE
endif

ifeq ($(PLATFORM),pqvexriscvup5k)
CFLAGS += -DVEXRISCV_VOLATILE -DVEXRISCV_RWMTVEC
endif

ifeq ($(PLATFORM),pqvexriscvicoboard)
CFLAGS += -DVEXRISCV_VOLATILE -DVEXRISCV_RWMTVEC
endif

ifeq ($(PLATFORM),pqvexriscvsim)
CFLAGS += -DVEXRISCV_RWMTVEC
endif

CROSS_PREFIX ?= riscv64-unknown-elf
CC = $(CROSS_PREFIX)-gcc
AR = $(CROSS_PREFIX)-ar

SRCS = init.c start.s
OBJS = $(PLATFORM)_init.o $(PLATFORM)_start.o
TARGET = lib$(PLATFORM)bsp.a

all: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) rcs $@ $^

$(PLATFORM)_%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(PLATFORM)_%.o: %.S
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJS) $(TARGET)
