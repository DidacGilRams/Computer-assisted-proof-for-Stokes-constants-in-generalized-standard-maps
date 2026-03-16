# ============================================================
# Project: Computer-assisted-proof-for-Stokes-constants...
# Location expected: repository_root/CODE/Makefile
# ============================================================

SHELL := /bin/bash

# ----- Project programs -----
PROGS ?= PI1_N4 PI1_N5 PI2_3_N9 TI1_N4
OTHERS ?=

# ----- Toolchain -----
CXX ?= g++
BUILD_TYPE ?= Release
NPROC ?= $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)

# ----- Paths -----
CURDIR_ABS  := $(abspath $(CURDIR))
PROJECT_ROOT := $(abspath $(CURDIR_ABS)/..)

DEPS_DIR       ?= $(PROJECT_ROOT)/.deps
CAPD_SRC_DIR   ?= $(DEPS_DIR)/CAPD
CAPD_BUILD_DIR ?= $(CAPD_SRC_DIR)/build
CAPD_PREFIX    ?= $(DEPS_DIR)/capd-install

OBJDIR ?= .obj
BINDIR ?= bin

# ----- CAPD source -----
CAPD_REPO ?= https://github.com/CAPDGroup/CAPD.git
CAPD_REF  ?= master
CAPD_CONFIG := $(CAPD_PREFIX)/bin/capd-config

# Optional extra CMake flags for CAPD, e.g.:
# make setup CAPD_CMAKE_FLAGS="-DCAPD_ENABLE_MULTIPRECISION=ON"
CAPD_CMAKE_FLAGS ?=

# ----- Flags obtained from installed CAPD -----
CAPD_CFLAGS = $(shell [ -x "$(CAPD_CONFIG)" ] && "$(CAPD_CONFIG)" --cflags)
CAPD_LIBS   = $(shell [ -x "$(CAPD_CONFIG)" ] && "$(CAPD_CONFIG)" --libs)

CXXFLAGS += -O2 -Wall $(CAPD_CFLAGS)
LDFLAGS  += -Wl,-rpath,$(CAPD_PREFIX)/lib -Wl,-rpath,$(CAPD_PREFIX)/lib64

# ----- Derived object lists -----
OTHERS_OBJ = $(OTHERS:%=$(OBJDIR)/%.o)
OBJ_FILES  = $(OTHERS_OBJ) $(PROGS:%=$(OBJDIR)/%.o)

.DEFAULT_GOAL := all

.PHONY: all setup build clean distclean rebuild help \
        check-tools capd-clone capd-fetch capd-configure capd-build capd-install \
        capd-clean capd-purge info

all: $(PROGS:%=$(BINDIR)/%)

build: all

help:
	@echo "Targets disponibles:"
	@echo "  make setup      -> clones, compiles and installs CAPD in $(CAPD_PREFIX)"
	@echo "  make            -> compile the programs of the project"
	@echo "  make clean      -> cleans binary and objets from the project"
	@echo "  make distclean  -> cleans the project + builds/installs local of CAPD"
	@echo "  make info       -> shows the active configuration"
	@echo ""
	@echo "Useful variables:"
	@echo "  CAPD_REF=<ref>                (by defaut: $(CAPD_REF))"
	@echo "  CAPD_CMAKE_FLAGS='<flags>'    (ej: -DCAPD_ENABLE_MULTIPRECISION=ON)"
	@echo "  BUILD_TYPE=<Debug|Release>"
	@echo "  CXX=<compiler>"

info:
	@echo "PROJECT_ROOT   = $(PROJECT_ROOT)"
	@echo "DEPS_DIR       = $(DEPS_DIR)"
	@echo "CAPD_SRC_DIR   = $(CAPD_SRC_DIR)"
	@echo "CAPD_BUILD_DIR = $(CAPD_BUILD_DIR)"
	@echo "CAPD_PREFIX    = $(CAPD_PREFIX)"
	@echo "CAPD_CONFIG    = $(CAPD_CONFIG)"
	@echo "CAPD_REF       = $(CAPD_REF)"
	@echo "CXX            = $(CXX)"
	@echo "BUILD_TYPE     = $(BUILD_TYPE)"
	@echo "PROGS          = $(PROGS)"

check-tools:
	@command -v git   >/dev/null 2>&1 || { echo "Error: git not found"; exit 1; }
	@command -v cmake >/dev/null 2>&1 || { echo "Error: cmake not found"; exit 1; }
	@command -v $(CXX) >/dev/null 2>&1 || { echo "Error: compiler $(CXX) not found"; exit 1; }
	@command -v make  >/dev/null 2>&1 || { echo "Error: make not found"; exit 1; }

setup: check-tools capd-install

capd-clone:
	@mkdir -p "$(DEPS_DIR)"
	@if [ ! -d "$(CAPD_SRC_DIR)/.git" ]; then \
		echo "[CAPD] clonaning repository in $(CAPD_SRC_DIR)"; \
		git clone "$(CAPD_REPO)" "$(CAPD_SRC_DIR)"; \
	else \
		echo "[CAPD] repository already in $(CAPD_SRC_DIR)"; \
	fi

capd-fetch: capd-clone
	@echo "[CAPD] upgreating refs and checkout $(CAPD_REF)"
	@git -C "$(CAPD_SRC_DIR)" fetch --all --tags --prune
	@git -C "$(CAPD_SRC_DIR)" checkout "$(CAPD_REF)"

capd-configure: capd-fetch
	@echo "[CAPD] configurating with CMake"
	@cmake -S "$(CAPD_SRC_DIR)" -B "$(CAPD_BUILD_DIR)" \
		-DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" \
		-DCMAKE_INSTALL_PREFIX="$(CAPD_PREFIX)" \
		-DCAPD_BUILD_TESTS=OFF \
		-DCAPD_BUILD_EXAMPLES=OFF \
		$(CAPD_CMAKE_FLAGS)

capd-build: capd-configure
	@echo "[CAPD] compiling"
	@$(MAKE) -C "$(CAPD_BUILD_DIR)" -j"$(NPROC)"

capd-install: capd-build
	@echo "[CAPD] installing in $(CAPD_PREFIX)"
	@$(MAKE) -C "$(CAPD_BUILD_DIR)" install
	@[ -x "$(CAPD_CONFIG)" ] || { echo "Error: $(CAPD_CONFIG) not found after installing"; exit 1; }

$(CAPD_CONFIG):
	@$(MAKE) setup

$(BINDIR):
	@mkdir -p "$(BINDIR)"

$(OBJDIR):
	@mkdir -p "$(OBJDIR)"

# ----- Link -----
$(PROGS:%=$(BINDIR)/%): $(BINDIR)/%: $(OBJDIR)/%.o $(OTHERS_OBJ) | $(BINDIR) $(CAPD_CONFIG)
	$(CXX) $(LDFLAGS) -o $@ $< $(OTHERS_OBJ) $(CAPD_LIBS)

# ----- Compile with dependency generation -----
-include $(OBJ_FILES:%=%.d)

$(OBJ_FILES): $(OBJDIR)/%.o: %.cpp | $(OBJDIR) $(CAPD_CONFIG)
	$(CXX) $(CXXFLAGS) -MT $@ -MD -MP -MF $@.d -c -o $@ $<

clean:
	rm -rf "$(OBJDIR)" "$(BINDIR)"

capd-clean:
	rm -rf "$(CAPD_BUILD_DIR)"

capd-purge:
	rm -rf "$(CAPD_SRC_DIR)" "$(CAPD_PREFIX)"

distclean: clean capd-clean
	rm -rf "$(CAPD_PREFIX)"

rebuild: clean all