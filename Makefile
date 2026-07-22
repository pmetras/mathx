# Set config to debug or release
config ?= debug

# The name of the Pony package
PACKAGE := mathx

# The compiler to use. In that case, a self-compiled compiler in release mode.
#To see verbose compilation, try:
#COMPILE_WITH := ponyc --verbose=4
COMPILE_WITH := ponyc

# Where the executable are placed: build/debug or build/release
BUILD_DIR ?= build/$(config)

# Package main sources are in a directory with the same name as the package, i.e. 'limits'
SRC_DIR := $(PACKAGE)

# Additional libraries and packages
LIB_PATH := .:../assertx:../pony_testx

# Examples are in the example directory
EXAMPLES_DIR := examples

# Unit tests are kept apart in the tests directory
TESTS_DIR := tests
TESTS_EXE := tests

# Bignum tests are in the tests_bignum directory
TESTS_BIGNUM_DIR := tests_bignum
TESTS_BIGNUM_EXE := tests_bignum
BIGNUM_DIR := bignum

# Bignum/GMP tests require libgmp and libmpfr (compiled separately to avoid name clashes)
TESTS_BIGNUM_GMP_DIR := tests_bignum_gmp
TESTS_BIGNUM_GMP_EXE := tests_bignum_gmp

# The name of the tests binary: build/debug/limits or build/release/limits
tests_binary := $(BUILD_DIR)/$(PACKAGE)
tests_bignum_binary := $(BUILD_DIR)/$(TESTS_BIGNUM_EXE)
tests_bignum_gmp_binary := $(BUILD_DIR)/$(TESTS_BIGNUM_GMP_EXE)

# API doc is generated into build/limits-docs
DOCS_DIR := build/$(PACKAGE)-docs

# Check that config is one of debug or release
ifdef config
	ifeq (,$(filter $(config),debug release))
		$(error Unknown configuration "$(config)")
	endif
endif

# You can set compiler options here depending on config build.
ifeq ($(config),release)
	PONYC = $(COMPILE_WITH) --path=$(LIB_PATH)
else
	PONYC = $(COMPILE_WITH) --debug --path=$(LIB_PATH)
endif

# Find all source files from the 'package' directory.
SOURCE_FILES := $(shell find $(SRC_DIR) -name '*.pony')

# Collect all examples from Pony files in example directory
EXAMPLE_SOURCE_FILES := $(shell find $(EXAMPLES_DIR) -name '*.pony')

# Collect all test files from tests directory
TEST_FILES := $(shell find $(TESTS_DIR) -name '*.pony')

# Collect all test files from tests_bignum directory
TEST_BIGNUM_FILES := $(shell find $(TESTS_BIGNUM_DIR) -name '*.pony')

# Collect all test files from tests_bignum_gmp directory (exclude symlinks; tracked via explicit dep below)
TEST_BIGNUM_GMP_FILES := $(shell find $(TESTS_BIGNUM_GMP_DIR) -name '*.pony' -not -type l)

# Collect all bignum source files
BIGNUM_SOURCE_FILES := $(shell find $(BIGNUM_DIR) -name '*.pony')


# Full tests is unit tests + examples
test: unit-tests bignum-tests build-examples ## Run unit tests, bignum tests, and examples

# To run the unit tests, run the tests binary with arguments --exclude=integration --sequential
unit-tests: $(tests_binary) ## Run mathx unit tests
	$(BUILD_DIR)/$(TESTS_EXE) --exclude=integration --sequential

# Run bignum unit tests
bignum-tests: $(tests_bignum_binary) ## Run bignum unit tests
	$(BUILD_DIR)/$(TESTS_BIGNUM_EXE) --exclude=integration --sequential

# Run bignum/gmp unit tests (requires libgmp and libmpfr)
bignum-gmp-tests: $(tests_bignum_gmp_binary) ## Run bignum/gmp unit tests (requires libgmp/libmpfr)
	$(BUILD_DIR)/$(TESTS_BIGNUM_GMP_EXE) --sequential

# How to create binary for tests: compile content of SRC_DIR into BUILD_DIR.
# I add a dependency on source files in case pre-processing is required
$(tests_binary): $(TEST_FILES) $(SOURCE_FILES) | $(BUILD_DIR)
	$(PONYC) -o ${BUILD_DIR} $(TESTS_DIR)

$(tests_bignum_binary): $(TEST_BIGNUM_FILES) $(BIGNUM_SOURCE_FILES) $(SOURCE_FILES) | $(BUILD_DIR)
	$(PONYC) -o ${BUILD_DIR} $(TESTS_BIGNUM_DIR)

$(tests_bignum_gmp_binary): $(TEST_BIGNUM_GMP_FILES) $(TESTS_BIGNUM_DIR)/_tests_bigreal_suites.pony $(BIGNUM_SOURCE_FILES) $(SOURCE_FILES) | $(BUILD_DIR)
	$(PONYC) -o ${BUILD_DIR} $(TESTS_BIGNUM_GMP_DIR)

# Build all examples: on all sub-directories in the example directory that contain Pony code,
# filtering on sub-directories that contains only FFI.
build-examples: $(SOURCE_FILES) $(EXAMPLES_SOURCE_FILES) | $(BUILD_DIR) ## Build all examples
	find $(EXAMPLES_DIR)/*/* -name '*.pony' -print | \
		xargs -n 1 dirname | \
		sort -u | \
		grep -v ffi- | \
		xargs -n 1 -I {} $(PONYC) -s --checktree -o $(BUILD_DIR) {}

clean: ## Clean build executables
	rm -rf $(BUILD_DIR)

realclean: ## Clean all build executables and documentation
	rm -rf build

# Build the documentation
$(DOCS_DIR): $(SOURCE_FILES)
	rm -rf $(DOCS_DIR)
	$(PONYC) --docs-public --pass=docs --output build $(SRC_DIR)

docs: $(DOCS_DIR) ## Build documentation

TAGS: ## Run ctags on project sources
	ctags --recurse=yes $(SRC_DIR)

# Build all
all: test ## Build all: source + tests + examples

# Create the build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

.PHONY: all clean realclean TAGS test unit-tests bignum-tests bignum-gmp-tests help

help: ## Print help on Make targets
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
