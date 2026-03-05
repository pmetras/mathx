# Set config to debug or release
config ?= debug

# The name of the Pony package
PACKAGE := mathx

# The compiler to use. In that case, a self-compiled compiler in release mode.
COMPILE_WITH := ponyc --verbose=4

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

# The name of the tests binary: build/debug/limits or build/release/limits
tests_binary := $(BUILD_DIR)/$(PACKAGE)

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


# Full tests is unit tests + examples
test: unit-tests build-examples ## Run unit tests and examples

# To run the unit tests, run the tests binary with arguments --exclude=integration --sequential
unit-tests: $(tests_binary) ## Run unit tests
	$(BUILD_DIR)/$(TESTS_EXE) --exclude=integration --sequential

# How to create binary for tests: compile content of SRC_DIR into BUILD_DIR.
# I add a dependency on source files in case pre-processing is required
$(tests_binary): $(TEST_FILES) $(SOURCE_FILES) | $(BUILD_DIR)
	$(PONYC) -o ${BUILD_DIR} $(TESTS_DIR)

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

.PHONY: all clean realclean TAGS test help

help: ## Print help on Make targets
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
