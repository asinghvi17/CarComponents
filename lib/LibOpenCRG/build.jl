#!/usr/bin/env julia
# ===================================================================
# Build script for LibOpenCRG
# -------------------------------------------------------------------
# Compiles the vendored ASAM OpenCRG "baselib" C sources (csrc/src/*.c,
# csrc/inc/*.h) into a shared library at build/libLibOpenCRG.<ext>.
#
# The baselib sources have no interdependencies beyond their shared
# private header (crgBaseLibPrivate.h) and the C standard library, so a
# single flat `cc` invocation compiling all translation units at once is
# simpler and more reliable than driving CMake for this purpose.
#
# Usage:
#   julia lib/LibOpenCRG/build.jl
# or, from within the package directory:
#   julia build.jl
# or via include(...) from another Julia session/script.
# ===================================================================

const PKG_DIR = @__DIR__
const CSRC_INC_DIR = joinpath(PKG_DIR, "csrc", "inc")
const CSRC_SRC_DIR = joinpath(PKG_DIR, "csrc", "src")
const BUILD_DIR = joinpath(PKG_DIR, "build")

# Shared-library extension conventions per platform.
const LIB_EXT = Sys.isapple() ? "dylib" : Sys.islinux() ? "so" : Sys.iswindows() ? "dll" : "so"
const LIB_NAME = "libLibOpenCRG.$(LIB_EXT)"
const LIB_PATH = joinpath(BUILD_DIR, LIB_NAME)

function find_cc()
    # Respect a user-specified CC environment variable, otherwise fall back
    # to the system `cc` (Apple clang on macOS, gcc/clang elsewhere).
    cc = get(ENV, "CC", "cc")
    return cc
end

function build()
    mkpath(BUILD_DIR)

    sources = sort(filter(f -> endswith(f, ".c"), readdir(CSRC_SRC_DIR; join = true)))
    isempty(sources) && error("No vendored C sources found in $(CSRC_SRC_DIR)")

    cc = find_cc()

    # -shared -fPIC works uniformly across the clang/gcc invocations we care
    # about (Apple clang emits a Mach-O dylib, gcc/clang on Linux emit an
    # ELF .so). -std=c11 matches the upstream project's own build settings
    # (see c-api/cmake/OpenCRGCompilerSettings.cmake in the upstream repo).
    cmd = `$cc -shared -fPIC -O2 -std=c11 -I$(CSRC_INC_DIR) $(sources) -o $(LIB_PATH)`

    println("Building LibOpenCRG shared library...")
    println("  compiler: ", cc)
    println("  sources:  ", length(sources), " files from ", CSRC_SRC_DIR)
    println("  output:   ", LIB_PATH)
    println()
    println("Running: ", cmd)

    run(cmd)

    if isfile(LIB_PATH)
        println()
        println("Build succeeded: ", LIB_PATH, " (", filesize(LIB_PATH), " bytes)")
    else
        error("Build appeared to succeed but $(LIB_PATH) was not created")
    end

    return LIB_PATH
end

# Only run automatically when executed as a script (`julia build.jl`), not
# when `include`d for its constants/functions.
if abspath(PROGRAM_FILE) == @__FILE__
    build()
end
