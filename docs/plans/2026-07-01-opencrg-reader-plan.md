# OpenCRG Reader Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement a pure-Julia, dependency-free parser for ASAM OpenCRG `.crg` road-surface files in `lib/OpenCRG`, producing world-frame `(u, v, X, Y, Z)` grids via a single batched forward transform.

**Architecture:** A layered parser (byte-safe header/payload split → section grouping → typed header structs → ASCII/binary payload decode → channel assembly) produces a `CRGData` struct. A separate transform stage (reference-line integration, lateral "miter normal" offset, banking/slope composition, `$ROAD_CRG_MODS` application) turns that into world-frame grids via `road_surface_grid`. Every numerically tricky step is cross-validated against the `lib/LibOpenCRG` C-library oracle (already built, 65/65 tests passing) rather than trusted on inspection alone.

**Tech Stack:** Julia (stdlib only: `Base`, `Test`), the sibling `lib/LibOpenCRG` package as a test-only dependency, three vendored ASAM OpenCRG example `.crg` files (Apache-2.0) as fixtures.

**Context for whoever executes this plan:** see `docs/plans/2026-07-01-opencrg-reader-design.md` for the design rationale. All format details below (field widths, NaN rules, integration formulas, MODS semantics) were verified against the upstream spec text and the reference C source (not guessed) — citations are in the design doc's git history / this session's research. Three fixture files already sit at `lib/OpenCRG/test/data/`: `handmade_curved_minimalist.crg` and `handmade_curved_banked_sloped.crg` (small ASCII/LRFI, no `$ROAD_CRG_MODS`, no reference-line end-anchoring), and `belgian_block.crg` (1.35MB real binary/KRBI scan, 1001×342 grid, also no MODS/end-anchoring). **None of the three exercise the backward+blend integration case (Task 11) or `$ROAD_CRG_MODS` (Task 16) — those need a small hand-constructed synthetic fixture, described in their tasks.** `lib/LibOpenCRG`'s UUID is `415980e8-35ba-4db9-9363-7da5222f1119`; `Test` stdlib's UUID (for `Project.toml`) is `8dfed614-e22c-5e08-85e1-65c5234f0b40`.

**Known limitation of the per-task "Run to verify" commands, Tasks 2 through 10:** `lib/OpenCRG/src/OpenCRG.jl` (Task 1) unconditionally `include`s `header.jl`, `payload.jl`, `read.jl`, `transform.jl` in that fixed order. Julia's `include` is not lazy or conditional, so `using OpenCRG` cannot succeed until ALL FOUR exist — which only first becomes true partway through Task 11 (`transform.jl` is the last one created). This means every task's literal `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; ...'` verification command, as written below for Tasks 2-10, will fail with a `SystemError` about whichever of those four files doesn't exist yet — **this is expected, not a sign anything is broken.** To verify a task's own test file in isolation before Task 11 lands, wrap only the source files that exist so far in a throwaway module instead of using `using OpenCRG`, e.g. (adjust the `include` list to whatever of `header.jl`/`payload.jl`/`read.jl` exists at that point in the plan):
```
julia -e '
using Test
module OpenCRG
    include("lib/OpenCRG/src/header.jl")
    include("lib/OpenCRG/src/payload.jl")   # only if it exists yet
end
const DATA = joinpath("lib/OpenCRG/test", "data")
include("lib/OpenCRG/test/test_XXX.jl")
'
```
Do not "fix" this by adding `isfile` guards around the `include`s in `OpenCRG.jl` itself — that would hide a real missing-file error in the finished package behind a silent no-op, for no benefit once Task 11 lands and all four files permanently exist.

---

### Task 1: Package skeleton and test harness

**Files:**
- Modify: `lib/OpenCRG/Project.toml`
- Modify: `lib/OpenCRG/src/OpenCRG.jl`
- Create: `lib/OpenCRG/test/runtests.jl`
- Verify present: `lib/OpenCRG/test/data/handmade_curved_minimalist.crg`, `lib/OpenCRG/test/data/handmade_curved_banked_sloped.crg`, `lib/OpenCRG/test/data/belgian_block.crg` (already vendored)

**Step 1: Update `lib/OpenCRG/Project.toml`**

Replace its contents with:

```toml
name = "OpenCRG"
uuid = "40a2bf50-15af-408a-8ff9-1cb7955b13ad"
version = "0.1.0"
authors = ["Anshul Singhvi <anshulsinghvi@gmail.com>"]

[extras]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
LibOpenCRG = "415980e8-35ba-4db9-9363-7da5222f1119"

[targets]
test = ["Test", "LibOpenCRG"]

[sources]
LibOpenCRG = {path = "../LibOpenCRG"}
```

**Step 2: Replace the stub module**

Replace `lib/OpenCRG/src/OpenCRG.jl` with:

```julia
module OpenCRG

include("header.jl")
include("payload.jl")
include("read.jl")
include("transform.jl")

export read_crg, road_surface_grid, CRGData

end # module OpenCRG
```

(`header.jl`, `payload.jl`, `read.jl`, `transform.jl` don't exist yet — later tasks create them. Julia will error on `include` until Task 2 creates `header.jl`; that's expected and resolves as soon as the next task's file exists.)

**Step 3: Create the test harness**

```julia
# lib/OpenCRG/test/runtests.jl
using Test
using OpenCRG

const DATA = joinpath(@__DIR__, "data")

@testset "OpenCRG.jl" begin
    include("test_header.jl")
    include("test_payload.jl")
    include("test_read.jl")
    include("test_transform.jl")
    include("test_crossvalidate.jl")
end
```

(These `test_*.jl` files don't exist yet either — each later task creates the one it needs, and `runtests.jl` will only fully pass once all tasks are done. That's expected mid-plan.)

**Step 4: Confirm the fixture files are present**

Run: `ls lib/OpenCRG/test/data/`
Expected: `belgian_block.crg`, `handmade_curved_banked_sloped.crg`, `handmade_curved_minimalist.crg`

**Step 5: Commit**

```bash
git add lib/OpenCRG/Project.toml lib/OpenCRG/src/OpenCRG.jl lib/OpenCRG/test/runtests.jl lib/OpenCRG/test/data/
git commit -m "OpenCRG: package skeleton, test harness, fixtures"
```

---

### Task 2: Header/payload byte split and section grouping

This is the lowest-level tokenizer. It must work on raw bytes, not a Julia `String`, because a binary (KRBI/KDBI) payload immediately follows the header and is not valid UTF-8 in general — decoding the *whole* file as a `String` first would be unsafe.

**Files:**
- Create: `lib/OpenCRG/src/header.jl`
- Create: `lib/OpenCRG/test/test_header.jl`

**Step 1: Write the failing test**

```julia
# lib/OpenCRG/test/test_header.jl
@testset "header tokenizing" begin
    @testset "find_header_end / split_lines on a real file" begin
        bytes = read(joinpath(DATA, "handmade_curved_minimalist.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        @test header_end <= length(bytes)
        header_lines = OpenCRG.split_lines(bytes[1:header_end-1])
        @test any(l -> startswith(l, "\$CT"), header_lines)
        # The line at header_end-1's start must be the "$$...." terminator
        term_line = OpenCRG.split_lines(bytes[1:header_end-1])[end] # not the terminator itself; sanity only
        @test !isempty(header_lines)
    end

    @testset "group_sections strips comments and groups by keyword" begin
        lines = [
            "\$ROAD_CRG",
            "* this whole line is a comment",
            "REFERENCE_LINE_INCREMENT = 1.0   ! inline comment here",
            "\$KD_DEFINITION",
            "#:LRFI",
        ]
        sections = OpenCRG.group_sections(lines)
        @test sections["ROAD_CRG"] == ["REFERENCE_LINE_INCREMENT = 1.0"]
        @test sections["KD_DEFINITION"] == ["#:LRFI"]
    end

    @testset "binary file: header/payload boundary is byte-exact" begin
        bytes = read(joinpath(DATA, "belgian_block.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        payload = bytes[header_end:end]
        # 1001 rows x 342 channels x 4 bytes, rounded up to a multiple of 80
        expected_padded = cld(1001 * 342 * 4, 80) * 80
        @test length(payload) == expected_padded
    end
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test; include("lib/OpenCRG/test/runtests.jl")'`
Expected: FAIL — `OpenCRG.find_header_end` (etc.) not defined, since `header.jl` doesn't exist and isn't `include`d successfully yet.

**Step 3: Write the implementation**

```julia
# lib/OpenCRG/src/header.jl

"""
    find_header_end(bytes) -> Int

Return the 1-based byte index of the first character after the header's
terminating `\$\$...` line (2+ leading `\$` characters) — i.e. where the
road-data payload begins. Returns `length(bytes) + 1` if no such line is
found (header-only file).
"""
function find_header_end(bytes::Vector{UInt8})
    pos, n = 1, length(bytes)
    while pos <= n
        nl = findnext(==(UInt8('\n')), bytes, pos)
        line_end = nl === nothing ? n : nl - 1
        next_pos = nl === nothing ? n + 1 : nl + 1
        if line_end - pos + 1 >= 2 && bytes[pos] == UInt8('$') && bytes[pos+1] == UInt8('$')
            return next_pos
        end
        pos = next_pos
    end
    return n + 1
end

"""
    split_lines(bytes) -> Vector{String}

Split raw bytes into LF-terminated lines (tolerating a preceding CR).
"""
function split_lines(bytes::AbstractVector{UInt8})
    lines = String[]
    pos, n = 1, length(bytes)
    while pos <= n
        nl = findnext(==(UInt8('\n')), bytes, pos)
        line_end = nl === nothing ? n : nl - 1
        real_end = (line_end >= pos && bytes[line_end] == UInt8('\r')) ? line_end - 1 : line_end
        push!(lines, String(bytes[pos:real_end]))
        pos = nl === nothing ? n + 1 : nl + 1
    end
    return lines
end

"""
    strip_inline_comment(line) -> String

Remove a `!` inline comment (everything from the first `!` onward).
"""
function strip_inline_comment(line::AbstractString)
    idx = findfirst('!', line)
    idx === nothing && return line
    return line[1:prevind(line, idx)]
end

"""Is `line`'s first non-space character `*` (a whole-line comment)?"""
function is_block_comment(line::AbstractString)
    s = lstrip(line)
    return !isempty(s) && s[1] == '*'
end

"""
    section_marker(line) -> Union{String,Nothing}

`nothing` if `line` isn't a section marker; otherwise the upper-cased
keyword (empty string for a bare `\$` section-end marker).
"""
function section_marker(line::AbstractString)
    (isempty(line) || line[1] != '$') && return nothing
    return uppercase(strip(line[2:end]))
end

"""
    group_sections(header_lines) -> Dict{String,Vector{String}}

Group header lines into sections keyed by keyword (e.g. `"ROAD_CRG"`),
stripping block (`*`) and inline (`!`) comments and blank lines. A bare `\$`
line closes the current section without opening a new one.
"""
function group_sections(header_lines::Vector{String})
    sections = Dict{String,Vector{String}}()
    current = ""
    for raw in header_lines
        kw = section_marker(raw)
        if kw !== nothing
            current = kw
            isempty(kw) && continue
            haskey(sections, current) || (sections[current] = String[])
            continue
        end
        is_block_comment(raw) && continue
        content = strip(strip_inline_comment(raw))
        (isempty(content) || isempty(current)) && continue
        push!(sections[current], content)
    end
    return sections
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test; include("lib/OpenCRG/test/runtests.jl")'`
Expected: The `"header tokenizing"` testset passes (later testsets will still error since `test_payload.jl` etc. don't exist — that's expected until their tasks land; for now run just this file directly to confirm in isolation: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_header.jl")'`)

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/header.jl lib/OpenCRG/test/test_header.jl
git commit -m "OpenCRG: byte-safe header/payload split and section grouping"
```

---

### Task 3: `$ROAD_CRG` key-value parsing → `ReferenceLineParams`

**Files:**
- Modify: `lib/OpenCRG/src/header.jl`
- Modify: `lib/OpenCRG/test/test_header.jl`

**Step 1: Write the failing test**

Append to `lib/OpenCRG/test/test_header.jl` (inside the existing `@testset "header tokenizing"` block, add a new nested testset):

```julia
    @testset "parse_road_crg on a real file" begin
        bytes = read(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        sections = OpenCRG.group_sections(OpenCRG.split_lines(bytes[1:header_end-1]))
        r = OpenCRG.parse_road_crg(sections["ROAD_CRG"])
        @test r.start_u == 0.0
        @test r.end_u == 22.0
        @test r.increment == 1.0
        @test r.start_y == 0.0
        @test r.start_phi == 0.0
        @test r.end_x === nothing   # this file has no explicit end position
        @test r.v_right == -1.5
        @test r.v_left == 1.5
        @test r.start_slope == 0.0   # default when REFERENCE_LINE_START_S is absent
        @test r.start_banking == 0.0
    end

    @testset "parse_keyvalues / parse_keyvalue_strings" begin
        @test OpenCRG.parse_keyvalues(["FOO = 1.5", "BAR=2"]) == Dict("FOO"=>1.5, "BAR"=>2.0)
        @test OpenCRG.parse_keyvalue_strings(["PROJ_NM = UTM", "PROJ_ZONE = 32"]) ==
            Dict("PROJ_NM"=>"UTM", "PROJ_ZONE"=>"32")
    end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_header.jl")'`
Expected: FAIL — `parse_road_crg`/`parse_keyvalues`/`parse_keyvalue_strings` not defined.

**Step 3: Write the implementation**

Append to `lib/OpenCRG/src/header.jl`:

```julia
"""
    parse_keyvalues(lines) -> Dict{String,Float64}

Parse `KEY = value` lines into a dict keyed by upper-cased key. Lines whose
value doesn't parse as a float are skipped (see `parse_keyvalue_strings` for
sections that mix numeric and string fields).
"""
function parse_keyvalues(lines::Vector{String})
    dict = Dict{String,Float64}()
    for line in lines
        eq = findfirst('=', line)
        eq === nothing && continue
        key = uppercase(strip(line[1:prevind(line, eq)]))
        v = tryparse(Float64, strip(line[nextind(line, eq):end]))
        v === nothing || (dict[key] = v)
    end
    return dict
end

"""
    parse_keyvalue_strings(lines) -> Dict{String,String}

Like `parse_keyvalues` but keeps every value as a raw string.
"""
function parse_keyvalue_strings(lines::Vector{String})
    dict = Dict{String,String}()
    for line in lines
        eq = findfirst('=', line)
        eq === nothing && continue
        key = uppercase(strip(line[1:prevind(line, eq)]))
        dict[key] = strip(line[nextind(line, eq):end])
    end
    return dict
end

"""
Reference-line and long-section parameters from `\$ROAD_CRG`. Optional
fields (anchoring end position/heading/elevation, explicit constant
slope/banking) are `Union{Float64,Nothing}`, defaulting to `nothing`
except where the spec gives a numeric default (slope/banking default 0.0).
"""
struct ReferenceLineParams
    start_u::Float64
    end_u::Float64
    increment::Float64
    start_x::Float64
    start_y::Float64
    start_phi::Float64
    end_x::Union{Float64,Nothing}
    end_y::Union{Float64,Nothing}
    end_phi::Union{Float64,Nothing}
    start_z::Float64
    end_z::Union{Float64,Nothing}
    v_right::Union{Float64,Nothing}
    v_left::Union{Float64,Nothing}
    v_increment::Union{Float64,Nothing}
    start_slope::Float64
    end_slope::Union{Float64,Nothing}
    start_banking::Float64
    end_banking::Union{Float64,Nothing}
end

function parse_road_crg(lines::Vector{String})
    d = parse_keyvalues(lines)
    g(k, default=0.0) = get(d, k, default)
    go(k) = get(d, k, nothing)
    return ReferenceLineParams(
        g("REFERENCE_LINE_START_U"), g("REFERENCE_LINE_END_U"), g("REFERENCE_LINE_INCREMENT"),
        g("REFERENCE_LINE_START_X"), g("REFERENCE_LINE_START_Y"), g("REFERENCE_LINE_START_PHI"),
        go("REFERENCE_LINE_END_X"), go("REFERENCE_LINE_END_Y"), go("REFERENCE_LINE_END_PHI"),
        g("REFERENCE_LINE_START_Z"), go("REFERENCE_LINE_END_Z"),
        go("LONG_SECTION_V_RIGHT"), go("LONG_SECTION_V_LEFT"), go("LONG_SECTION_V_INCREMENT"),
        g("REFERENCE_LINE_START_S"), go("REFERENCE_LINE_END_S"),
        g("REFERENCE_LINE_START_B"), go("REFERENCE_LINE_END_B"),
    )
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_header.jl")'`
Expected: PASS (all testsets in the file green)

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/header.jl lib/OpenCRG/test/test_header.jl
git commit -m "OpenCRG: parse \$ROAD_CRG into ReferenceLineParams"
```

---

### Task 4: `$KD_DEFINITION` parsing

**Files:**
- Modify: `lib/OpenCRG/src/header.jl`
- Create: `lib/OpenCRG/test/test_kddef.jl`
- Modify: `lib/OpenCRG/test/runtests.jl` (add `include("test_kddef.jl")`)

**Step 1: Write the failing test**

```julia
# lib/OpenCRG/test/test_kddef.jl
@testset "\$KD_DEFINITION parsing" begin
    @testset "position-form (at v = ...), U: lines ignored" begin
        lines = [
            "#:LRFI",
            "U:reference line u,m,0,1.0",
            "D:reference line phi,rad",
            "D:reference line banking,m/m",
            "D:reference line slope,m/m",
            "D:long section at v = -1.500,m",
            "D:long section at v =  0.000,m",
            "D:long section at v =  1.500,m",
        ]
        format_code, channels = OpenCRG.parse_kd_definition(lines)
        @test format_code == :LRFI
        @test length(channels) == 6   # the U: line contributes nothing
        @test [c.kind for c in channels] == [:phi, :banking, :slope, :long_section, :long_section, :long_section]
        @test [c.v for c in channels if c.kind == :long_section] == [-1.5, 0.0, 1.5]
    end

    @testset "index-form (bare N), default format code" begin
        lines = ["D:long section 1,m", "D:long section 2,m", "D:long section 3,m"]
        format_code, channels = OpenCRG.parse_kd_definition(lines)
        @test format_code == :KRBI   # default when no #: line is present
        @test [c.index for c in channels] == [1, 2, 3]
    end

    @testset "mixed position/index form is an error" begin
        lines = ["D:long section at v = 0.0,m", "D:long section 2,m"]
        @test_throws Exception OpenCRG.parse_kd_definition(lines)
    end

    @testset "v_axis: index-form uniform spacing" begin
        r = OpenCRG.parse_road_crg(["LONG_SECTION_V_RIGHT = -1.0", "LONG_SECTION_V_LEFT = 1.0", "LONG_SECTION_V_INCREMENT = 1.0"])
        channels = [OpenCRG.ChannelDef(:long_section, nothing, i) for i in 1:3]
        @test OpenCRG.v_axis(channels, r) == [-1.0, 0.0, 1.0]
    end

    @testset "v_axis: position-form reads channel v directly, including non-uniform spacing" begin
        r = OpenCRG.parse_road_crg(String[])
        channels = [OpenCRG.ChannelDef(:long_section, v, nothing) for v in [-1.0, 0.0, 0.3]]
        @test OpenCRG.v_axis(channels, r) == [-1.0, 0.0, 0.3]
    end

    @testset "real file: handmade_curved_banked_sloped.crg" begin
        bytes = read(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        sections = OpenCRG.group_sections(OpenCRG.split_lines(bytes[1:header_end-1]))
        format_code, channels = OpenCRG.parse_kd_definition(sections["KD_DEFINITION"])
        @test format_code == :LRFI
        @test length(channels) == 10   # phi, banking, slope, 7 long sections
    end
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_kddef.jl")'`
Expected: FAIL — `ChannelDef`/`parse_kd_definition`/`v_axis` not defined.

**Step 3: Write the implementation**

Append to `lib/OpenCRG/src/header.jl`:

```julia
"""
One `D:` channel from `\$KD_DEFINITION`. `v`/`index` are only set for
`kind == :long_section`, in the position-form / index-form respectively
(mutually exclusive — see `parse_kd_definition`). `U:` (virtual) channels
are discarded entirely during parsing and never produce a `ChannelDef`,
matching the reference implementation's `decodeIndependent`, which is a
documented no-op.
"""
struct ChannelDef
    kind::Symbol   # :phi, :banking, :slope, :long_section
    v::Union{Float64,Nothing}
    index::Union{Int,Nothing}
end

const CRG_FORMAT_CODES = Dict("LRFI"=>:LRFI, "LDFI"=>:LDFI, "KRBI"=>:KRBI, "KDBI"=>:KDBI)

"""
    parse_kd_definition(lines) -> (format_code::Symbol, channels::Vector{ChannelDef})

Parse the `#:` format-code line (default `:KRBI` if absent, per spec) and
`D:` channel lines, in declaration order. `D:long section at v = X,unit`
(position-form) and `D:long section N,unit` (index-form) cannot be mixed in
one file — the reference loader treats that as fatal, and so do we.
"""
function parse_kd_definition(lines::Vector{String})
    format_code = :KRBI
    channels = ChannelDef[]
    v_mode = nothing
    for line in lines
        if startswith(line, "#:")
            code = uppercase(strip(line[3:end]))
            haskey(CRG_FORMAT_CODES, code) || error("unknown CRG data format code: $code")
            format_code = CRG_FORMAT_CODES[code]
        elseif startswith(line, "U:")
            continue
        elseif startswith(line, "D:")
            label = line[3:end]
            comma = findfirst(',', label)
            name = comma === nothing ? label : label[1:prevind(label, comma)]
            name_lower = lowercase(strip(name))
            if occursin("reference line phi", name_lower)
                push!(channels, ChannelDef(:phi, nothing, nothing))
            elseif occursin("reference line bank", name_lower)
                push!(channels, ChannelDef(:banking, nothing, nothing))
            elseif occursin("reference line slope", name_lower)
                push!(channels, ChannelDef(:slope, nothing, nothing))
            elseif occursin("long section", name_lower)
                atidx = findfirst("at v", name_lower)
                if atidx !== nothing
                    v_mode == :index && error("\$KD_DEFINITION mixes explicit and implicit long-section v definitions")
                    v_mode = :position
                    eqidx = findfirst('=', name)
                    v = parse(Float64, strip(name[nextind(name, eqidx):end]))
                    push!(channels, ChannelDef(:long_section, v, nothing))
                else
                    v_mode == :position && error("\$KD_DEFINITION mixes explicit and implicit long-section v definitions")
                    v_mode = :index
                    idx = parse(Int, strip(replace(name_lower, "long section" => "")))
                    push!(channels, ChannelDef(:long_section, nothing, idx))
                end
            else
                error("unrecognized \$KD_DEFINITION channel label: $name")
            end
        end
    end
    return format_code, channels
end

"""
    v_axis(channels, refline) -> Vector{Float64}

The v-coordinate for each `:long_section` channel, in declaration order.
Index-form derives a uniform axis from `LONG_SECTION_V_RIGHT/_LEFT/_INCREMENT`;
position-form reads each channel's own parsed `v` directly (including truly
non-uniform spacing — unlike the reference C implementation, which snaps
near-uniform spacings and leaves a genuinely non-uniform axis in a code path
its own authors flag as possibly incomplete, we don't need to worry about
this at all: OpenCRG.jl does no interpolation itself, so non-uniform v is
just data, not a numerical hazard).
"""
function v_axis(channels::Vector{ChannelDef}, refline::ReferenceLineParams)
    long_sections = filter(c -> c.kind == :long_section, channels)
    isempty(long_sections) && return Float64[]
    if long_sections[1].index !== nothing
        refline.v_right === nothing && error("index-form long sections require LONG_SECTION_V_RIGHT")
        n = length(long_sections)
        inc = refline.v_increment
        if inc === nothing
            n == 1 && error("LONG_SECTION_V_INCREMENT required when there is only one long section (spacing can't be derived from a single point)")
            inc = (refline.v_left - refline.v_right) / (n - 1)
        end
        return [refline.v_right + (c.index - 1) * inc for c in long_sections]
    else
        return [c.v for c in long_sections]
    end
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_kddef.jl")'`
Expected: PASS

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/header.jl lib/OpenCRG/test/test_kddef.jl lib/OpenCRG/test/runtests.jl
git commit -m "OpenCRG: parse \$KD_DEFINITION and derive the v-axis"
```

---

### Task 5: `$ROAD_CRG_OPTS` / `$ROAD_CRG_MPRO` (parsed, not applied)

**Files:**
- Modify: `lib/OpenCRG/test/test_header.jl` (add a testset — no new src code needed, these reuse Task 3's generic parsers)

**Step 1: Write the failing test**

Append a testset to `lib/OpenCRG/test/test_header.jl`:

```julia
    @testset "\$ROAD_CRG_OPTS / \$ROAD_CRG_MPRO parse as generic dicts (not applied)" begin
        opts = OpenCRG.parse_keyvalues(["BORDER_MODE_U = 2", "BORDER_MODE_V = 0"])
        @test opts["BORDER_MODE_U"] == 2.0

        mpro = OpenCRG.parse_keyvalue_strings(["PROJ_NM = UTM", "GELL_A = 6378137.0"])
        @test mpro["PROJ_NM"] == "UTM"
        @test mpro["GELL_A"] == "6378137.0"   # kept as a string; no geodesy math is implemented
    end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_header.jl")'`
Expected: this specific testset should actually already PASS, since `parse_keyvalues`/`parse_keyvalue_strings` already exist from Task 3 — there is no new implementation here. This task exists only to pin down, with an explicit test, that OPTS/MPRO deliberately get no dedicated struct or applied behavior (a scope boundary worth a regression test, not just a comment).

**Step 3: (no implementation needed — reusing Task 3's generic parsers)**

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_header.jl")'`
Expected: PASS

**Step 5: Commit**

```bash
git add lib/OpenCRG/test/test_header.jl
git commit -m "OpenCRG: pin down that \$ROAD_CRG_OPTS/\$ROAD_CRG_MPRO are parsed but inert"
```

---

### Task 6: `$ROAD_CRG_MODS` structured parsing

**Files:**
- Modify: `lib/OpenCRG/src/header.jl`
- Create: `lib/OpenCRG/test/test_mods.jl`
- Modify: `lib/OpenCRG/test/runtests.jl` (add `include("test_mods.jl")`)

**Step 1: Write the failing test**

```julia
# lib/OpenCRG/test/test_mods.jl
@testset "\$ROAD_CRG_MODS parsing" begin
    lines = [
        "SCALE_Z_GRID = 2.0",
        "REFLINE_OFFSET_PHI = 1.57",
        "REFLINE_OFFSET_X = 100.0",
        "REFPOINT_PHI = 0.0",   # presence alone should flip has_refpoint on, even though it's 0.0
    ]
    m = OpenCRG.parse_road_crg_mods(lines)
    @test m.scale_z_grid == 2.0
    @test m.scale_slope === nothing
    @test m.refline_offset_phi == 1.57
    @test m.refline_offset_x == 100.0
    @test m.refpoint_phi == 0.0

    empty_mods = OpenCRG.parse_road_crg_mods(String[])
    @test empty_mods.scale_z_grid === nothing
    @test empty_mods.refpoint_phi === nothing
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_mods.jl")'`
Expected: FAIL — `RoadCrgMods`/`parse_road_crg_mods` not defined.

**Step 3: Write the implementation**

Append to `lib/OpenCRG/src/header.jl`:

```julia
"""
Every field from `\$ROAD_CRG_MODS`. All optional (`nothing` if absent from
the file) — `nothing` here means "not present," which matters: e.g.
`has_refpoint` in `apply_mods` (Task 16) checks presence of *any*
`refpoint_*` field, not whether it's zero.
"""
Base.@kwdef struct RoadCrgMods
    scale_z_grid::Union{Float64,Nothing} = nothing
    scale_slope::Union{Float64,Nothing} = nothing
    scale_banking::Union{Float64,Nothing} = nothing
    scale_length::Union{Float64,Nothing} = nothing
    scale_width::Union{Float64,Nothing} = nothing
    scale_curvature::Union{Float64,Nothing} = nothing
    grid_nan_mode::Union{Int,Nothing} = nothing
    grid_nan_offset::Union{Float64,Nothing} = nothing
    refpoint_u::Union{Float64,Nothing} = nothing
    refpoint_u_fraction::Union{Float64,Nothing} = nothing
    refpoint_u_offset::Union{Float64,Nothing} = nothing
    refpoint_v::Union{Float64,Nothing} = nothing
    refpoint_v_fraction::Union{Float64,Nothing} = nothing
    refpoint_v_offset::Union{Float64,Nothing} = nothing
    refpoint_x::Union{Float64,Nothing} = nothing
    refpoint_y::Union{Float64,Nothing} = nothing
    refpoint_z::Union{Float64,Nothing} = nothing
    refpoint_phi::Union{Float64,Nothing} = nothing
    refline_rotcenter_x::Union{Float64,Nothing} = nothing
    refline_rotcenter_y::Union{Float64,Nothing} = nothing
    refline_offset_x::Union{Float64,Nothing} = nothing
    refline_offset_y::Union{Float64,Nothing} = nothing
    refline_offset_z::Union{Float64,Nothing} = nothing
    refline_offset_phi::Union{Float64,Nothing} = nothing
end

"""
    parse_road_crg_mods(lines) -> RoadCrgMods

Parse `\$ROAD_CRG_MODS` section lines into a `RoadCrgMods`.
"""
function parse_road_crg_mods(lines::Vector{String})
    d = parse_keyvalues(lines)
    go(k) = get(d, k, nothing)
    goi(k) = (v = get(d, k, nothing); v === nothing ? nothing : Int(v))
    return RoadCrgMods(;
        scale_z_grid=go("SCALE_Z_GRID"), scale_slope=go("SCALE_SLOPE"), scale_banking=go("SCALE_BANKING"),
        scale_length=go("SCALE_LENGTH"), scale_width=go("SCALE_WIDTH"), scale_curvature=go("SCALE_CURVATURE"),
        grid_nan_mode=goi("GRID_NAN_MODE"), grid_nan_offset=go("GRID_NAN_OFFSET"),
        refpoint_u=go("REFPOINT_U"), refpoint_u_fraction=go("REFPOINT_U_FRACTION"), refpoint_u_offset=go("REFPOINT_U_OFFSET"),
        refpoint_v=go("REFPOINT_V"), refpoint_v_fraction=go("REFPOINT_V_FRACTION"), refpoint_v_offset=go("REFPOINT_V_OFFSET"),
        refpoint_x=go("REFPOINT_X"), refpoint_y=go("REFPOINT_Y"), refpoint_z=go("REFPOINT_Z"), refpoint_phi=go("REFPOINT_PHI"),
        refline_rotcenter_x=go("REFLINE_ROTCENTER_X"), refline_rotcenter_y=go("REFLINE_ROTCENTER_Y"),
        refline_offset_x=go("REFLINE_OFFSET_X"), refline_offset_y=go("REFLINE_OFFSET_Y"),
        refline_offset_z=go("REFLINE_OFFSET_Z"), refline_offset_phi=go("REFLINE_OFFSET_PHI"),
    )
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_mods.jl")'`
Expected: PASS

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/header.jl lib/OpenCRG/test/test_mods.jl lib/OpenCRG/test/runtests.jl
git commit -m "OpenCRG: parse \$ROAD_CRG_MODS into a structured RoadCrgMods"
```

---

### Task 7: ASCII payload decoding (LRFI/LDFI)

**Files:**
- Create: `lib/OpenCRG/src/payload.jl`
- Create: `lib/OpenCRG/test/test_payload.jl`

**Step 1: Write the failing test**

```julia
# lib/OpenCRG/test/test_payload.jl
@testset "ASCII payload decoding" begin
    @testset "decode_ascii_field" begin
        @test OpenCRG.decode_ascii_field("**unused**") == 0.0   # exact literal only
        @test isnan(OpenCRG.decode_ascii_field("*missing*"))
        @test OpenCRG.decode_ascii_field(" 0.0111111") == 0.0111111
        @test OpenCRG.decode_ascii_field("-1.500000") == -1.5
        @test OpenCRG.decode_ascii_field("1.0D+02") == 100.0   # Fortran D-exponent normalization
        @test isnan(OpenCRG.decode_ascii_field("**unused**" * " "^10))   # 20-char LDFI-width field: never matches the 10-char literal
    end

    @testset "decode_ascii_payload: exact multiple of per_record, no line wrap (8 channels at LRFI)" begin
        # Naive single-space-joined fields don't add up to 8x10 = 80 chars; each
        # field must be individually padded to its full 10-char width.
        lines = [" -1.000000  0.000000  1.000000  2.000000  3.000000  4.000000  5.000000  6.000000"]
        m = OpenCRG.decode_ascii_payload(lines, :LRFI, 8)
        @test size(m) == (1, 8)
        @test m[1, :] == [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    end

    @testset "decode_ascii_payload: non-whole-row remainder is an error, not a silent truncation" begin
        lines = [
            "**unused** 0.0000000**unused** 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000",
            " 0.0000000 0.0000000",
            " 0.0000000 0.0000000",   # a stray 3rd line: 3 lines don't divide evenly by 2 lines/row
        ]
        @test_throws Exception OpenCRG.decode_ascii_payload(lines, :LRFI, 10)
    end

    @testset "decode_ascii_payload: row-wrapping, 10 channels at LRFI (8/record)" begin
        # Row 0 from handmade_curved_banked_sloped.crg, 10 channels, wraps onto 2 lines.
        lines = [
            "**unused** 0.0000000**unused** 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000",
            " 0.0000000 0.0000000",
        ]
        m = OpenCRG.decode_ascii_payload(lines, :LRFI, 10)
        @test size(m) == (1, 10)
        @test m[1,1] == 0.0          # **unused** -> 0.0
        @test m[1,2] == 0.0
        @test m[1,3] == 0.0          # slope, row 0: also the exact **unused** literal in this fixture, not *missing* (see note)
    end

    @testset "real file end-to-end shape" begin
        bytes = read(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        payload_lines = OpenCRG.split_lines(bytes[header_end:end])
        m = OpenCRG.decode_ascii_payload(payload_lines, :LRFI, 10)
        @test size(m) == (23, 10)    # nu is DERIVED from payload size: 46 lines / 2 lines-per-row = 23
        @test all(isfinite, m[2:end, 1])   # phi is defined for every row except row 1 (index 1 in Julia)
    end
end
```

Note on the second testset: this file's channel order is phi, banking, slope, then 7 long-sections — channel 3 is slope, not banking. Confirmed by reading the real bytes: row 0/channel 3 is the exact literal `**unused**` (not `*missing*`), hence `== 0.0` above rather than `isnan`.

Add a synthetic LDFI test proving the `**unused**` quirk's width-sensitivity, to the same `@testset "decode_ascii_field"` block:
```julia
        @test isnan(OpenCRG.decode_ascii_field("**unused**" * " "^10))   # 20-char LDFI-width field: never matches the 10-char literal
```
This is the one case in the whole task that's easy to get backwards: it's tempting to `strip` a field before checking it against `"**unused**"`, but that destroys exactly the width information that makes this quirk correctly NOT apply to LDFI's wider fields in the reference implementation (a 20-byte buffer can never `strcmp`-equal a 10-character literal, no matter its padding) — see `decode_ascii_field`'s docstring below for why the strip has to happen *inside* the function, after the exact-literal check, not by the caller beforehand.

**Why `decode_ascii_payload` doesn't take `nu` as a parameter:** an earlier draft of this plan had `read_crg` (Task 10) compute `nu` from `\$ROAD_CRG`'s `REFERENCE_LINE_START_U`/`END_U`/`INCREMENT` and pass it in. That's wrong in general — `handmade_curved_minimalist.crg` (already vendored, used in Task 10's own tests) has NO `REFERENCE_LINE_START_U`/`END_U` at all (its `\$ROAD_CRG` section only sets `REFERENCE_LINE_INCREMENT`), yet its payload has 23 real data rows. `start_u`/`end_u` are optional, sometimes-redundant descriptive metadata — the payload's actual row count is the only ground truth for `nu`. So `nu` is derived here, from `length(payload_lines)`, not threaded in from the header.

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_payload.jl")'`
Expected: FAIL — `decode_ascii_field`/`decode_ascii_payload` not defined.

**Step 3: Write the implementation**

```julia
# lib/OpenCRG/src/payload.jl

const ASCII_FIELD_WIDTH = Dict(:LRFI => 10, :LDFI => 20)
const FIELDS_PER_RECORD = Dict(:LRFI => 8, :LDFI => 4, :KRBI => 20, :KDBI => 10)
const BINARY_ELEM_TYPE = Dict(:KRBI => Float32, :KDBI => Float64)

"""
    decode_ascii_field(raw) -> Float64

Decode one fixed-width ASCII field, exactly as sliced from the payload line
— NOT pre-stripped by the caller, since the exact-width comparison below
depends on seeing any trailing padding. The literal 10-character content
`**unused**`, compared against the RAW (unstripped) field, decodes to `0.0`.
This quirk structurally cannot apply to LDFI's 20-character fields: a
20-byte raw field can never string-equal the 10-character literal (their
lengths differ), so it always falls through to the numeric path below —
matching the reference C decoder's `strcmp` check, which compares the full
fixed-width buffer, not a pre-trimmed one. (An earlier version of this
function had the caller `strip` before calling it, which silently broke
this width-sensitivity — a stripped 20-char `"**unused**"` + padding
collapses to the bare 10-character literal and would incorrectly match.)
Once the exact-unused check fails, the field is `strip`ped and parsed as a
float; any character outside `0123456789+-.eEdD` after stripping decodes to
NaN (word-agnostic: `*missing*`, `*none*`, etc. all become NaN).
Fortran-style `D`/`d` exponent markers are normalized to `e` before `parse`.
"""
function decode_ascii_field(raw::AbstractString)
    raw == "**unused**" && return 0.0
    s = strip(raw)
    if !all(c -> c in "0123456789+-.eEdD", s)
        return NaN
    end
    normalized = replace(replace(s, 'D' => 'e'), 'd' => 'e')
    v = tryparse(Float64, normalized)
    return v === nothing ? NaN : v
end

"""
    decode_ascii_payload(payload_lines, format_code, nchannels) -> Matrix{Float64}

Decode `nchannels`-wide ASCII payload rows (`:LRFI`: 10-char fields,
8/record; `:LDFI`: 20-char fields, 4/record). The row count `nu` is
DERIVED as `length(payload_lines) ÷ lines_per_row` — see the note in Task
7 above for why this can't come from `\$ROAD_CRG`'s `start_u`/`end_u`. It's
an error for the payload to have a number of lines that isn't a whole
multiple of `lines_per_row` — silently floor-dividing would quietly drop a
truncated/corrupted file's leftover partial row instead of surfacing it.
Each row always starts on a fresh line — per spec §6.4.8, plain-text rows
never pack a previous row's leftover fields onto the same line, unlike
binary (Task 8).
"""
function decode_ascii_payload(payload_lines::Vector{<:AbstractString}, format_code::Symbol, nchannels::Int)
    width = ASCII_FIELD_WIDTH[format_code]
    per_record = FIELDS_PER_RECORD[format_code]
    lines_per_row = cld(nchannels, per_record)
    rem(length(payload_lines), lines_per_row) == 0 ||
        error("ASCII payload has $(length(payload_lines)) lines, not a whole multiple of $lines_per_row lines/row for $nchannels channels at $format_code")
    nu = length(payload_lines) ÷ lines_per_row
    data = Matrix{Float64}(undef, nu, nchannels)
    line_idx = 0
    for row in 1:nu
        ch = 0
        for _ in 1:lines_per_row
            line_idx += 1
            line = payload_lines[line_idx]
            nfields_here = min(per_record, nchannels - ch)
            for f in 1:nfields_here
                lo = (f - 1) * width + 1
                hi = min(f * width, ncodeunits(line))
                field = lo <= ncodeunits(line) ? line[lo:hi] : ""
                ch += 1
                data[row, ch] = decode_ascii_field(field)
            end
        end
    end
    return data
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_payload.jl")'`
Expected: PASS (after correcting the `**unused**`-vs-`*missing*` assertion per the note in Step 1 if needed)

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/payload.jl lib/OpenCRG/test/test_payload.jl
git commit -m "OpenCRG: decode ASCII (LRFI/LDFI) payload rows"
```

---

### Task 8: Binary payload decoding (KRBI/KDBI)

**Files:**
- Modify: `lib/OpenCRG/src/payload.jl`
- Modify: `lib/OpenCRG/test/test_payload.jl`

**Step 1: Write the failing test**

Append to `lib/OpenCRG/test/test_payload.jl`:

```julia
@testset "binary payload decoding" begin
    @testset "synthetic 2-row x 3-channel KRBI blob catches row/column transposition" begin
        # Row 0 = [1.0, 2.0, 3.0], Row 1 = [4.0, 5.0, 6.0], packed tightly,
        # big-endian Float32, no padding between rows (unlike ASCII).
        vals = Float32[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        bytes = UInt8[]
        for v in vals
            append!(bytes, reverse(reinterpret(UInt8, [v])))  # host is little-endian -> reverse for big-endian
        end
        m = OpenCRG.decode_binary_payload(bytes, :KRBI, 3)
        @test size(m) == (2, 3)
        @test m[1, :] == [1.0, 2.0, 3.0]
        @test m[2, :] == [4.0, 5.0, 6.0]
    end

    @testset "synthetic 2-row x 3-channel KDBI blob (Float64, no upstream fixture covers this)" begin
        vals = Float64[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        bytes = UInt8[]
        for v in vals
            append!(bytes, reverse(reinterpret(UInt8, [v])))
        end
        m = OpenCRG.decode_binary_payload(bytes, :KDBI, 3)
        @test size(m) == (2, 3)
        @test m[1, :] == [1.0, 2.0, 3.0]
        @test m[2, :] == [4.0, 5.0, 6.0]
    end

    @testset "truncated row is an error, not a silent drop" begin
        # 25 channels so row_bytes = 100 > 80: a real remainder can never
        # exceed 80 in a valid file, but a small nchannels count (e.g. the
        # 3-channel case above, row_bytes=12) could never produce a remainder
        # >= 80 either way, so this needs a wider row to actually exercise
        # the check.
        row = Float32.(1:25)
        bytes = UInt8[]
        for v in row
            append!(bytes, reverse(reinterpret(UInt8, [v])))
        end
        partial_second_row = bytes[1:90]   # 90 leftover bytes: not plausible padding (>= 80)
        @test_throws Exception OpenCRG.decode_binary_payload(vcat(bytes, partial_second_row), :KRBI, 25)
    end

    @testset "real binary file: shape and no per-row 80-byte alignment" begin
        bytes = read(joinpath(DATA, "belgian_block.crg"))
        header_end = OpenCRG.find_header_end(bytes)
        payload = bytes[header_end:end]
        m = OpenCRG.decode_binary_payload(payload, :KRBI, 342)
        @test size(m) == (1001, 342)   # nu is DERIVED: 1369440 bytes ÷ 1368 bytes/row = 1001 (72 trailing padding bytes ignored)
        @test isnan(m[1, 1])                    # row 0 phi placeholder, per spec/research
        @test m[2, 1] ≈ 2.6527974605560303       # bit-exact match to REFERENCE_LINE_START_PHI in this file's header
        @test all(isnan, m[1, 2:22])             # first cross-section's missing left-border samples (verified byte-for-byte: channels 2-22 are NaN)
        @test isfinite(m[1, 23])                 # first real elevation sample in row 0 (verified: channel 23 == 2.1214599609375)
    end
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_payload.jl")'`
Expected: FAIL — `decode_binary_payload` not defined.

**Step 3: Write the implementation**

Append to `lib/OpenCRG/src/payload.jl`:

```julia
"""
    decode_binary_payload(payload, format_code, nchannels) -> Matrix{Float64}

Decode `nchannels`-wide big-endian binary payload rows (`Float32` for
`:KRBI`, `Float64` for `:KDBI`). The row count `nu` is DERIVED via floor
division, `length(payload) ÷ (nchannels * sizeof(elem))` — trailing padding
bytes (the payload is padded to a multiple of 80 bytes overall, per spec)
are simply leftover and ignored, for the same reason `decode_ascii_payload`
(Task 7) derives its row count from the payload rather than trusting
`\$ROAD_CRG`'s `start_u`/`end_u`. **Rows are packed tightly with no per-row
padding or alignment** — unlike ASCII, binary row stride is exactly
`nchannels * sizeof(elem)` bytes; only the very last record of the *entire*
payload is padded (with NaNs, per spec) to a multiple of 80 bytes. This was
confirmed against a real file: row 1 begins at byte offset `342*4 = 1368`
into `belgian_block.crg`'s payload, and `1368 / 80` is not an integer — rows
do not start on fresh 80-byte records the way ASCII rows do.

Any IEEE-754 NaN bit pattern (not one specific sentinel) decodes as NaN,
matching the reference decoder's plain `isnan()` check.

The leftover remainder after floor division must be less than 80 bytes (the
spec's end-of-payload padding quantum) — anything bigger means a truncated
or corrupted file is silently dropping most of a real row, the same
truncation risk `decode_ascii_payload` (Task 7) guards against.
"""
function decode_binary_payload(payload::AbstractVector{UInt8}, format_code::Symbol, nchannels::Int)
    T = BINARY_ELEM_TYPE[format_code]
    esize = sizeof(T)
    row_bytes = nchannels * esize
    nu = length(payload) ÷ row_bytes
    needed = nu * row_bytes
    length(payload) - needed < 80 ||
        error("binary payload has $(length(payload) - needed) leftover bytes after $nu whole rows of $row_bytes bytes each — expected less than 80 bytes of end-of-payload padding, this looks like a truncated row")
    raw = reinterpret(T, Vector{UInt8}(payload[1:needed]))
    raw_be = ntoh.(raw)
    # File layout is row-major (channel fastest within a row); Julia's reshape
    # fills column-major, so reshape(flat, nchannels, nu) puts each file row
    # into one *column* first — transpose to get our (nu, nchannels) convention.
    return Float64.(permutedims(reshape(collect(raw_be), nchannels, nu)))
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_payload.jl")'`
Expected: PASS. If the transposition test fails with rows/columns swapped, that confirms exactly the `reshape`/`permutedims` gotcha flagged above — fix by adjusting which one is transposed, not by reordering the byte read.

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/payload.jl lib/OpenCRG/test/test_payload.jl
git commit -m "OpenCRG: decode binary (KRBI/KDBI) payload rows"
```

---

### Task 9: Channel assembly

**Files:**
- Modify: `lib/OpenCRG/src/payload.jl`
- Modify: `lib/OpenCRG/test/test_payload.jl`

**Step 1: Write the failing test**

Append to `lib/OpenCRG/test/test_payload.jl`:

```julia
@testset "assemble_channels" begin
    r = OpenCRG.parse_road_crg(["REFERENCE_LINE_START_PHI = 0.25"])
    channels = [
        OpenCRG.ChannelDef(:phi, nothing, nothing),
        OpenCRG.ChannelDef(:banking, nothing, nothing),
        OpenCRG.ChannelDef(:long_section, 1.0, nothing),   # declared out of ascending order on purpose
        OpenCRG.ChannelDef(:long_section, -1.0, nothing),
    ]
    raw = [  # 2 rows x 4 channels, in DECLARATION order (phi, banking, v=1.0, v=-1.0)
        99.0  0.1  10.0  20.0
        0.5   0.2  11.0  21.0
    ]
    phi, banking, slope, v, z = OpenCRG.assemble_channels(raw, channels, r)
    @test phi[1] == 0.25            # row-0 placeholder overwritten with REFERENCE_LINE_START_PHI...
    @test phi[2] == 0.5             # ...but row 1 keeps its real stored value
    @test banking == [0.1, 0.2]
    @test slope === nothing
    @test v == [-1.0, 1.0]          # sorted ascending, regardless of declaration order
    @test z == [20.0 10.0; 21.0 11.0]   # columns reordered to match the sorted v
end

@testset "assemble_channels: no long_section channels is an error" begin
    r = OpenCRG.parse_road_crg(String[])
    channels = [OpenCRG.ChannelDef(:phi, nothing, nothing)]
    raw = reshape([0.0, 0.5], 2, 1)
    @test_throws Exception OpenCRG.assemble_channels(raw, channels, r)
end

@testset "assemble_channels: no phi channel declared -- NaN propagates, doesn't error" begin
    r = OpenCRG.parse_road_crg(["REFERENCE_LINE_START_PHI = 0.25"])
    channels = [OpenCRG.ChannelDef(:long_section, 0.0, nothing)]
    raw = reshape([10.0, 11.0], 2, 1)
    phi, banking, slope, v, z = OpenCRG.assemble_channels(raw, channels, r)
    @test phi[1] == 0.25      # row-1 placeholder still gets overwritten...
    @test isnan(phi[2])       # ...but there's no real data to recover the rest from
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_payload.jl")'`
Expected: FAIL — `assemble_channels` not defined.

**Step 3: Write the implementation**

Append to `lib/OpenCRG/src/payload.jl`:

```julia
"""
    assemble_channels(raw, channels, refline) -> (phi, banking, slope, v, z)

Split the raw `(nu, nchannels)` matrix into named channels per their role in
`channels` (from `parse_kd_definition`). Long-section columns are sorted by
ascending `v` — the position-defined form's declaration order isn't
guaranteed ascending by spec, only observed as such in example files.
Row 1's `phi` (1-based; "row 0" in the file/spec's 0-based convention) is
always overwritten with `refline.start_phi`, regardless of what's stored —
the reference implementation does this unconditionally, since that row's
stored value is a documented placeholder, never actually used as an arrival
heading.

Requires at least one `:long_section` channel — a CRG file with none would
have no elevation data at all, which is caught here rather than silently
producing an empty `(nu, 0)` z-matrix that would only fail confusingly much
later (in Task 11+). If no `:phi` channel is declared at all (both vendored
fixtures always declare one, so this is untested against real data),
`phi` stays `fill(NaN, nu)` except for the row-1 placeholder overwrite —
deliberately NOT an error, since NaN already poisons Task 11's integration
loudly rather than silently defaulting to a fabricated straight/`phi=0`
road.
"""
function assemble_channels(raw::Matrix{Float64}, channels::Vector{ChannelDef}, refline::ReferenceLineParams)
    nu = size(raw, 1)
    phi = fill(NaN, nu)
    banking = nothing
    slope = nothing
    long_idxs = Int[]
    for (i, c) in enumerate(channels)
        if c.kind == :phi
            phi = raw[:, i]
        elseif c.kind == :banking
            banking = raw[:, i]
        elseif c.kind == :slope
            slope = raw[:, i]
        elseif c.kind == :long_section
            push!(long_idxs, i)
        end
    end
    isempty(long_idxs) && error("no :long_section channels declared in \$KD_DEFINITION -- a CRG file must have elevation data")
    v_all = v_axis(channels, refline)
    order = sortperm(v_all)
    v = v_all[order]
    z = raw[:, long_idxs[order]]
    isempty(phi) || (phi = copy(phi); phi[1] = refline.start_phi)
    return phi, banking, slope, v, z
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_payload.jl")'`
Expected: PASS

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/payload.jl lib/OpenCRG/test/test_payload.jl
git commit -m "OpenCRG: assemble parsed rows into named channels + sorted v-axis"
```

---

### Task 10: `CRGData` struct and `read_crg`

**Files:**
- Create: `lib/OpenCRG/src/read.jl`
- Create: `lib/OpenCRG/test/test_read.jl`

**Step 1: Write the failing test**

```julia
# lib/OpenCRG/test/test_read.jl
@testset "read_crg" begin
    @testset "small ASCII file end-to-end" begin
        data = OpenCRG.read_crg(joinpath(DATA, "handmade_curved_minimalist.crg"))
        @test data isa OpenCRG.CRGData
        @test data.format_code == :LRFI
        @test length(data.phi) == 23      # derived from the payload's actual row count — this file has no REFERENCE_LINE_END_U at all
        @test size(data.z, 1) == 23
        @test size(data.z, 2) == length(data.v)
    end

    @testset "banked/sloped ASCII file has banking and slope channels" begin
        data = OpenCRG.read_crg(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        @test data.banking !== nothing
        @test data.slope !== nothing
        @test length(data.v) == 7
    end

    @testset "real binary file end-to-end" begin
        data = OpenCRG.read_crg(joinpath(DATA, "belgian_block.crg"))
        @test data.format_code == :KRBI
        @test size(data.z) == (1001, 341)   # 342 channels total minus 1 phi channel
        @test !isempty(data.opts)   # this file has a real $ROAD_CRG_OPTS section -- guards against a swapped section-key typo
    end

    @testset "comment/opts/mods/mpro are wired to the right sections, not silently swapped" begin
        # A prior review caught (via mutation testing) that no existing test would
        # notice if e.g. `opts`/`mods` were accidentally wired to each other's
        # section keys inside read_crg -- these two checks close that gap.
        data = OpenCRG.read_crg(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        @test data.mods isa OpenCRG.RoadCrgMods   # present even when $ROAD_CRG_MODS is absent -- all-nothing default, not a crash
        @test all(f -> getfield(data.mods, f) === nothing, fieldnames(OpenCRG.RoadCrgMods))   # none of these 3 fixtures declare $ROAD_CRG_MODS
    end

    @testset "no channels declared is a clear error, not a DivideError" begin
        path = tempname()
        write(path, "\$CT\nempty file, no \$KD_DEFINITION at all\n\$\$\$\$\n")
        try
            @test_throws Exception OpenCRG.read_crg(path)
        finally
            rm(path)
        end
    end
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_read.jl")'`
Expected: FAIL — `CRGData`/`read_crg` not defined.

**Step 3: Write the implementation**

```julia
# lib/OpenCRG/src/read.jl

"""
A fully parsed OpenCRG file: header metadata plus raw channels. No
interpolation/evaluation lives here — see `transform.jl`'s
`road_surface_grid` for the batched forward transform into world-frame
grids.
"""
struct CRGData
    comment::String
    refline::ReferenceLineParams
    format_code::Symbol
    opts::Dict{String,Float64}
    mods::RoadCrgMods
    mpro::Dict{String,String}
    phi::Vector{Float64}
    banking::Union{Vector{Float64},Nothing}
    slope::Union{Vector{Float64},Nothing}
    v::Vector{Float64}
    z::Matrix{Float64}
end

"""
    read_crg(path) -> CRGData

Parse an ASAM OpenCRG `.crg` file at `path`. The number of longitudinal
grid rows is never computed from `\$ROAD_CRG`'s `start_u`/`end_u`/
`increment` — those fields are sometimes absent entirely (e.g. a
"minimalist" file that only sets `REFERENCE_LINE_INCREMENT`) and are only
descriptive/redundant metadata when present. The true row count comes from
`decode_ascii_payload`/`decode_binary_payload`, which derive it directly
from the payload's actual size (see Tasks 7-8).

Errors immediately, with a clear message, if no channels are declared at
all (a missing or empty `\$KD_DEFINITION` — the sign of a non-CRG file, or
one with no data at all) — without this check, `nchannels == 0` would reach
`decode_ascii_payload`/`decode_binary_payload` and fail with a bare, cryptic
`DivideError` instead.
"""
function read_crg(path::AbstractString)
    bytes = read(path)
    header_end = find_header_end(bytes)
    header_lines = split_lines(bytes[1:header_end-1])
    payload = bytes[header_end:end]
    sections = group_sections(header_lines)

    comment = join(get(sections, "CT", String[]), " ")
    refline = parse_road_crg(get(sections, "ROAD_CRG", String[]))
    format_code, channels = parse_kd_definition(get(sections, "KD_DEFINITION", String[]))
    opts = parse_keyvalues(get(sections, "ROAD_CRG_OPTS", String[]))
    mods = parse_road_crg_mods(get(sections, "ROAD_CRG_MODS", String[]))
    mpro = parse_keyvalue_strings(get(sections, "ROAD_CRG_MPRO", String[]))

    nchannels = length(channels)
    nchannels == 0 && error("no channels declared in \$KD_DEFINITION (section missing or empty) -- not a valid CRG file")
    raw = if format_code in (:LRFI, :LDFI)
        decode_ascii_payload(split_lines(payload), format_code, nchannels)
    else
        decode_binary_payload(payload, format_code, nchannels)
    end
    phi, banking, slope, v, z = assemble_channels(raw, channels, refline)

    return CRGData(comment, refline, format_code, opts, mods, mpro, phi, banking, slope, v, z)
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_read.jl")'`
Expected: PASS

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/read.jl lib/OpenCRG/test/test_read.jl
git commit -m "OpenCRG: CRGData struct and top-level read_crg"
```

---

### Task 11: Reference-line `(x,y)` integration

**Files:**
- Create: `lib/OpenCRG/src/transform.jl`
- Create: `lib/OpenCRG/test/test_transform.jl`
- Create: `lib/OpenCRG/test/data/synthetic_end_anchored.crg` (hand-built fixture for the backward+blend case)

**Step 1: Write the failing test**

First, hand-build a tiny synthetic fixture that exercises the end-anchored (Case B) path — a straight line (`phi` constant `0.0`) where the declared `REFERENCE_LINE_END_X` deliberately does *not* match what simple forward integration would give, so the blend correction is forced to do something non-trivial and checkable by hand:

```
$CT
synthetic straight-line fixture for testing end-anchored integration
$
$ROAD_CRG
REFERENCE_LINE_START_U   = 0.0
REFERENCE_LINE_END_U     = 4.0
REFERENCE_LINE_INCREMENT = 1.0
REFERENCE_LINE_START_X   = 0.0
REFERENCE_LINE_START_Y   = 0.0
REFERENCE_LINE_START_PHI = 0.0
REFERENCE_LINE_END_X     = 4.4
REFERENCE_LINE_END_Y     = 0.0
LONG_SECTION_V_RIGHT     = 0.0
LONG_SECTION_V_LEFT      = 0.0
$
$KD_Definition
#:LRFI
D:reference line phi,rad
D:long section at v = 0.000,m
$
$$$$$$$$10$$$$$$$$20$$$$$$$$30$$$$$$$$40$$$$$$$$50$$$$$$$$60$$$$$$$$70$$$$$$$$80
**unused** 0.0000000
 0.0000000 0.0000000
 0.0000000 0.0000000
 0.0000000 0.0000000
 0.0000000 0.0000000
```

Save this as `lib/OpenCRG/test/data/synthetic_end_anchored.crg` (5 rows: u = 0,1,2,3,4; phi = 0 throughout, i.e. dead straight; but the declared end x is 4.4, not the 4.0 that simple forward integration of a straight line would give — forcing 0.4m of accumulated "error" to be redistributed linearly across the 4 segments).

```julia
# lib/OpenCRG/test/test_transform.jl
@testset "integrate_reference_line" begin
    @testset "Case A: no end anchoring, simple forward Euler" begin
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0", "REFERENCE_LINE_START_X=0.0", "REFERENCE_LINE_START_Y=0.0"])
        phi = [0.0, 0.0, 0.0, 0.0]   # straight line along +x; phi[1] unused
        x, y = OpenCRG.integrate_reference_line(r, phi)
        @test x ≈ [0.0, 1.0, 2.0, 3.0]
        @test y ≈ [0.0, 0.0, 0.0, 0.0]
    end

    @testset "Case B: end-anchored, error redistributed linearly" begin
        data = OpenCRG.read_crg(joinpath(DATA, "synthetic_end_anchored.crg"))
        x, y = OpenCRG.integrate_reference_line(data.refline, data.phi)
        @test x[1] ≈ 0.0
        @test x[end] ≈ 4.4          # true end position is hit exactly
        # 0.4m of error over 4 segments, redistributed with fraction i/(n-1):
        # x[i+1] = (1-frac)*(x[i]+du) + frac*xb[i+1]; verify node 3 (index 3, u=2) by hand:
        @test x[3] ≈ 2.1333333333333333 atol=1e-12
        @test all(y .≈ 0.0)
    end
end
```

(The hand-computed expected value for `x[3]` should be double-checked by actually running the backward pass by hand or in a scratch REPL once the implementation exists — don't take the number above on faith; if it disagrees with your own hand calculation from the formula in Step 3, trust your calculation and fix the test.)

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_transform.jl")'`
Expected: FAIL — `integrate_reference_line` not defined, and `synthetic_end_anchored.crg` doesn't exist yet (create it first, per Step 1).

**Step 3: Write the implementation**

```julia
# lib/OpenCRG/src/transform.jl

"""
    integrate_reference_line(refline, phi) -> (x, y)

Integrate heading `phi` (`phi[i]` = heading of the segment ARRIVING at node
`i`; `phi[1]` is unused as an arrival heading — it's `REFERENCE_LINE_START_PHI`,
already substituted in by `assemble_channels`) into world-frame positions,
replicating `calcRefLine` in the reference implementation's `crgLoader.c`:

- No end position given (`refline.end_x === nothing`): simple forward Euler,
  `x[i+1] = x[i] + du*cos(phi[i+1])`.
- Both `end_x`/`end_y` given: integrate backward from the true end point
  using the same arrival-phi convention, then blend the forward-continuation
  value with the backward-derived value linearly across the whole line
  (fraction 0 at the start, 1 at the end) — redistributing integration error
  instead of leaving it all as a discontinuity at the very end.
"""
function integrate_reference_line(refline::ReferenceLineParams, phi::Vector{Float64})
    nu = length(phi)
    du = refline.increment
    x, y = Vector{Float64}(undef, nu), Vector{Float64}(undef, nu)
    x[1], y[1] = refline.start_x, refline.start_y
    if refline.end_x === nothing || refline.end_y === nothing
        for i in 1:nu-1
            x[i+1] = x[i] + du * cos(phi[i+1])
            y[i+1] = y[i] + du * sin(phi[i+1])
        end
        return x, y
    end
    xb, yb = Vector{Float64}(undef, nu), Vector{Float64}(undef, nu)
    xb[nu], yb[nu] = refline.end_x, refline.end_y
    for i in nu:-1:2
        xb[i-1] = xb[i] - du * cos(phi[i])
        yb[i-1] = yb[i] - du * sin(phi[i])
    end
    for i in 1:nu-1
        fraction = i / (nu - 1)
        x[i+1] = (1 - fraction) * (x[i] + du * cos(phi[i+1])) + fraction * xb[i+1]
        y[i+1] = (1 - fraction) * (y[i] + du * sin(phi[i+1])) + fraction * yb[i+1]
    end
    return x, y
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_transform.jl")'`
Expected: PASS (adjust the hand-computed `x[3]` literal in the test if your own hand calculation gives a different number — see the caveat in Step 1)

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/transform.jl lib/OpenCRG/test/test_transform.jl lib/OpenCRG/test/data/synthetic_end_anchored.crg
git commit -m "OpenCRG: reference-line (x,y) integration, forward and end-anchored"
```

---

### Task 12: Reference-line `z` integration (slope)

**Files:**
- Modify: `lib/OpenCRG/src/transform.jl`
- Modify: `lib/OpenCRG/test/test_transform.jl`

**Step 1: Write the failing test**

Append to `lib/OpenCRG/test/test_transform.jl`:

```julia
@testset "integrate_reference_z" begin
    @testset "no slope channel, zero start_z: early-out to all zeros" begin
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0"])
        @test OpenCRG.integrate_reference_z(r, nothing, 4) == zeros(4)
    end

    @testset "constant slope (no channel, nonzero REFERENCE_LINE_START_S)" begin
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0", "REFERENCE_LINE_START_S=0.1"])
        z_ref = OpenCRG.integrate_reference_z(r, nothing, 4)
        @test z_ref ≈ [0.0, 0.1, 0.2, 0.3]
    end

    @testset "per-row slope channel, no end anchoring" begin
        r = OpenCRG.parse_road_crg(["REFERENCE_LINE_INCREMENT=1.0", "REFERENCE_LINE_START_Z=1.0"])
        slope = [0.0, 0.1, 0.2, 0.3]   # slope[1] unused, matching the phi convention
        z_ref = OpenCRG.integrate_reference_z(r, slope, 4)
        @test z_ref ≈ [1.0, 1.1, 1.3, 1.6]
    end
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_transform.jl")'`
Expected: FAIL — `integrate_reference_z` not defined.

**Step 3: Write the implementation**

Append to `lib/OpenCRG/src/transform.jl`:

```julia
"""
    integrate_reference_z(refline, slope, nu) -> z_ref

1D analogue of `integrate_reference_line`, integrating longitudinal `slope`
(dz/du) into a reference-line elevation profile. If there's no `slope`
channel at all, `refline.start_slope` (default 0.0, from
`REFERENCE_LINE_START_S`) is used as a constant slope for every step. If
there's no slope channel AND `start_slope == 0.0`, the reference
implementation (`calcRefLineZ`) skips this early — we represent that as a
constant-zero vector, since `z`-grid values are always added on top
downstream regardless (see `assemble_z_grid`).
"""
function integrate_reference_z(refline::ReferenceLineParams, slope::Union{Vector{Float64},Nothing}, nu::Int)
    if slope === nothing && refline.start_slope == 0.0
        return zeros(nu)
    end
    du = refline.increment
    slope_at(i) = slope === nothing ? refline.start_slope : slope[i]
    z_ref = Vector{Float64}(undef, nu)
    z_ref[1] = refline.start_z
    if refline.end_z === nothing
        for i in 1:nu-1
            z_ref[i+1] = z_ref[i] + slope_at(i+1) * du
        end
        return z_ref
    end
    zb = Vector{Float64}(undef, nu)
    zb[nu] = refline.end_z
    for i in nu:-1:2
        zb[i-1] = zb[i] - slope_at(i) * du
    end
    for i in 1:nu-1
        fraction = i / (nu - 1)
        z_ref[i+1] = (1 - fraction) * (z_ref[i] + slope_at(i+1) * du) + fraction * zb[i+1]
    end
    return z_ref
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_transform.jl")'`
Expected: PASS

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/transform.jl lib/OpenCRG/test/test_transform.jl
git commit -m "OpenCRG: reference-line z (elevation) integration from slope"
```

---

### Task 13: Lateral offset grid (miter-normal)

This is the trickiest geometry in the whole plan. Read the docstring below carefully before touching the code — the formula is derived from the reference implementation's `crgEvaluv2xy.c`, not from first-principles guessing, and its correctness is only fully confirmed by Task 17's comprehensive cross-validation, not by this task's smaller unit tests alone.

**Files:**
- Modify: `lib/OpenCRG/src/transform.jl`
- Modify: `lib/OpenCRG/test/test_transform.jl`

**Step 1: Write the failing test**

Append to `lib/OpenCRG/test/test_transform.jl`:

```julia
@testset "lateral_offset_grid" begin
    @testset "straight line: pure perpendicular offset" begin
        x = [0.0, 1.0, 2.0, 3.0]   # straight along +x
        y = [0.0, 0.0, 0.0, 0.0]
        v = [-1.0, 0.0, 1.0]
        X, Y = OpenCRG.lateral_offset_grid(x, y, v)
        @test X ≈ repeat(x, 1, 3)                     # offsetting perpendicular to +x doesn't move x
        @test Y ≈ [-1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0] || Y ≈ -[-1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0; -1.0 0.0 1.0]
        # (whichever sign convention `perp` uses; either is "correct" in isolation —
        # Task 17's cross-validation against the C library is the real arbiter.)
    end

    @testset "shape sanity: v=0 reproduces the reference line exactly" begin
        x = [0.0, 1.0, 2.5, 2.5]   # includes a kink, to exercise the miter-normal interior formula
        y = [0.0, 0.5, 1.0, 2.0]
        X, Y = OpenCRG.lateral_offset_grid(x, y, [0.0])
        @test X[:, 1] ≈ x
        @test Y[:, 1] ≈ y
    end
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_transform.jl")'`
Expected: FAIL — `lateral_offset_grid` not defined.

**Step 3: Write the implementation**

```julia
"""
    lateral_offset_grid(x, y, v) -> (X, Y)

Offset the reference line `(x,y)` laterally by each `v` to build the full
`(nu, nv)` world-frame grid, replicating the reference implementation's
"miter normal" scheme (`crgEvaluv2xy.c`) evaluated exactly at grid nodes.

At each INTERIOR node `i`, the true C function (queried continuously between
nodes) computes an offset point by bisecting the incoming/outgoing segment
directions via the chord skipping over node `i` (from node `i-1` directly to
node `i+1`), then rescales that bisector so its projection onto the
FOLLOWING segment's own normal reproduces the true perpendicular distance
`v` (a standard "miter join" construction — this is what keeps offset
polylines from gapping/overlapping at kinks in the reference line). At
`u = u_i` exactly, only this one node's offset point is used (no
interpolation between the segment's two endpoints is needed, since
evaluating exactly at a grid node has interpolation fraction 0). The first
and last nodes have no bisector partner and fall back to their single
adjacent segment's own normal.

Perpendicular convention here: `perp(dx, dy) = (-dy, dx)` (90° CCW). If
Task 17's full-grid cross-validation shows every `v` mismatched by a sign
flip (i.e. `X`/`Y` match with `v` negated), that's this convention being
mirrored relative to the C library — fix by negating `perp`'s output here,
not by negating `v` itself, since `v`'s sign also matters for the banking
term in `assemble_z_grid` (Task 14) and must stay consistent with the
parsed `v` axis.
"""
function lateral_offset_grid(x::Vector{Float64}, y::Vector{Float64}, v::Vector{Float64})
    nu = length(x)
    perp(d) = (-d[2], d[1])
    normalize2(d) = (n = hypot(d[1], d[2]); (d[1]/n, d[2]/n))

    seg = [normalize2((x[i+1]-x[i], y[i+1]-y[i])) for i in 1:nu-1]
    n12 = perp.(seg)

    offset_dir = Vector{Tuple{Float64,Float64}}(undef, nu)
    if nu == 1
        offset_dir[1] = (0.0, 1.0)   # degenerate single-node "line"; arbitrary but consistent
    else
        offset_dir[1] = n12[1]
        offset_dir[nu] = n12[nu-1]
        for i in 2:nu-1
            chord = normalize2((x[i+1]-x[i-1], y[i+1]-y[i-1]))
            n1 = perp(chord)
            denom = n1[1]*n12[i][1] + n1[2]*n12[i][2]
            offset_dir[i] = (n1[1]/denom, n1[2]/denom)
        end
    end

    nv = length(v)
    X, Y = Matrix{Float64}(undef, nu, nv), Matrix{Float64}(undef, nu, nv)
    for i in 1:nu, j in 1:nv
        X[i,j] = x[i] + v[j] * offset_dir[i][1]
        Y[i,j] = y[i] + v[j] * offset_dir[i][2]
    end
    return X, Y
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_transform.jl")'`
Expected: PASS

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/transform.jl lib/OpenCRG/test/test_transform.jl
git commit -m "OpenCRG: lateral miter-normal offset grid (x,y)"
```

---

### Task 14: Z grid assembly (elevation + banking)

**Files:**
- Modify: `lib/OpenCRG/src/transform.jl`
- Modify: `lib/OpenCRG/test/test_transform.jl`

**Step 1: Write the failing test**

Append to `lib/OpenCRG/test/test_transform.jl`:

```julia
@testset "assemble_z_grid" begin
    r = OpenCRG.parse_road_crg(String[])
    z_grid = [0.0 0.1 0.2; 1.0 1.1 1.2]   # 2 rows x 3 v-columns
    z_ref = [10.0, 20.0]
    v = [-1.0, 0.0, 1.0]

    @testset "no banking: pure additive z_grid + z_ref" begin
        Z = OpenCRG.assemble_z_grid(z_grid, z_ref, nothing, r, v)
        @test Z ≈ [10.0 10.1 10.2; 21.0 21.1 21.2]
    end

    @testset "banking adds v * bank(u), clipped to [v_min, v_max]" begin
        banking = [0.05, -0.05]
        Z = OpenCRG.assemble_z_grid(z_grid, z_ref, banking, r, v)
        @test Z[1, :] ≈ [10.0 - 0.05, 10.1, 10.2 + 0.05]
        @test Z[2, :] ≈ [21.0 + 0.05, 21.1, 21.2 - 0.05]
    end

    @testset "banking's v-reach is clipped to the road's actual width" begin
        banking = [0.1, 0.1]
        wide_v = [-5.0, 0.0, 5.0]   # beyond [v[1], v[end]] = [-1,1] in this synthetic case
        Z = OpenCRG.assemble_z_grid(z_grid, z_ref, banking, r, wide_v)
        @test Z[1, 1] ≈ z_grid[1,1] + z_ref[1] + 0.1 * (-1.0)   # clipped to v[1], not -5.0
        @test Z[1, 3] ≈ z_grid[1,3] + z_ref[1] + 0.1 * 1.0      # clipped to v[end], not 5.0
    end
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_transform.jl")'`
Expected: FAIL — `assemble_z_grid` not defined.

**Step 3: Write the implementation**

Append to `lib/OpenCRG/src/transform.jl`:

```julia
"""
    assemble_z_grid(z_grid, z_ref, banking, refline, v) -> Z

`Z[i,j] = z_grid[i,j] + z_ref[i] + bank(i) * clamp(v[j], v[1], v[end])`.
Banking's lateral reach is clamped to the road's actual v-range regardless
of what `v` values are being queried — matching `crgEvalz.c`, which clips
`v` for the banking term specifically (but not for the grid-z lookup
itself). If there's no `banking` channel, `refline.start_banking` (default
0.0, from `REFERENCE_LINE_START_B`) is used as a constant cross-slope for
every row.
"""
function assemble_z_grid(z_grid::Matrix{Float64}, z_ref::Vector{Float64}, banking::Union{Vector{Float64},Nothing}, refline::ReferenceLineParams, v::Vector{Float64})
    nu, nv = size(z_grid)
    Z = Matrix{Float64}(undef, nu, nv)
    vmin, vmax = first(v), last(v)
    for i in 1:nu
        bank_i = banking === nothing ? refline.start_banking : banking[i]
        for j in 1:nv
            vc = clamp(v[j], vmin, vmax)
            Z[i,j] = z_grid[i,j] + z_ref[i] + bank_i * vc
        end
    end
    return Z
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_transform.jl")'`
Expected: PASS

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/transform.jl lib/OpenCRG/test/test_transform.jl
git commit -m "OpenCRG: assemble Z grid from elevation + z_ref + clipped banking"
```

---

### Task 15: Extend `lib/LibOpenCRG` with modifier bindings

The MODS cross-validation in Task 16 needs the C library's own modifier-setting/apply functions, which weren't in scope when `lib/LibOpenCRG` was first built. This task only touches `lib/LibOpenCRG`, not `lib/OpenCRG`.

**Files:**
- Modify: `lib/LibOpenCRG/src/LibOpenCRG.jl`

**Step 1: Write the failing test**

There's no separate test file for this — add directly to `lib/LibOpenCRG/test/runtests.jl` (append a new `@testset`):

```julia
@testset "modifier bindings" begin
    dsId = crgLoaderReadFile(joinpath(@__DIR__, "data", "handmade_curved_minimalist.crg"))
    @test dsId != 0
    @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineOffsetX, 100.0) != 0
    @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineOffsetPhi, 1.57) != 0
    crgDataSetModifiersApply(dsId)   # returns Nothing; success is implicit if it doesn't crash

    cpId = crgContactPointCreate(dsId)
    r = crgEvaluv2xy(cpId, 0.0, 0.0)
    @test r.status != 0
    @test r.x ≈ 100.0 atol=1e-6   # the reference line's own start point, after +100 in x

    crgContactPointDelete(cpId)
    crgDataSetRelease(dsId)
    crgMemRelease()
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/LibOpenCRG lib/LibOpenCRG/test/runtests.jl`
Expected: FAIL — `crgDataSetModifierSetDouble`/`crgDataSetModifiersApply`/`dCrgModRefLineOffsetX`/`dCrgModRefLineOffsetPhi` not defined.

**Step 3: Write the implementation**

In `lib/LibOpenCRG/src/LibOpenCRG.jl`, add the new names to the existing `export` list (alongside the current ones):

```julia
       crgDataSetModifierSetInt,
       crgDataSetModifierSetDouble,
       crgDataSetModifierRemoveAll,
       crgDataSetModifiersApply,
```

Then append a new section near the other `crgMgr.c`-sourced functions:

```julia
# =====================================================================
# ====== METHODS in crgMgr.c (modifiers) ======
# =====================================================================

# Modifier IDs, from crgBaseLib.h (values MUST match exactly — these select
# which field a Set*/Get* call targets).
const dCrgModScaleZ            = 21
const dCrgModScaleSlope        = 22
const dCrgModScaleBank         = 23
const dCrgModScaleLength       = 24
const dCrgModScaleWidth        = 25
const dCrgModScaleCurvature    = 26
const dCrgModGridNaNMode       = 27
const dCrgModGridNaNOffset     = 28
const dCrgModRefPointV         = 29
const dCrgModRefPointVFrac     = 30
const dCrgModRefPointVOffset   = 31
const dCrgModRefPointU         = 32
const dCrgModRefPointUFrac     = 33
const dCrgModRefPointUOffset   = 34
const dCrgModRefPointX         = 35
const dCrgModRefPointY         = 36
const dCrgModRefPointZ         = 37
const dCrgModRefPointPhi       = 38
const dCrgModRefLineOffsetX    = 39
const dCrgModRefLineOffsetY    = 40
const dCrgModRefLineOffsetZ    = 41
const dCrgModRefLineOffsetPhi  = 42
const dCrgModRefLineRotCenterX = 43
const dCrgModRefLineRotCenterY = 44

"""
    crgDataSetModifierSetInt(dataSetId, optionId, optionValue) -> Cint

Set/add an integer-valued modifier (e.g. `dCrgModGridNaNMode`), to be
applied the next time `crgDataSetModifiersApply` is called. Returns 1 on
success.
"""
function crgDataSetModifierSetInt(dataSetId::Integer, optionId::Integer, optionValue::Integer)
    return ccall((:crgDataSetModifierSetInt, LIBOPENCRG_PATH), Cint,
                 (Cint, Cuint, Cint), dataSetId, optionId, optionValue)
end

"""
    crgDataSetModifierSetDouble(dataSetId, optionId, optionValue) -> Cint

Set/add a double-valued modifier (e.g. `dCrgModRefLineOffsetPhi`), to be
applied the next time `crgDataSetModifiersApply` is called. Returns 1 on
success.
"""
function crgDataSetModifierSetDouble(dataSetId::Integer, optionId::Integer, optionValue::Real)
    return ccall((:crgDataSetModifierSetDouble, LIBOPENCRG_PATH), Cint,
                 (Cint, Cuint, Cdouble), dataSetId, optionId, optionValue)
end

"""
    crgDataSetModifierRemoveAll(dataSetId) -> Cint

Remove all modifiers from the data set's modifier list without applying
them. Returns 1 on success.
"""
function crgDataSetModifierRemoveAll(dataSetId::Integer)
    return ccall((:crgDataSetModifierRemoveAll, LIBOPENCRG_PATH), Cint, (Cint,), dataSetId)
end

"""
    crgDataSetModifiersApply(dataSetId) -> Nothing

Apply all currently-set modifiers to the data set once, then clear them
from the modifier list.
"""
function crgDataSetModifiersApply(dataSetId::Integer)
    ccall((:crgDataSetModifiersApply, LIBOPENCRG_PATH), Cvoid, (Cint,), dataSetId)
    return nothing
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/LibOpenCRG lib/LibOpenCRG/test/runtests.jl`
Expected: PASS (all testsets, including the new one)

**Step 5: Commit**

```bash
git add lib/LibOpenCRG/src/LibOpenCRG.jl lib/LibOpenCRG/test/runtests.jl
git commit -m "LibOpenCRG: bind modifier-setting functions for MODS cross-validation"
```

---

### Task 16: Apply `$ROAD_CRG_MODS`

**Files:**
- Modify: `lib/OpenCRG/src/transform.jl`
- Modify: `lib/OpenCRG/test/test_transform.jl`
- Modify: `lib/OpenCRG/Project.toml` (no change needed — `LibOpenCRG` already added as a test dep in Task 1)

**Step 1: Write the failing test**

Append to `lib/OpenCRG/test/test_transform.jl`:

```julia
using LibOpenCRG

@testset "apply_mods" begin
    @testset "no mods: identity (road_surface_grid unaffected)" begin
        data = OpenCRG.read_crg(joinpath(DATA, "handmade_curved_minimalist.crg"))
        d2 = OpenCRG.apply_mods(data)
        @test d2.phi == data.phi
        @test d2.refline.start_x == data.refline.start_x
    end

    @testset "REFLINE_OFFSET_*: rotate then translate, cross-validated against LibOpenCRG" begin
        path = joinpath(DATA, "handmade_curved_minimalist.crg")
        data = OpenCRG.read_crg(path)
        mods = OpenCRG.RoadCrgMods(refline_offset_phi=1.57, refline_offset_x=100.0, refline_offset_y=50.0)
        data_with_mods = OpenCRG.CRGData(data.comment, data.refline, data.format_code, data.opts, mods,
                                          data.mpro, data.phi, data.banking, data.slope, data.v, data.z)
        u, v, X, Y, Z = OpenCRG.road_surface_grid(data_with_mods)

        dsId = crgLoaderReadFile(path)
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineOffsetPhi, 1.57) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineOffsetX, 100.0) != 0
        @test crgDataSetModifierSetDouble(dsId, LibOpenCRG.dCrgModRefLineOffsetY, 50.0) != 0
        crgDataSetModifiersApply(dsId)
        cpId = crgContactPointCreate(dsId)
        for i in eachindex(u), j in eachindex(v)
            ref = crgEvaluv2xy(cpId, u[i], v[j])
            @test ref.status != 0
            @test X[i,j] ≈ ref.x atol=1e-6
            @test Y[i,j] ≈ ref.y atol=1e-6
        end
        crgContactPointDelete(cpId); crgDataSetRelease(dsId); crgMemRelease()
    end

    @testset "SCALE_Z_GRID doubles elevation" begin
        data = OpenCRG.read_crg(joinpath(DATA, "handmade_curved_banked_sloped.crg"))
        mods = OpenCRG.RoadCrgMods(scale_z_grid=2.0)
        data_with_mods = OpenCRG.CRGData(data.comment, data.refline, data.format_code, data.opts, mods,
                                          data.mpro, data.phi, data.banking, data.slope, data.v, data.z)
        d2 = OpenCRG.apply_mods(data_with_mods)
        @test d2.z ≈ 2.0 .* data.z
    end
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_transform.jl")'`
Expected: FAIL — `apply_mods` not defined (and `road_surface_grid`, added in Task 17, not yet either — if Task 17 hasn't landed yet when running this in isolation, temporarily stub `road_surface_grid` or run this task's tests after Task 17; the plan presents them in this order for narrative clarity, but Tasks 16 and 17 are tightly coupled and may be easiest to implement together in one sitting).

**Step 3: Write the implementation**

Append to `lib/OpenCRG/src/transform.jl`:

```julia
"""
    apply_mods(data::CRGData) -> CRGData

Apply `\$ROAD_CRG_MODS`, returning a new `CRGData` with adjusted raw
channels. Order matches the reference implementation
(`crgDataSetModifiersApply` in `crgMgr.c`): scale channels first (z-grid,
slope, banking, length, width, curvature), then a single rotate+translate.

Key insight that keeps this simple: rotating every `phi` value by a
constant angle rotates the *shape* of the eventually-integrated reference
line by that same angle, because `(cos(a+θ), sin(a+θ)) = R(θ)·(cos a, sin a)`
— so only `refline.start_x/start_y` (rotated about the pivot, then
translated) and `phi` itself need to change here. There's no need to
integrate the reference line inside this function at all, EXCEPT for the
`REFPOINT_*` case, which needs to evaluate the CURRENT (pre-transform)
position at a specific `(u,v)` to know what point is being pinned down.

`REFPOINT_*` (if ANY such field is set) takes over the rotate+translate
step entirely, ignoring `REFLINE_OFFSET_*`/`REFLINE_ROTCENTER_*` completely
— they do not compose, matching `crgDataApplyTransformations` in `crgMgr.c`.

Deliberately NOT applied: `mods.grid_nan_mode`/`mods.grid_nan_offset` (parsed
by Task 6, real behavior in `crgMgr.c` around line 601) control how NaN gaps
in the z-grid get filled/replaced — a data-cleaning concern, not the
scale/offset/rotate geometry transform this task's scope was explicitly
limited to. Same deferral treatment as `\$ROAD_CRG_MPRO` (Task 5's design
doc note) — parsed and available on `RoadCrgMods`, never auto-applied.
"""
function apply_mods(data::CRGData)
    mods = data.mods
    r = data.refline
    phi = copy(data.phi)
    z = copy(data.z)
    slope = data.slope === nothing ? nothing : copy(data.slope)
    banking = data.banking === nothing ? nothing : copy(data.banking)
    v = copy(data.v)
    u_increment, start_slope, end_slope, start_banking, end_banking =
        r.increment, r.start_slope, r.end_slope, r.start_banking, r.end_banking

    mods.scale_z_grid === nothing || (z .*= mods.scale_z_grid)
    if mods.scale_slope !== nothing
        slope === nothing || (slope .*= mods.scale_slope)
        start_slope *= mods.scale_slope
        end_slope = end_slope === nothing ? nothing : end_slope * mods.scale_slope
    end
    if mods.scale_banking !== nothing
        banking === nothing || (banking .*= mods.scale_banking)
        start_banking *= mods.scale_banking
        end_banking = end_banking === nothing ? nothing : end_banking * mods.scale_banking
    end
    mods.scale_length === nothing || (u_increment *= mods.scale_length)
    mods.scale_width === nothing || (v .*= mods.scale_width)
    if mods.scale_curvature !== nothing
        base = phi[1]
        for i in 2:length(phi)
            phi[i] = base + mods.scale_curvature * (phi[i] - base)
        end
    end

    r2 = ReferenceLineParams(r.start_u, r.end_u, u_increment, r.start_x, r.start_y, r.start_phi,
        r.end_x, r.end_y, r.end_phi, r.start_z, r.end_z, r.v_right, r.v_left, r.v_increment,
        start_slope, end_slope, start_banking, end_banking)

    has_refpoint = any(f -> getfield(mods, f) !== nothing, (:refpoint_u, :refpoint_u_fraction,
        :refpoint_u_offset, :refpoint_v, :refpoint_v_fraction, :refpoint_v_offset,
        :refpoint_x, :refpoint_y, :refpoint_z, :refpoint_phi))

    if has_refpoint
        x0, y0 = integrate_reference_line(r2, phi)
        u_frac, u_off = something(mods.refpoint_u_fraction, 0.0), something(mods.refpoint_u_offset, 0.0)
        u_pos = something(mods.refpoint_u, r2.start_u + u_frac * (r2.end_u - r2.start_u) + u_off)
        idx = clamp(round(Int, (u_pos - r2.start_u) / u_increment) + 1, 1, length(phi))
        from_x, from_y, from_phi = x0[idx], y0[idx], phi[idx]
        rot_center = (from_x, from_y)
        rot_angle = something(mods.refpoint_phi, 0.0) - from_phi
        translation = (something(mods.refpoint_x, 0.0) - from_x, something(mods.refpoint_y, 0.0) - from_y)
    else
        rot_center = (something(mods.refline_rotcenter_x, r2.start_x), something(mods.refline_rotcenter_y, r2.start_y))
        rot_angle = something(mods.refline_offset_phi, 0.0)
        translation = (something(mods.refline_offset_x, 0.0), something(mods.refline_offset_y, 0.0))
    end

    c, s = cos(rot_angle), sin(rot_angle)
    dx, dy = r2.start_x - rot_center[1], r2.start_y - rot_center[2]
    new_start_x = rot_center[1] + dx*c - dy*s + translation[1]
    new_start_y = rot_center[2] + dx*s + dy*c + translation[2]
    new_start_phi = r2.start_phi + rot_angle
    phi .+= rot_angle
    new_start_z = r2.start_z + something(mods.refline_offset_z, 0.0)
    new_end_z = r2.end_z === nothing ? nothing : r2.end_z + something(mods.refline_offset_z, 0.0)
    new_end_x = r2.end_x === nothing ? nothing : rot_center[1] + (r2.end_x-rot_center[1])*c - (r2.end_y-rot_center[2])*s + translation[1]
    new_end_y = r2.end_y === nothing ? nothing : rot_center[2] + (r2.end_x-rot_center[1])*s + (r2.end_y-rot_center[2])*c + translation[2]
    new_end_phi = r2.end_phi === nothing ? nothing : r2.end_phi + rot_angle

    r3 = ReferenceLineParams(r2.start_u, r2.end_u, r2.increment, new_start_x, new_start_y, new_start_phi,
        new_end_x, new_end_y, new_end_phi, new_start_z, new_end_z,
        r2.v_right, r2.v_left, r2.v_increment, r2.start_slope, r2.end_slope, r2.start_banking, r2.end_banking)

    return CRGData(data.comment, r3, data.format_code, data.opts, mods, data.mpro, phi, banking, slope, v, z)
end
```

Note the added `new_end_x/new_end_y/new_end_phi` handling (rotating/translating the anchor point too, when present) — this wasn't in the earlier design sketch and is necessary so that `integrate_reference_line`'s Case B (Task 11) still anchors to the *correctly transformed* end point, not the pre-rotation one, when both end-anchoring and MODS rotation are present in the same file.

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG, LibOpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_transform.jl")'`
Expected: PASS. If the cross-validation sub-testset fails, the mismatch is almost certainly in the rotate/translate math above (sign of `rot_angle`, or `rot_center`/`translation` mixed up) — re-read `crgDataApplyTransformations` in the vendored `lib/LibOpenCRG/csrc/src/crgMgr.c` at the specific point of disagreement rather than guessing a fix.

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/transform.jl lib/OpenCRG/test/test_transform.jl
git commit -m "OpenCRG: apply \$ROAD_CRG_MODS (scale chain + rotate/translate)"
```

---

### Task 17: `road_surface_grid` and full-grid cross-validation

**Files:**
- Modify: `lib/OpenCRG/src/transform.jl`
- Create: `lib/OpenCRG/test/test_crossvalidate.jl`
- Modify: `lib/OpenCRG/Project.toml` (no change — already done in Task 1)

**Step 1: Write the failing test**

```julia
# lib/OpenCRG/test/test_crossvalidate.jl
using LibOpenCRG

@testset "road_surface_grid, cross-validated against the LibOpenCRG oracle" begin
    for fname in ["handmade_curved_minimalist.crg", "handmade_curved_banked_sloped.crg"]
        @testset "$fname" begin
            path = joinpath(DATA, fname)
            data = OpenCRG.read_crg(path)
            u, v, X, Y, Z = OpenCRG.road_surface_grid(data)
            @test size(X) == (length(u), length(v))
            @test size(Y) == size(X)
            @test size(Z) == size(X)

            dsId = crgLoaderReadFile(path)
            @test dsId != 0
            @test crgCheck(dsId) != 0
            cpId = crgContactPointCreate(dsId)
            @test cpId != -1

            mismatches = 0
            for i in eachindex(u), j in eachindex(v)
                ref_xy = crgEvaluv2xy(cpId, u[i], v[j])
                if ref_xy.status != 0
                    if !(X[i,j] ≈ ref_xy.x atol=1e-6) || !(Y[i,j] ≈ ref_xy.y atol=1e-6)
                        mismatches += 1
                    end
                end
                ref_z = crgEvaluv2z(cpId, u[i], v[j])
                if ref_z.status != 0 && !isnan(ref_z.z) && !isnan(Z[i,j])
                    if !(Z[i,j] ≈ ref_z.z atol=1e-6)
                        mismatches += 1
                    end
                end
            end
            @test mismatches == 0

            crgContactPointDelete(cpId)
            crgDataSetRelease(dsId)
            crgMemRelease()
        end
    end
end
```

(Deliberately counting mismatches into one variable rather than asserting inside the loop — with a per-node `@test` inside a double loop over a real file's full grid, a systematic bug produces hundreds of near-identical failure messages that bury the actual pattern. One summary assertion is easier to act on; if it fails, drop a `@show i, j, X[i,j], ref_xy.x, Y[i,j], ref_xy.y` right before the mismatch-counting `if` to see exactly where and how it diverges.)

**Step 2: Run test to verify it fails**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG, LibOpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_crossvalidate.jl")'`
Expected: FAIL — `road_surface_grid` not defined.

**Step 3: Write the implementation**

Append to `lib/OpenCRG/src/transform.jl`:

```julia
"""
    road_surface_grid(data::CRGData) -> (u, v, X, Y, Z)

The batched forward transform: applies `\$ROAD_CRG_MODS` (a no-op if none
are set), integrates the reference line and its elevation profile, and
offsets the whole grid laterally — producing world-frame `(u, v, X, Y, Z)`,
all `(nu, nv)` except the `u`/`v` axis vectors themselves.
"""
function road_surface_grid(data::CRGData)
    d = apply_mods(data)
    nu = length(d.phi)
    u = [d.refline.start_u + (i-1)*d.refline.increment for i in 1:nu]
    x, y = integrate_reference_line(d.refline, d.phi)
    z_ref = integrate_reference_z(d.refline, d.slope, nu)
    X, Y = lateral_offset_grid(x, y, d.v)
    Z = assemble_z_grid(d.z, z_ref, d.banking, d.refline, d.v)
    return u, d.v, X, Y, Z
end
```

Also update `lib/OpenCRG/Project.toml`'s `[deps]`... actually no change needed there (LibOpenCRG is already a test-only `[extras]`/`[targets]` dependency from Task 1) — but do add `using LibOpenCRG` is only inside test files, never inside `src/`, since production `OpenCRG.jl` must stay dependency-free per the design doc.

**Step 4: Run test to verify it passes**

Run: `julia --project=lib/OpenCRG -e 'using Test, OpenCRG, LibOpenCRG; const DATA=joinpath("lib/OpenCRG/test","data"); include("lib/OpenCRG/test/test_crossvalidate.jl")'`

If this fails with every point systematically off by a constant sign flip in `Y` (or `X`) relative to `v`, that confirms the `perp` convention in Task 13 is mirrored relative to the C library — negate `perp`'s output in `lateral_offset_grid` and rerun. If it fails only at the FIRST or LAST `u` node (not throughout), that points at the boundary fallback (`offset_dir[1]`/`offset_dir[nu]`) needing a different formula than "just the adjacent segment's own normal" — re-read `crgEvaluv2xy.c`'s boundary handling directly (`lib/LibOpenCRG/csrc/src/crgEvaluv2xy.c`) at that point rather than guessing further. If it fails throughout but by varying (not constant) amounts, suspect the miter-normal formula's algebra itself (re-derive from the `crgEvaluv2xy.c` source directly, checking the `a`/`b`/`n1`/`n2`/`n12` construction line by line against what's implemented here).

Once green, also run the FULL suite for the whole package: `julia --project=lib/OpenCRG lib/OpenCRG/test/runtests.jl`
Expected: all testsets across all files pass.

**Step 5: Commit**

```bash
git add lib/OpenCRG/src/transform.jl lib/OpenCRG/test/test_crossvalidate.jl
git commit -m "OpenCRG: road_surface_grid top-level transform, cross-validated end-to-end"
```

---

## After all 17 tasks

Run the complete suite once more from a clean Julia session to make sure nothing was left in a stale/order-dependent state:

```bash
julia --project=lib/OpenCRG -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Then consider whether `OpenCRG` should be wired into the root `CarComponents` `dyad/` model (replacing or augmenting `RoadData`'s synthetic sine bump with a real `road_surface_grid` fed through `DataInterpolationsND`/`BlockComponents.Tables.InterpolatedTable`) — that integration is intentionally out of scope for this plan (per the design doc's package-boundary decision) and would be its own follow-up plan.
