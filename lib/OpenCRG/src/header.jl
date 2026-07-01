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
keyword (empty string for a bare `\$` section-end marker). A trailing
inline (`!`) comment on the marker line itself — e.g. `\$ROAD_CRG ! crg road
parameters`, seen in real vendored fixtures — is stripped first, the same
way it already is for ordinary content lines; otherwise the comment text
would leak into the keyword and no section would ever be found by name.
"""
function section_marker(line::AbstractString)
    (isempty(line) || line[1] != '$') && return nothing
    return uppercase(strip(strip_inline_comment(line[2:end])))
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
            if isempty(kw) || all(==('$'), kw)   # bare "$" OR an all-"$" terminator line
                current = ""
                continue
            end
            current = kw
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

"""
    parse_road_crg(lines) -> ReferenceLineParams

Parse `\$ROAD_CRG` section lines into a `ReferenceLineParams`. Deliberately
does NOT extract `REFERENCE_LINE_START_LON`/`_LAT`/`_ALT` or the
corresponding `_END_*` geodetic fields (present in some real files, e.g.
`belgian_block.crg`) — placing the road on a real-world map is
`\$ROAD_CRG_MPRO`'s job (geospatial projection metadata), which this
package parses but never applies; see the design doc.
"""
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
