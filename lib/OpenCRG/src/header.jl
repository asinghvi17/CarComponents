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
