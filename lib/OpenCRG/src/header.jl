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
