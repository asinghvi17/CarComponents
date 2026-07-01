# lib/OpenCRG/src/payload.jl

const ASCII_FIELD_WIDTH = Dict(:LRFI => 10, :LDFI => 20)
const FIELDS_PER_RECORD = Dict(:LRFI => 8, :LDFI => 4, :KRBI => 20, :KDBI => 10)
const BINARY_ELEM_TYPE = Dict(:KRBI => Float32, :KDBI => Float64)

"""
    decode_ascii_field(s) -> Float64

Decode one fixed-width ASCII field (already `strip`ped of surrounding
whitespace by the caller). The literal 10-character token `**unused**`
decodes to `0.0` (an exact-string match only — this quirk does not apply to
LDFI's 20-character fields, where a padded `**unused**...` never matches
this exact 10-character literal and correctly falls through to NaN below,
matching the reference C decoder's `strcmp` check). Any other field
containing a character outside `0123456789+-.eEdD` decodes to NaN
(word-agnostic: `*missing*`, `*none*`, etc. all become NaN). Fortran-style
`D`/`d` exponent markers are normalized to `e` before `parse`.
"""
function decode_ascii_field(s::AbstractString)
    s == "**unused**" && return 0.0
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
7 above for why this can't come from `\$ROAD_CRG`'s `start_u`/`end_u`.
Each row always starts on a fresh line — per spec §6.4.8, plain-text rows
never pack a previous row's leftover fields onto the same line, unlike
binary (Task 8).
"""
function decode_ascii_payload(payload_lines::Vector{<:AbstractString}, format_code::Symbol, nchannels::Int)
    width = ASCII_FIELD_WIDTH[format_code]
    per_record = FIELDS_PER_RECORD[format_code]
    lines_per_row = cld(nchannels, per_record)
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
                data[row, ch] = decode_ascii_field(strip(field))
            end
        end
    end
    return data
end
