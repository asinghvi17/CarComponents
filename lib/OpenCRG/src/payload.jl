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
