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
all (a missing or empty `\$KD_DEFINITION` -- the sign of a non-CRG file, or
one with no data at all) -- without this check, `nchannels == 0` would reach
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
