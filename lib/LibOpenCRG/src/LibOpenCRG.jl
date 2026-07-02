"""
    LibOpenCRG

Thin `ccall` wrapper around the upstream ASAM OpenCRG C reference
implementation ("baselib"), vendored and compiled locally from source (see
`csrc/` and `build.jl`).

This package exists **only** to serve as a test oracle for cross-validating
the pure-Julia `OpenCRG.jl` package. It is not intended for use in production
vehicle models: the pure-Julia implementation is the primary implementation.

Before this module can be loaded, the vendored C sources must be compiled
into a shared library by running:

    julia lib/LibOpenCRG/build.jl

which produces `lib/LibOpenCRG/build/libLibOpenCRG.<dylib|so|dll>`.

Function names and signatures mirror `csrc/inc/crgBaseLib.h` exactly. Where
the C function has pointer "out" parameters, the Julia wrapper allocates the
`Ref` internally and returns a `NamedTuple` of `(status, <outputs...>)`
instead of requiring the caller to pass in `Ref`s.
"""
module LibOpenCRG

export crgDataSetRelease,
       crgDataPrintHeader,
       crgDataPrintChannelInfo,
       crgDataPrintRoadInfo,
       crgDataSetGetURange,
       crgDataSetGetVRange,
       crgDataSetGetIncrements,
       crgDataSetGetUtilityDataClosedTrack,
       crgMemRelease,
       crgGetReleaseInfo,
       crgMsgSetLevel,
       crgCheck,
       crgLoaderReadFile,
       crgLoaderSuppressFileNotFoundFatalMsg,
       crgContactPointCreate,
       crgContactPointDelete,
       crgContactPointDeleteAll,
       crgEvalxy2uv,
       crgEvaluv2xy,
       crgEvaluv2z,
       crgEvalxy2z,
       crgEvaluv2pk,
       crgEvalxy2pk,
       crgDataSetModifierSetInt,
       crgDataSetModifierSetDouble,
       crgDataSetModifierRemoveAll,
       crgDataSetModifiersApply

const PKG_DIR = joinpath(@__DIR__, "..")

const LIB_EXT = Sys.isapple() ? "dylib" : Sys.islinux() ? "so" : Sys.iswindows() ? "dll" : "so"

"""
Path to the compiled shared library, resolved relative to the package
directory at load time. Built by `build.jl`.
"""
const LIBOPENCRG_PATH = joinpath(PKG_DIR, "build", "libLibOpenCRG.$(LIB_EXT)")

function __init__()
    if !isfile(LIBOPENCRG_PATH)
        error("""
        LibOpenCRG: shared library not found at:
            $(LIBOPENCRG_PATH)

        It must be compiled from the vendored C sources before this module \
        can be used. Build it with:

            julia $(joinpath(PKG_DIR, "build.jl"))

        (or `include("$(joinpath(PKG_DIR, "build.jl"))")` from a running Julia session).
        """)
    end
    return nothing
end

# =====================================================================
# ====== METHODS in crgMgr.c ======
# =====================================================================

"""
    crgDataSetRelease(dataSetId::Integer) -> Cint

Destroy the data of the given data set. Returns 1 if successful, 0 if failed.
"""
function crgDataSetRelease(dataSetId::Integer)
    return ccall((:crgDataSetRelease, LIBOPENCRG_PATH), Cint, (Cint,), dataSetId)
end

"""
    crgDataPrintHeader(dataSetId::Integer) -> Nothing

Print information contained in the CRG file's header (to stdout, via the
library's own message system).
"""
function crgDataPrintHeader(dataSetId::Integer)
    ccall((:crgDataPrintHeader, LIBOPENCRG_PATH), Cvoid, (Cint,), dataSetId)
    return nothing
end

"""
    crgDataPrintChannelInfo(dataSetId::Integer) -> Nothing

Print information about the CRG file's channels.
"""
function crgDataPrintChannelInfo(dataSetId::Integer)
    ccall((:crgDataPrintChannelInfo, LIBOPENCRG_PATH), Cvoid, (Cint,), dataSetId)
    return nothing
end

"""
    crgDataPrintRoadInfo(dataSetId::Integer) -> Nothing

Print information about the CRG road.
"""
function crgDataPrintRoadInfo(dataSetId::Integer)
    ccall((:crgDataPrintRoadInfo, LIBOPENCRG_PATH), Cvoid, (Cint,), dataSetId)
    return nothing
end

"""
    crgDataSetGetURange(dataSetId::Integer) -> (status, uMin, uMax)

Get the u co-ordinate range of a CRG data set. `status` is 1 upon success,
otherwise 0.
"""
function crgDataSetGetURange(dataSetId::Integer)
    uMin = Ref{Cdouble}(NaN)
    uMax = Ref{Cdouble}(NaN)
    status = ccall((:crgDataSetGetURange, LIBOPENCRG_PATH), Cint,
                    (Cint, Ref{Cdouble}, Ref{Cdouble}), dataSetId, uMin, uMax)
    return (status = status, uMin = uMin[], uMax = uMax[])
end

"""
    crgDataSetGetVRange(dataSetId::Integer) -> (status, vMin, vMax)

Get the v co-ordinate range of a CRG data set. `status` is 1 upon success,
otherwise 0.
"""
function crgDataSetGetVRange(dataSetId::Integer)
    vMin = Ref{Cdouble}(NaN)
    vMax = Ref{Cdouble}(NaN)
    status = ccall((:crgDataSetGetVRange, LIBOPENCRG_PATH), Cint,
                    (Cint, Ref{Cdouble}, Ref{Cdouble}), dataSetId, vMin, vMax)
    return (status = status, vMin = vMin[], vMax = vMax[])
end

"""
    crgDataSetGetIncrements(dataSetId::Integer) -> (status, uInc, vInc)

Get the u and v co-ordinate increments. `vInc` is 0 if explicit v sections
have been defined. `status` is 1 upon success, otherwise 0.
"""
function crgDataSetGetIncrements(dataSetId::Integer)
    uInc = Ref{Cdouble}(NaN)
    vInc = Ref{Cdouble}(NaN)
    status = ccall((:crgDataSetGetIncrements, LIBOPENCRG_PATH), Cint,
                    (Cint, Ref{Cdouble}, Ref{Cdouble}), dataSetId, uInc, vInc)
    return (status = status, uInc = uInc[], vInc = vInc[])
end

"""
    crgDataSetGetUtilityDataClosedTrack(dataSetId::Integer) -> (status, uIsClosed, uCloseMin, uCloseMax)

Get closed track utility data. `uIsClosed` records whether the reference line
can be closed; `uCloseMin`/`uCloseMax` are NaN if `uIsClosed == 0`. `status`
is 1 upon success, otherwise 0.
"""
function crgDataSetGetUtilityDataClosedTrack(dataSetId::Integer)
    uIsClosed = Ref{Cint}(0)
    uCloseMin = Ref{Cdouble}(NaN)
    uCloseMax = Ref{Cdouble}(NaN)
    status = ccall((:crgDataSetGetUtilityDataClosedTrack, LIBOPENCRG_PATH), Cint,
                    (Cint, Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}),
                    dataSetId, uIsClosed, uCloseMin, uCloseMax)
    return (status = status, uIsClosed = uIsClosed[], uCloseMin = uCloseMin[], uCloseMax = uCloseMax[])
end

"""
    crgMemRelease() -> Nothing

Release all data held by the crg library.
"""
function crgMemRelease()
    ccall((:crgMemRelease, LIBOPENCRG_PATH), Cvoid, ())
    return nothing
end

"""
    crgGetReleaseInfo() -> String

Get the release string indicating the current version of the underlying C
library.
"""
function crgGetReleaseInfo()
    ptr = ccall((:crgGetReleaseInfo, LIBOPENCRG_PATH), Cstring, ())
    return unsafe_string(ptr)
end

# =====================================================================
# ====== METHODS in crgMsg.c ======
# =====================================================================

"""
    crgMsgSetLevel(level::Integer) -> Nothing

Set the maximum level of messages that will be handled/printed by the
library (e.g. `dCrgMsgLevelNone = 0` .. `dCrgMsgLevelDebug = 5`, per
`crgBaseLib.h`).
"""
function crgMsgSetLevel(level::Integer)
    ccall((:crgMsgSetLevel, LIBOPENCRG_PATH), Cvoid, (Cint,), level)
    return nothing
end

# =====================================================================
# ====== METHODS in crgLoader.c ======
# =====================================================================

"""
    crgCheck(dataSetId::Integer) -> Cint

Check CRG data for consistency and accuracy. Returns truthy (nonzero) if
`crgData` is valid.
"""
function crgCheck(dataSetId::Integer)
    return ccall((:crgCheck, LIBOPENCRG_PATH), Cint, (Cint,), dataSetId)
end

"""
    crgLoaderReadFile(filename::AbstractString) -> Cint

Load CRG data from an existing IPL-formatted file. Returns the identifier of
the resulting data set, or 0 if not successful.
"""
function crgLoaderReadFile(filename::AbstractString)
    return ccall((:crgLoaderReadFile, LIBOPENCRG_PATH), Cint, (Cstring,), filename)
end

"""
    crgLoaderSuppressFileNotFoundFatalMsg(suppress::Bool) -> Nothing

Suppress the fatal message if the file was not found in `crgLoaderAddFile`
and `crgLoaderReadFile`.
"""
function crgLoaderSuppressFileNotFoundFatalMsg(suppress::Bool)
    ccall((:crgLoaderSuppressFileNotFoundFatalMsg, LIBOPENCRG_PATH), Cvoid, (Bool,), suppress)
    return nothing
end

# =====================================================================
# ====== METHODS in crgContactPoint.c ======
# =====================================================================

"""
    crgContactPointCreate(dataSetId::Integer) -> Cint

Create a new contact point working on the indicated data set. Returns the id
of the new contact point, or -1 if an error occurred.
"""
function crgContactPointCreate(dataSetId::Integer)
    return ccall((:crgContactPointCreate, LIBOPENCRG_PATH), Cint, (Cint,), dataSetId)
end

"""
    crgContactPointDelete(cpId::Integer) -> Cint

Delete a contact point and its associated data (not: crgData). Returns 1 if
successful, otherwise 0.
"""
function crgContactPointDelete(cpId::Integer)
    return ccall((:crgContactPointDelete, LIBOPENCRG_PATH), Cint, (Cint,), cpId)
end

"""
    crgContactPointDeleteAll(dataSetId::Integer) -> Nothing

Delete all contact points and associated data for a given data set (not:
crgData). If `dataSetId == -1`, all contact points of all data sets are
deleted.
"""
function crgContactPointDeleteAll(dataSetId::Integer)
    ccall((:crgContactPointDeleteAll, LIBOPENCRG_PATH), Cvoid, (Cint,), dataSetId)
    return nothing
end

# =====================================================================
# ====== METHODS in crgEvalxy2uv.c / crgEvaluv2xy.c ======
# =====================================================================

"""
    crgEvalxy2uv(cpId::Integer, x::Real, y::Real) -> (status, u, v)

Convert a given (x,y) position into the corresponding (u,v) position.
`status` is 1 if successful, otherwise 0.
"""
function crgEvalxy2uv(cpId::Integer, x::Real, y::Real)
    u = Ref{Cdouble}(NaN)
    v = Ref{Cdouble}(NaN)
    status = ccall((:crgEvalxy2uv, LIBOPENCRG_PATH), Cint,
                    (Cint, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
                    cpId, x, y, u, v)
    return (status = status, u = u[], v = v[])
end

"""
    crgEvaluv2xy(cpId::Integer, u::Real, v::Real) -> (status, x, y)

Convert a given (u,v) position into the corresponding (x,y) position.
`status` is 1 if successful, otherwise 0.
"""
function crgEvaluv2xy(cpId::Integer, u::Real, v::Real)
    x = Ref{Cdouble}(NaN)
    y = Ref{Cdouble}(NaN)
    status = ccall((:crgEvaluv2xy, LIBOPENCRG_PATH), Cint,
                    (Cint, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
                    cpId, u, v, x, y)
    return (status = status, x = x[], y = y[])
end

# =====================================================================
# ====== METHODS in crgEvalz.c ======
# =====================================================================

"""
    crgEvaluv2z(cpId::Integer, u::Real, v::Real) -> (status, z)

Compute the z value at a given (u,v) position using bilinear interpolation.
`status` is 1 if successful, otherwise 0.
"""
function crgEvaluv2z(cpId::Integer, u::Real, v::Real)
    z = Ref{Cdouble}(NaN)
    status = ccall((:crgEvaluv2z, LIBOPENCRG_PATH), Cint,
                    (Cint, Cdouble, Cdouble, Ref{Cdouble}),
                    cpId, u, v, z)
    return (status = status, z = z[])
end

"""
    crgEvalxy2z(cpId::Integer, x::Real, y::Real) -> (status, z)

Compute the z value at a given (x,y) position using bilinear interpolation.
`status` is 1 if successful, otherwise 0.
"""
function crgEvalxy2z(cpId::Integer, x::Real, y::Real)
    z = Ref{Cdouble}(NaN)
    status = ccall((:crgEvalxy2z, LIBOPENCRG_PATH), Cint,
                    (Cint, Cdouble, Cdouble, Ref{Cdouble}),
                    cpId, x, y, z)
    return (status = status, z = z[])
end

# =====================================================================
# ====== METHODS in crgEvalpk.c ======
# =====================================================================

"""
    crgEvaluv2pk(cpId::Integer, u::Real, v::Real) -> (status, phi, curv)

Compute the heading and curvature value at a given (u,v) position and store
it in the contact point structure. `status` is 1 if successful, otherwise 0.
"""
function crgEvaluv2pk(cpId::Integer, u::Real, v::Real)
    phi = Ref{Cdouble}(NaN)
    curv = Ref{Cdouble}(NaN)
    status = ccall((:crgEvaluv2pk, LIBOPENCRG_PATH), Cint,
                    (Cint, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
                    cpId, u, v, phi, curv)
    return (status = status, phi = phi[], curv = curv[])
end

"""
    crgEvalxy2pk(cpId::Integer, x::Real, y::Real) -> (status, phi, curv)

Compute the heading and curvature value at a given (x,y) position and store
it in the contact point structure. `status` is 1 if successful, otherwise 0.
"""
function crgEvalxy2pk(cpId::Integer, x::Real, y::Real)
    phi = Ref{Cdouble}(NaN)
    curv = Ref{Cdouble}(NaN)
    status = ccall((:crgEvalxy2pk, LIBOPENCRG_PATH), Cint,
                    (Cint, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
                    cpId, x, y, phi, curv)
    return (status = status, phi = phi[], curv = curv[])
end

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

end # module LibOpenCRG
