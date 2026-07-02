using Test
using LibOpenCRG

const DATA_DIR = joinpath(@__DIR__, "data")

# Keep the underlying C library quiet during tests: only warning/fatal
# messages (dCrgMsgLevelWarn = 2, see csrc/inc/crgBaseLib.h) will be
# printed; at its default level the library dumps a verbose per-file
# NOTICE-level statistics report on every successful load.
crgMsgSetLevel(2)

@testset "LibOpenCRG" begin
    @testset "crgGetReleaseInfo" begin
        info = crgGetReleaseInfo()
        @test info isa String
        @test !isempty(info)
        @test occursin("OpenCRG", info)
    end

    @testset "smoke test: $(file)" for file in (
        "handmade_curved_minimalist.crg",
        "handmade_curved_banked_sloped.crg",
    )
        path = joinpath(DATA_DIR, file)
        @test isfile(path)

        dataSetId = crgLoaderReadFile(path)
        @test dataSetId > 0

        try
            # crgCheck must report success (nonzero / no errors) on a valid
            # vendored example file.
            checkStatus = crgCheck(dataSetId)
            @test checkStatus != 0

            uRangeRes = crgDataSetGetURange(dataSetId)
            @test uRangeRes.status != 0
            @test uRangeRes.uMax > uRangeRes.uMin

            vRangeRes = crgDataSetGetVRange(dataSetId)
            @test vRangeRes.status != 0
            @test vRangeRes.vMax >= vRangeRes.vMin

            incRes = crgDataSetGetIncrements(dataSetId)
            @test incRes.status != 0
            @test incRes.uInc > 0

            closedRes = crgDataSetGetUtilityDataClosedTrack(dataSetId)
            @test closedRes.status != 0

            cpId = crgContactPointCreate(dataSetId)
            @test cpId >= 0

            try
                # Evaluate elevation at three interior (u, v) points spread
                # across the data set's valid range.
                uSpan = uRangeRes.uMax - uRangeRes.uMin
                vSpan = vRangeRes.vMax - vRangeRes.vMin
                testPoints = [
                    (uRangeRes.uMin + f * uSpan, vRangeRes.vMin + f * vSpan)
                    for f in (0.25, 0.5, 0.75)
                ]

                for (u, v) in testPoints
                    res = crgEvaluv2z(cpId, u, v)
                    @test res.status != 0
                    @test res.z isa Real
                    @test isfinite(res.z)
                end

                # Bonus coverage: heading/curvature and the xy<->uv coordinate
                # maps should also produce finite results on the same contact
                # point, and a uv->xy->uv round trip should be self-consistent.
                uMid, vMid = testPoints[2]

                pk = crgEvaluv2pk(cpId, uMid, vMid)
                @test pk.status != 0
                @test isfinite(pk.phi)
                @test isfinite(pk.curv)

                xy = crgEvaluv2xy(cpId, uMid, vMid)
                @test xy.status != 0
                @test isfinite(xy.x)
                @test isfinite(xy.y)

                uv = crgEvalxy2uv(cpId, xy.x, xy.y)
                @test uv.status != 0
                @test isapprox(uv.u, uMid; atol = 1e-6)
                @test isapprox(uv.v, vMid; atol = 1e-6)
            finally
                @test crgContactPointDelete(cpId) != 0
            end
        finally
            @test crgDataSetRelease(dataSetId) != 0
        end
    end

    @testset "modifier bindings" begin
        dsId = crgLoaderReadFile(joinpath(@__DIR__, "data", "handmade_curved_minimalist.crg"))
        @test dsId != 0

        # crgDataSetCreate seeds the modifier list with DEFAULT
        # dCrgModRefPointX/Y/Z/Phi = 0.0 entries (crgOptionSetDefaultModifiers,
        # in the vendored crgOptionMgmt.c). crgDataApplyTransformations checks
        # for *any* dCrgModRefPoint* entry before it ever looks at
        # dCrgModRefLineOffset*, so those defaults being merely "set" (even
        # though their value is 0.0) silently take over the whole transform
        # and turn it into a same-point/zero-rotation no-op -- verified by
        # running this exact sequence without the removeAll call first: the
        # offset modifiers are accepted (status 1) but never applied, and
        # r.x comes back 0.0, not 100.0. crgDataSetModifierRemoveAll clears
        # those defaults first so the reference-line-offset path actually
        # runs.
        @test crgDataSetModifierRemoveAll(dsId) != 0
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
end

crgMemRelease()
