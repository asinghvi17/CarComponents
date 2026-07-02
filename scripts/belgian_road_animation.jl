using CarComponents
using ModelingToolkit
using OrdinaryDiffEqDefault
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit
using ADTypes
using SciMLBase
using MultibodyComponents
using DyadCompilerPasses
using GLMakie

const REPO_ROOT = dirname(@__DIR__)
const OUTDIR = joinpath(REPO_ROOT, "assets")
const VIDEO_PATH = joinpath(OUTDIR, "belgian_road_straight_car_with_suspension_dashboard.mp4")
const PREVIEW_PATH = joinpath(OUTDIR, "belgian_road_straight_car_dashboard_preview.png")
const G0 = 9.80665

function compile_controlled_belgian_model(road_profile)
    @named model = CarComponents.ControlledBelgianRoadStraightLineCar(
        wheel_elastic_contact = true,
        steer_limit = 0.25,
        heading_gain = 0.8,
        lateral_gain = 0.25,
        road_profile = road_profile,
    )

    reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm(
        inline_linear_sccs = true,
        analytical_linear_scc_limit = 1,
    )

    sys = ModelingToolkit.mtkcompile(
        model;
        additional_passes = [],
        reassemble_alg,
        optimize = [DyadCompilerPasses.LDIV_RULE],
    )

    return model, sys
end

function initial_conditions(sys, road_profile; speed = 0.8, baseline = nothing)
    corners = [
        sys.excited_suspension_fr,
        sys.excited_suspension_fl,
        sys.excited_suspension_br,
        sys.excited_suspension_bl,
    ]
    wheel_dirs = [1.0, -1.0, 1.0, -1.0]

    defs = Pair[]
    for (corner, dir) in zip(corners, wheel_dirs)
        push!(defs, corner.wheel.wheeljoint.v_small => 1e-3)
        push!(defs, corner.suspension.ks => 5 * 44000)
        push!(defs, corner.suspension.cs => 5 * 4000)
        push!(defs, corner.suspension.r2.phi => 5.932380614359173)
        push!(defs, corner.wheel_rotation.phi => 0.0)
        push!(defs, corner.wheel_rotation.w => -dir * speed / 0.2)
    end

    body_y0 = baseline === nothing ? 0.193 + road_profile(0.0) : baseline.body_y

    append!(defs, [
        sys.back_front.body.r_0[1] => 0.0,
        sys.back_front.body.r_0[2] => body_y0,
        sys.back_front.body.r_0[3] => 0.0,
        sys.back_front.body.v_0[1] => speed,
        sys.back_front.body.v_0[2] => 0.0,
        sys.back_front.body.v_0[3] => 0.0,
        sys.back_front.body.phi[2] => 0.0,
        sys.excited_suspension_fr.steer_rotation.phi => 0.0,
        sys.excited_suspension_fl.steer_rotation.phi => 0.0,
        sys.excited_suspension_fr.steering_position.w => 0.0,
        sys.excited_suspension_fl.steering_position.w => 0.0,
        sys.world.render => true,
        sys.world.axis_radius => 0.01,
    ])

    if baseline !== nothing
        append!(defs, [
            sys.excited_suspension_fr.suspension.springdamper.s => baseline.fr_s,
            sys.excited_suspension_fl.suspension.springdamper.s => baseline.fl_s,
            sys.excited_suspension_br.suspension.springdamper.s => baseline.br_s,
            sys.excited_suspension_bl.suspension.springdamper.s => baseline.bl_s,
        ])
    end

    return defs
end

function simulate_belgian_road(sys, road_profile; speed = 0.8, tstop = 12.0, save_dt = 1 / 60, baseline = nothing)
    prob = ODEProblem(
        sys,
        initial_conditions(sys, road_profile; speed, baseline),
        (0.0, tstop);
        optimize = :basic,
        saveat = save_dt,
    )

    sol = solve(
        prob,
        Rodas5P(autodiff = AutoFiniteDiff());
        initializealg = BrownFullBasicInit(),
        abstol = 1e-6,
        reltol = 1e-6,
        dtmax = 0.01,
    )

    if !SciMLBase.successful_retcode(sol)
        error("Belgian road simulation failed: retcode=$(sol.retcode), t_end=$(sol.t[end])")
    end

    return sol
end

function static_baseline(sys, road_profile; settle_time = 4.0)
    sol = simulate_belgian_road(sys, road_profile; speed = 0.0, tstop = settle_time, save_dt = 1 / 60)
    t = sol.t[end]

    fr_x = sol(t; idxs = sys.wheel_position_fr)
    fl_x = sol(t; idxs = sys.wheel_position_fl)
    br_x = sol(t; idxs = sys.wheel_position_br)
    bl_x = sol(t; idxs = sys.wheel_position_bl)
    body_x = sol(t; idxs = sys.back_front.body.r_0[1])

    return (
        body_y = sol(t; idxs = sys.back_front.body.r_0[2]),
        body_ay = sol(t; idxs = sys.back_front.body.a_0[2]),
        body_x = body_x,
        road_body_center = road_profile(body_x),
        road_wheel_average = (road_profile(fr_x) + road_profile(fl_x) + road_profile(br_x) + road_profile(bl_x)) / 4,
        fr_s = sol(t; idxs = sys.excited_suspension_fr.suspension.springdamper.s),
        fl_s = sol(t; idxs = sys.excited_suspension_fl.suspension.springdamper.s),
        br_s = sol(t; idxs = sys.excited_suspension_br.suspension.springdamper.s),
        bl_s = sol(t; idxs = sys.excited_suspension_bl.suspension.springdamper.s),
    )
end

function road_mesh(road_profile; xspan = (-0.5, 10.5), zspan = (-1.2, 1.2), nx = 180, nz = 18)
    xs = range(xspan[1], xspan[2]; length = nx)
    zs = range(zspan[1], zspan[2]; length = nz)

    points = GLMakie.Point3f[]
    for z in zs, x in xs
        push!(points, GLMakie.Point3f(Float32(x), Float32(road_profile(x)), Float32(z)))
    end

    faces = GLMakie.GeometryBasics.GLTriangleFace[]
    for j in 1:(nz - 1), i in 1:(nx - 1)
        a = (j - 1) * nx + i
        b = a + 1
        c = j * nx + i
        d = c + 1
        push!(faces, GLMakie.GeometryBasics.GLTriangleFace(a, c, b))
        push!(faces, GLMakie.GeometryBasics.GLTriangleFace(b, c, d))
    end

    return GLMakie.GeometryBasics.Mesh(points, faces)
end

function dashboard_data(sys, sol, road_profile, baseline)
    body_x = sol[sys.back_front.body.r_0[1]]
    fr_x = sol[sys.wheel_position_fr]
    fl_x = sol[sys.wheel_position_fl]
    br_x = sol[sys.wheel_position_br]
    bl_x = sol[sys.wheel_position_bl]

    road_body_center_mm = [1000 * (road_profile(x) - baseline.road_body_center) for x in body_x]
    road_wheel_average_mm = [
        1000 * ((road_profile(fr_x[i]) + road_profile(fl_x[i]) + road_profile(br_x[i]) + road_profile(bl_x[i])) / 4 - baseline.road_wheel_average)
        for i in eachindex(body_x)
    ]

    return (
        body_x = body_x,
        body_z = sol[sys.back_front.body.r_0[3]],
        body_ay_g = (sol[sys.back_front.body.a_0[2]] .- baseline.body_ay) ./ G0,
        heading = -sol[sys.back_front.body.phi[2]],
        steer = sol[sys.controller.steer_angle],
        road_body_center_mm = road_body_center_mm,
        road_wheel_average_mm = road_wheel_average_mm,
        fr_s = 1000 .* (sol[sys.excited_suspension_fr.suspension.springdamper.s] .- baseline.fr_s),
        fl_s = 1000 .* (sol[sys.excited_suspension_fl.suspension.springdamper.s] .- baseline.fl_s),
        br_s = 1000 .* (sol[sys.excited_suspension_br.suspension.springdamper.s] .- baseline.br_s),
        bl_s = 1000 .* (sol[sys.excited_suspension_bl.suspension.springdamper.s] .- baseline.bl_s),
    )
end

function padded_limits(values; pad_fraction = 0.08, min_span = 1e-6)
    lo, hi = extrema(values)
    span = max(hi - lo, min_span)
    return lo - pad_fraction * span, hi + pad_fraction * span
end

function build_dashboard(model, sys, sol, road_profile, baseline; preview_time = 0.0)
    data = dashboard_data(sys, sol, road_profile, baseline)

    fig, time_obs, scene = MultibodyComponents.render(
        model,
        sol,
        preview_time;
        slider = false,
        x = 5.2,
        y = 3.2,
        z = 6.0,
        lookat = GLMakie.Vec3f(5.0, 0.0, 0.0),
        up = GLMakie.Vec3f(0.0, 1.0, 0.0),
        show_axis = true,
        traces = [sys.back_front.frame_cm],
        size = (1920, 1080),
    )

    side = fig[1, 2] = GLMakie.GridLayout()
    GLMakie.colsize!(fig.layout, 1, GLMakie.Relative(0.62))
    GLMakie.colsize!(fig.layout, 2, GLMakie.Relative(0.38))

    GLMakie.Label(
        fig[0, :],
        "Controlled straight-line car on Belgian-block road: physical render + suspension plots";
        fontsize = 26,
        tellwidth = false,
    )

    GLMakie.mesh!(scene, road_mesh(road_profile); color = (:gray65, 0.55), transparency = true)

    road_x = range(0.0, 10.2; length = 300)
    GLMakie.lines!(scene, road_x, road_profile.(road_x) .+ 0.025, zeros(length(road_x)); color = :deepskyblue, linewidth = 5)
    GLMakie.lines!(scene, data.body_x, road_profile.(data.body_x) .+ 0.04, data.body_z; color = :orange, linewidth = 4)

    status = GLMakie.Observable("t = 0.00 s    x = 0.00 m    ay = 0.00 g")
    GLMakie.Label(side[1, 1], status; fontsize = 18, tellwidth = false, halign = :left)

    ax1 = GLMakie.Axis(
        side[2, 1];
        xlabel = "distance along road x [m]",
        ylabel = "relative road height [mm]",
        title = "Road input relative to standing still",
    )
    GLMakie.lines!(ax1, data.body_x, data.road_wheel_average_mm; color = :black, linewidth = 2, label = "average road under wheels")
    GLMakie.lines!(ax1, data.body_x, data.road_body_center_mm; color = :gray45, linestyle = :dash, linewidth = 2, label = "road at chassis center")
    GLMakie.axislegend(ax1; position = :lb, labelsize = 12)
    GLMakie.xlims!(ax1, 0, 10.2)
    GLMakie.ylims!(ax1, padded_limits(vcat(data.road_wheel_average_mm, data.road_body_center_mm))...)

    ax2 = GLMakie.Axis(
        side[3, 1];
        xlabel = "distance along road x [m]",
        ylabel = "vertical acceleration [g]",
        title = "Chassis center vertical acceleration",
    )
    GLMakie.hlines!(ax2, [0.0]; color = :gray70, linewidth = 1)
    GLMakie.lines!(ax2, data.body_x, data.body_ay_g; color = :royalblue, linewidth = 2, label = "chassis ay")
    GLMakie.axislegend(ax2; position = :lb, labelsize = 12)
    GLMakie.xlims!(ax2, 0, 10.2)
    GLMakie.ylims!(ax2, padded_limits(data.body_ay_g)...)

    ax3 = GLMakie.Axis(
        side[4, 1];
        xlabel = "time [s]",
        ylabel = "spring travel relative to static [mm]",
        title = "Suspension travel check",
    )
    GLMakie.lines!(ax3, sol.t, data.fr_s; color = :red, linewidth = 2, label = "FR")
    GLMakie.lines!(ax3, sol.t, data.fl_s; color = :magenta, linewidth = 2, label = "FL")
    GLMakie.lines!(ax3, sol.t, data.br_s; color = :green, linewidth = 2, label = "BR")
    GLMakie.lines!(ax3, sol.t, data.bl_s; color = :cyan4, linewidth = 2, label = "BL")
    GLMakie.axislegend(ax3; position = :lb, nbanks = 2, labelsize = 12)
    GLMakie.xlims!(ax3, sol.t[1], sol.t[end])

    ax4 = GLMakie.Axis(
        side[5, 1];
        xlabel = "time [s]",
        ylabel = "z [m] / angles [rad]",
        title = "Straight-line tracking and steering effort",
    )
    GLMakie.lines!(ax4, sol.t, data.body_z; color = :orange, linewidth = 2, label = "lateral z")
    GLMakie.lines!(ax4, sol.t, data.heading; color = :purple, linewidth = 2, label = "heading")
    GLMakie.lines!(ax4, sol.t, data.steer; color = :brown, linewidth = 2, label = "steer")
    GLMakie.axislegend(ax4; position = :lb, nbanks = 3, labelsize = 12)
    GLMakie.xlims!(ax4, sol.t[1], sol.t[end])

    x_marker = GLMakie.Observable([data.body_x[1]])
    t_marker = GLMakie.Observable([sol.t[1]])
    road_marker = GLMakie.Observable([data.road_wheel_average_mm[1]])
    accel_marker = GLMakie.Observable([data.body_ay_g[1]])

    GLMakie.vlines!(ax1, x_marker; color = :gray30, linewidth = 2)
    GLMakie.scatter!(ax1, x_marker, road_marker; color = :black, markersize = 12)
    GLMakie.vlines!(ax2, x_marker; color = :gray30, linewidth = 2)
    GLMakie.scatter!(ax2, x_marker, accel_marker; color = :royalblue, markersize = 12)
    GLMakie.vlines!(ax3, t_marker; color = :gray30, linewidth = 2)
    GLMakie.vlines!(ax4, t_marker; color = :gray30, linewidth = 2)

    function update_markers!(time)
        x = sol(time; idxs = sys.back_front.body.r_0[1])
        fr_x = sol(time; idxs = sys.wheel_position_fr)
        fl_x = sol(time; idxs = sys.wheel_position_fl)
        br_x = sol(time; idxs = sys.wheel_position_br)
        bl_x = sol(time; idxs = sys.wheel_position_bl)
        road_average_mm = 1000 * ((road_profile(fr_x) + road_profile(fl_x) + road_profile(br_x) + road_profile(bl_x)) / 4 - baseline.road_wheel_average)
        ay_g = (sol(time; idxs = sys.back_front.body.a_0[2]) - baseline.body_ay) / G0

        x_marker[] = [x]
        t_marker[] = [time]
        road_marker[] = [road_average_mm]
        accel_marker[] = [ay_g]
        status[] = "t = $(round(time, digits = 2)) s    x = $(round(x, digits = 2)) m    ay = $(round(ay_g, digits = 2)) g"
    end

    update_markers!(preview_time)
    return fig, time_obs, update_markers!
end

function record_dashboard(fig, time_obs, update_markers!, sol; output = VIDEO_PATH, framerate = 30)
    times = collect(range(sol.t[1], sol.t[end]; step = 1 / framerate))
    times[end] < sol.t[end] && push!(times, sol.t[end])

    return GLMakie.record(fig, output, times; framerate) do time
        time_obs[] = time
        update_markers!(time)
    end
end

function main(; output = VIDEO_PATH, speed = 0.8, tstop = 12.0, static_settle_time = 4.0, framerate = 30)
    mkpath(OUTDIR)
    GLMakie.activate!()

    road_profile = CarComponents.belgian_block_centerline_profile_interpolator()

    println("Compiling controlled Belgian-road car...")
    model, sys = compile_controlled_belgian_model(road_profile)
    println("Compiled: $(length(unknowns(sys))) unknowns, $(length(equations(sys))) equations")

    println("Finding standing-still baseline...")
    baseline = static_baseline(sys, road_profile; settle_time = static_settle_time)

    println("Simulating moving run...")
    sol = simulate_belgian_road(sys, road_profile; speed, tstop, save_dt = 1 / 60, baseline)
    println("Simulation retcode: $(sol.retcode), samples: $(length(sol.t))")

    println("Building dashboard...")
    fig, time_obs, update_markers! = build_dashboard(model, sys, sol, road_profile, baseline)
    GLMakie.save(PREVIEW_PATH, fig)
    println("Preview saved to $PREVIEW_PATH")

    println("Recording video to $output ...")
    filename = record_dashboard(fig, time_obs, update_markers!, sol; output, framerate)
    println("Video saved to $filename")

    return filename
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
