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
const VIDEO_PATH = joinpath(OUTDIR, "country_road_powertrain_car_dashboard.mp4")
const PREVIEW_PATH = joinpath(OUTDIR, "country_road_powertrain_car_dashboard_preview.png")
const G0 = 9.80665

const COUNTRY_ROAD_TARGET_AXLE_SPACING = 2.3481016274718134
const COUNTRY_ROAD_TARGET_WHEEL_RADIUS = 0.30760131319880757

function compile_controlled_country_model(road_surface, center_z_profile, center_heading_profile;
        target_speed = 6.0, drive_speed_gain = 450.0, drive_integral_gain = 25.0,
        brake_speed_gain = 700.0, max_drive_torque = 900.0, max_brake_torque = 1200.0)
    @named model = CarComponents.ControlledCountryRoadCarWithPowertrain(
        wheel_elastic_contact = true,
        axle_spacing = COUNTRY_ROAD_TARGET_AXLE_SPACING,
        wheel_radius = COUNTRY_ROAD_TARGET_WHEEL_RADIUS,
        steer_limit = 0.45,
        heading_gain = 3.0,
        lateral_gain = 1.5,
        target_speed = target_speed,
        drive_speed_gain = drive_speed_gain,
        drive_integral_gain = drive_integral_gain,
        brake_speed_gain = brake_speed_gain,
        max_drive_torque = max_drive_torque,
        max_brake_torque = max_brake_torque,
        road_surface = road_surface,
        center_z_profile = center_z_profile,
        center_heading_profile = center_heading_profile,
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

function country_initial_conditions(sys, road_surface, center_z_profile, center_heading_profile;
        speed = 6.0, baseline = nothing, start_x = 2.0)
    corners = [
        sys.excited_suspension_fr,
        sys.excited_suspension_fl,
        sys.excited_suspension_br,
        sys.excited_suspension_bl,
    ]
    wheel_dirs = [1.0, -1.0, 1.0, -1.0]

    start_z = center_z_profile(start_x)
    start_heading = center_heading_profile(start_x)
    body_y0 = baseline === nothing ? 0.193 + road_surface(start_x, start_z) : baseline.body_y

    defs = Pair[]
    for (corner, dir) in zip(corners, wheel_dirs)
        push!(defs, corner.wheel.wheeljoint.v_small => 1e-3)
        push!(defs, corner.suspension.ks => 5 * 44000)
        push!(defs, corner.suspension.cs => 5 * 4000)
        push!(defs, corner.suspension.r2.phi => 5.932380614359173)
        push!(defs, corner.wheel_rotation.phi => 0.0)
        push!(defs, corner.wheel_rotation.w => -dir * speed / COUNTRY_ROAD_TARGET_WHEEL_RADIUS)
    end

    append!(defs, [
        sys.back_front.body.r_0[1] => start_x,
        sys.back_front.body.r_0[2] => body_y0,
        sys.back_front.body.r_0[3] => start_z,
        sys.back_front.body.v_0[1] => speed * cos(start_heading),
        sys.back_front.body.v_0[2] => 0.0,
        sys.back_front.body.v_0[3] => speed * sin(start_heading),
        sys.back_front.body.phi[2] => -start_heading,
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

function simulate_country_road(sys, road_surface, center_z_profile, center_heading_profile;
        speed = 6.0, target_speed = 6.0, tstop = 50.0, save_dt = 1 / 30, baseline = nothing, start_x = 2.0)
    defs = country_initial_conditions(sys, road_surface, center_z_profile, center_heading_profile;
        speed = speed, baseline = baseline, start_x = start_x)
    append!(defs, [
        sys.target_speed => target_speed,
        sys.controller.target_speed => target_speed,
    ])

    prob = ODEProblem(
        sys,
        defs,
        (0.0, tstop);
        optimize = :basic,
        saveat = save_dt,
        warn_initialize_determined = false,
    )

    sol = solve(
        prob,
        Rodas5P(autodiff = AutoFiniteDiff());
        initializealg = BrownFullBasicInit(),
        abstol = 1e-6,
        reltol = 1e-6,
        dtmax = 0.02,
    )

    if !SciMLBase.successful_retcode(sol)
        error("Country-road simulation failed: retcode=$(sol.retcode), t_end=$(sol.t[end])")
    end

    return sol
end

function static_baseline(sys, road_surface, center_z_profile, center_heading_profile;
        settle_time = 4.0, start_x = 2.0)
    sol = simulate_country_road(
        sys,
        road_surface,
        center_z_profile,
        center_heading_profile;
        speed = 0.0,
        target_speed = 0.0,
        tstop = settle_time,
        save_dt = 1 / 60,
        start_x = start_x,
    )
    t = sol.t[end]
    body_x = sol(t; idxs = sys.back_front.body.r_0[1])
    body_z = sol(t; idxs = sys.back_front.body.r_0[3])

    return (
        body_y = sol(t; idxs = sys.back_front.body.r_0[2]),
        body_ay = sol(t; idxs = sys.back_front.body.a_0[2]),
        body_x = body_x,
        body_z = body_z,
        road_body_center = road_surface(body_x, body_z),
        fr_s = sol(t; idxs = sys.excited_suspension_fr.suspension.springdamper.s),
        fl_s = sol(t; idxs = sys.excited_suspension_fl.suspension.springdamper.s),
        br_s = sol(t; idxs = sys.excited_suspension_br.suspension.springdamper.s),
        bl_s = sol(t; idxs = sys.excited_suspension_bl.suspension.springdamper.s),
    )
end

function simulated_road_ribbon_mesh(data, road_surface, center_z_profile, center_heading_profile)
    points = GLMakie.Point3f[]
    for j in eachindex(data.v_axis), i in eachindex(data.x_axis)
        x = data.x_axis[i]
        heading = center_heading_profile(x)
        z = center_z_profile(x) + data.v_axis[j] / cos(heading)
        y = road_surface(x, z)
        push!(points, GLMakie.Point3f(Float32(x), Float32(y), Float32(z)))
    end

    nx = length(data.x_axis)
    nz = length(data.v_axis)
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

function simulated_centerline(data, road_surface, center_z_profile)
    x = collect(data.x_axis)
    z = [center_z_profile(xi) for xi in x]
    y = [road_surface(x[i], z[i]) for i in eachindex(x)]
    return x, y, z
end

function dashboard_data(sys, sol, road_data, road_surface, center_z_profile, center_heading_profile, baseline, target_speed)
    body_x = sol[sys.back_front.body.r_0[1]]
    body_z = sol[sys.back_front.body.r_0[3]]
    heading = -sol[sys.back_front.body.phi[2]]
    steer = sol[sys.controller.steer_angle]

    body_vx = sol[sys.back_front.body.v_0[1]]
    body_vz = sol[sys.back_front.body.v_0[3]]
    center_z = [center_z_profile(x) for x in body_x]
    center_heading = [center_heading_profile(x) for x in body_x]
    lateral_error = [cos(center_heading[i]) * (body_z[i] - center_z[i]) for i in eachindex(body_x)]
    heading_error = [atan(sin(center_heading[i] - heading[i]), cos(center_heading[i] - heading[i])) for i in eachindex(body_x)]
    road_speed = [cos(center_heading[i]) * body_vx[i] + sin(center_heading[i]) * body_vz[i] for i in eachindex(body_x)]
    target_speed_series = fill(target_speed, length(body_x))
    drive_torque = sol[sys.controller.drive_torque]
    brake_torque = sol[sys.controller.brake_torque]

    fr_x = sol[sys.wheel_position_fr]
    fl_x = sol[sys.wheel_position_fl]
    br_x = sol[sys.wheel_position_br]
    bl_x = sol[sys.wheel_position_bl]
    fr_z = sol[sys.wheel_lateral_position_fr]
    fl_z = sol[sys.wheel_lateral_position_fl]
    br_z = sol[sys.wheel_lateral_position_br]
    bl_z = sol[sys.wheel_lateral_position_bl]
    road_wheel_average_mm = [
        1000 * ((road_surface(fr_x[i], fr_z[i]) + road_surface(fl_x[i], fl_z[i]) + road_surface(br_x[i], br_z[i]) + road_surface(bl_x[i], bl_z[i])) / 4 - baseline.road_body_center)
        for i in eachindex(body_x)
    ]

    return (
        body_x = body_x,
        body_z = body_z,
        heading = heading,
        steer = steer,
        center_z = center_z,
        center_heading = center_heading,
        lateral_error = lateral_error,
        heading_error = heading_error,
        road_speed = road_speed,
        target_speed = target_speed_series,
        drive_torque = drive_torque,
        brake_torque = brake_torque,
        body_ay_g = (sol[sys.back_front.body.a_0[2]] .- baseline.body_ay) ./ G0,
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

function camera_pose(x, z, heading, road_height)
    eye = GLMakie.Vec3f(
        Float32(x - 13 * cos(heading) - 5 * sin(heading)),
        Float32(road_height + 7.0),
        Float32(z - 13 * sin(heading) + 5 * cos(heading)),
    )
    lookat = GLMakie.Vec3f(
        Float32(x + 7 * cos(heading)),
        Float32(road_height + 0.4),
        Float32(z + 7 * sin(heading)),
    )
    return eye, lookat, GLMakie.Vec3f(0, 1, 0)
end

function update_camera!(scene, x, z, heading, road_height)
    eye, lookat, up = camera_pose(x, z, heading, road_height)
    GLMakie.update_cam!(scene.scene, GLMakie.cameracontrols(scene.scene), eye, lookat, up)
end

function build_dashboard(model, sys, sol, road_data, road_surface, center_z_profile, center_heading_profile, baseline, target_speed; preview_time = 0.0)
    data = dashboard_data(sys, sol, road_data, road_surface, center_z_profile, center_heading_profile, baseline, target_speed)
    x0 = sol(preview_time; idxs = sys.back_front.body.r_0[1])
    z0 = sol(preview_time; idxs = sys.back_front.body.r_0[3])
    heading0 = -sol(preview_time; idxs = sys.back_front.body.phi[2])
    road0 = road_surface(x0, z0)
    eye0, lookat0, up0 = camera_pose(x0, z0, heading0, road0)

    fig, time_obs, scene = MultibodyComponents.render(
        model,
        sol,
        preview_time;
        slider = false,
        x = eye0[1],
        y = eye0[2],
        z = eye0[3],
        lookat = lookat0,
        up = up0,
        show_axis = true,
        traces = [sys.back_front.frame_cm],
        size = (1920, 1080),
    )

    side = fig[1, 2] = GLMakie.GridLayout()
    GLMakie.colsize!(fig.layout, 1, GLMakie.Relative(0.62))
    GLMakie.colsize!(fig.layout, 2, GLMakie.Relative(0.38))

    GLMakie.Label(
        fig[0, :],
        "Powered full-car model on transformed OpenCRG country road";
        fontsize = 26,
        tellwidth = false,
    )

    ribbon = simulated_road_ribbon_mesh(road_data, road_surface, center_z_profile, center_heading_profile)
    centerline_x, centerline_y, centerline_z = simulated_centerline(road_data, road_surface, center_z_profile)
    GLMakie.mesh!(scene, ribbon; color = (:gray65, 0.55), transparency = true)
    GLMakie.lines!(scene, centerline_x, centerline_y .+ 0.035, centerline_z; color = :deepskyblue, linewidth = 5)
    GLMakie.lines!(scene, data.body_x, [road_surface(data.body_x[i], data.body_z[i]) + 0.05 for i in eachindex(data.body_x)], data.body_z; color = :orange, linewidth = 4)

    status = GLMakie.Observable("t = 0.00 s    x = 0.00 m    speed = 0.00 m/s    drive = 0.00 kN*m    brake = 0.00 kN*m")
    GLMakie.Label(side[1, 1], status; fontsize = 18, tellwidth = false, halign = :left)

    ax1 = GLMakie.Axis(
        side[2, 1];
        xlabel = "world x [m]",
        ylabel = "world z [m]",
        title = "Country-road route and vehicle path",
        aspect = GLMakie.DataAspect(),
    )
    GLMakie.lines!(ax1, centerline_x, centerline_z; color = :black, linewidth = 2, label = "simulated reference line")
    GLMakie.lines!(ax1, data.body_x, data.body_z; color = :orange, linewidth = 2, label = "chassis path")
    GLMakie.axislegend(ax1; position = :lb, labelsize = 12)

    ax2 = GLMakie.Axis(
        side[3, 1];
        xlabel = "distance x [m]",
        ylabel = "road input [mm] / ay [g]",
        title = "Road input and chassis vertical acceleration",
    )
    GLMakie.lines!(ax2, data.body_x, data.road_wheel_average_mm; color = :black, linewidth = 2, label = "average road under wheels [mm]")
    GLMakie.lines!(ax2, data.body_x, data.body_ay_g; color = :royalblue, linewidth = 2, label = "chassis ay [g]")
    GLMakie.axislegend(ax2; position = :lb, labelsize = 12)
    GLMakie.ylims!(ax2, padded_limits(vcat(data.road_wheel_average_mm, data.body_ay_g))...)

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
        ylabel = "speed [m/s] / torque [kN*m] / steer [rad]",
        title = "Speed control and wheel actuation",
    )
    GLMakie.lines!(ax4, sol.t, data.road_speed; color = :black, linewidth = 2, label = "road speed [m/s]")
    GLMakie.lines!(ax4, sol.t, data.target_speed; color = :gray45, linestyle = :dash, linewidth = 2, label = "target speed [m/s]")
    GLMakie.lines!(ax4, sol.t, data.drive_torque ./ 1000; color = :green, linewidth = 2, label = "drive torque [kN*m]")
    GLMakie.lines!(ax4, sol.t, data.brake_torque ./ 1000; color = :red, linewidth = 2, label = "brake torque/wheel [kN*m]")
    GLMakie.lines!(ax4, sol.t, data.steer; color = :brown, linewidth = 2, label = "steer [rad]")
    GLMakie.axislegend(ax4; position = :lb, nbanks = 2, labelsize = 12)
    GLMakie.xlims!(ax4, sol.t[1], sol.t[end])

    route_x_marker = GLMakie.Observable([data.body_x[1]])
    route_z_marker = GLMakie.Observable([data.body_z[1]])
    x_marker = GLMakie.Observable([data.body_x[1]])
    t_marker = GLMakie.Observable([sol.t[1]])
    road_marker = GLMakie.Observable([data.road_wheel_average_mm[1]])
    accel_marker = GLMakie.Observable([data.body_ay_g[1]])

    GLMakie.scatter!(ax1, route_x_marker, route_z_marker; color = :orange, markersize = 14)
    GLMakie.vlines!(ax2, x_marker; color = :gray30, linewidth = 2)
    GLMakie.scatter!(ax2, x_marker, road_marker; color = :black, markersize = 12)
    GLMakie.scatter!(ax2, x_marker, accel_marker; color = :royalblue, markersize = 12)
    GLMakie.vlines!(ax3, t_marker; color = :gray30, linewidth = 2)
    GLMakie.vlines!(ax4, t_marker; color = :gray30, linewidth = 2)

    function update_markers!(time)
        x = sol(time; idxs = sys.back_front.body.r_0[1])
        z = sol(time; idxs = sys.back_front.body.r_0[3])
        heading = -sol(time; idxs = sys.back_front.body.phi[2])
        road_height = road_surface(x, z)
        center_z = center_z_profile(x)
        center_heading = center_heading_profile(x)
        vx = sol(time; idxs = sys.back_front.body.v_0[1])
        vz = sol(time; idxs = sys.back_front.body.v_0[3])
        road_speed = cos(center_heading) * vx + sin(center_heading) * vz
        drive_torque = sol(time; idxs = sys.controller.drive_torque)
        brake_torque = sol(time; idxs = sys.controller.brake_torque)
        lateral = cos(center_heading) * (z - center_z)
        ay_g = (sol(time; idxs = sys.back_front.body.a_0[2]) - baseline.body_ay) / G0
        fr_x = sol(time; idxs = sys.wheel_position_fr)
        fl_x = sol(time; idxs = sys.wheel_position_fl)
        br_x = sol(time; idxs = sys.wheel_position_br)
        bl_x = sol(time; idxs = sys.wheel_position_bl)
        fr_z = sol(time; idxs = sys.wheel_lateral_position_fr)
        fl_z = sol(time; idxs = sys.wheel_lateral_position_fl)
        br_z = sol(time; idxs = sys.wheel_lateral_position_br)
        bl_z = sol(time; idxs = sys.wheel_lateral_position_bl)
        road_average_mm = 1000 * ((road_surface(fr_x, fr_z) + road_surface(fl_x, fl_z) + road_surface(br_x, br_z) + road_surface(bl_x, bl_z)) / 4 - baseline.road_body_center)

        route_x_marker[] = [x]
        route_z_marker[] = [z]
        x_marker[] = [x]
        t_marker[] = [time]
        road_marker[] = [road_average_mm]
        accel_marker[] = [ay_g]
        status[] = "t = $(round(time, digits = 2)) s    x = $(round(x, digits = 1)) m    speed = $(round(road_speed, digits = 2)) m/s    drive = $(round(drive_torque / 1000, digits = 2)) kN*m    brake = $(round(brake_torque / 1000, digits = 2)) kN*m    lateral = $(round(lateral, digits = 3)) m"
        update_camera!(scene, x, z, heading, road_height)
    end

    update_markers!(preview_time)
    return fig, time_obs, update_markers!
end

function record_dashboard(fig, time_obs, update_markers!, sol; output = VIDEO_PATH, framerate = 24)
    times = collect(range(sol.t[1], sol.t[end]; step = 1 / framerate))
    times[end] < sol.t[end] && push!(times, sol.t[end])

    return GLMakie.record(fig, output, times; framerate) do time
        time_obs[] = time
        update_markers!(time)
    end
end

function main(; output = VIDEO_PATH, speed = 0.0, target_speed = 6.0, tstop = 50.0, static_settle_time = 4.0,
        framerate = 24, road_spacing = 0.5, lateral_stride = 8, start_x = 2.0)
    mkpath(OUTDIR)
    GLMakie.activate!()

    println("Preparing transformed country-road data...")
    road_data = CarComponents.country_road_curved_data(longitudinal_spacing = road_spacing, lateral_stride = lateral_stride)
    road_surface = CarComponents.opencrg_curved_surface_interpolator(road_data)
    center_z_profile = CarComponents.opencrg_curved_center_z_interpolator(road_data)
    center_heading_profile = CarComponents.opencrg_curved_heading_interpolator(road_data)
    println("Road mesh: $(size(road_data.X)), x range $(extrema(road_data.x_axis)), z range $(extrema(road_data.Z))")

    println("Compiling powered country-road car...")
    model, sys = compile_controlled_country_model(road_surface, center_z_profile, center_heading_profile; target_speed = target_speed)
    println("Compiled: $(length(unknowns(sys))) unknowns, $(length(equations(sys))) equations")

    println("Finding standing-still baseline...")
    baseline = static_baseline(sys, road_surface, center_z_profile, center_heading_profile; settle_time = static_settle_time, start_x = start_x)

    println("Simulating moving run...")
    sol = simulate_country_road(sys, road_surface, center_z_profile, center_heading_profile; speed, target_speed, tstop, save_dt = 1 / framerate, baseline, start_x)
    println("Simulation retcode: $(sol.retcode), samples: $(length(sol.t))")

    println("Building dashboard...")
    fig, time_obs, update_markers! = build_dashboard(model, sys, sol, road_data, road_surface, center_z_profile, center_heading_profile, baseline, target_speed)
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
