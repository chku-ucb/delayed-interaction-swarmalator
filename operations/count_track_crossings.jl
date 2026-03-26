#!/usr/bin/env julia
# Count threshold crossings for TrackMate XML files produced by Fiji.
# A crossing is recorded when a particle (track) moves across the circle of
# radius `threshold = 0.85 * 133`. Moving from inside→outside counts as +1,
# outside→inside counts as -1.
# Uses Plots.jl for outputs (PNG).

# Adjust the path to environment and utils.jl as needed
script_dir = @__DIR__
code_root = joinpath(script_dir, "..") |> normpath
env_dir = joinpath(code_root, "env")
# -------------------------------------------------
using Pkg
Pkg.activate(env_dir)
include(joinpath(code_root,"utils.jl"))

const THRESHOLD = 130.0
const CENTER_X = 300.0
const CENTER_Y = 200.0
const TRACKS_DIR = "tracks_low_tif/res_track/"
const PLOT_WIDTH = 900
const PLOT_HEIGHT = 420
const MARGIN = 50
const OUT_ANIM_FPS = 12
const OUT_ANIM_SIZE = (600, 600)
const OUT_ANIM_FRAME_SKIP = 5
const OUT_COLOR = :red
const IN_COLOR = :blue
const BOTH_COLOR = :purple
const PATH_ALPHA = 0.5
const BLUE = (31, 119, 180)
const RED = (214, 39, 40)
const BLACK = (51, 51, 51)
const WHITE = (255, 255, 255)
const ENABLE_ANIMATION = true
const DEFAULT_XLABEL = "time (frame)"
const FIGURE_TIME_XLABEL = "time (image timestamp)"

cli_track_filter() = isempty(ARGS) ? nothing : ARGS[1]

"""
    parse_attr(line, key) -> Union{String,Nothing}

Return the value of `key="..."` from the XML line, or `nothing` if absent.
"""
function parse_attr(line::AbstractString, key::AbstractString)
    m = match(Regex("\\b$(key)=\"([^\"]+)\""), line)
    return m === nothing ? nothing : m.captures[1]
end

"""
    require_attr(line, key)

Fetch `key` attribute from `line` or throw a clear error if missing.
"""
function require_attr(line::AbstractString, key::AbstractString)
    val = parse_attr(line, key)
    val === nothing && error("Missing attribute $(key) in line:\n", line)
    return val
end

function find_figure_dir(xml_path::AbstractString)
    base = replace(basename(xml_path), r"\.xml$" => "")
    xml_dir = dirname(xml_path)
    tracks_dir_abs = try
        abspath(TRACKS_DIR)
    catch
        nothing
    end
    candidates = String[
        base,
        joinpath(xml_dir, base),
        joinpath(@__DIR__, base)
    ]
    if tracks_dir_abs !== nothing
        push!(candidates, joinpath(dirname(tracks_dir_abs), base))
    end
    seen = Set{String}()
    for cand in candidates
        cand_abs = abspath(cand)
        cand_abs in seen && continue
        push!(seen, cand_abs)
        if isdir(cand_abs)
            return cand_abs
        end
    end
    return nothing
end

function collect_time_axis(fig_dir::AbstractString, expected_len::Int)
    expected_len == 0 && return nothing
    times = Int[]
    for entry in readdir(fig_dir)
        m = match(r"_time_(\d+)", entry)
        m === nothing && continue
        push!(times, parse(Int, m.captures[1]))
    end
    isempty(times) && return nothing
    sort!(times)
    unique!(times)
    length(times) < expected_len && return nothing
    return times[1:expected_len]
end

function figure_time_axis(xml_path::AbstractString, expected_len::Int)
    fig_dir = find_figure_dir(xml_path)
    fig_dir === nothing && return nothing
    return collect_time_axis(fig_dir, expected_len)
end

function track_suffix_id(path::AbstractString)
    base = basename(path)
    m = match(r"_(\d+)_stack\.xml$", base)
    # Fallback to older pattern or just numeric suffix if _stack not found?
    # For now, let's also try the generic numeric suffix if the specific one fails
    if m === nothing
        m = match(r"(\d+)(?=\_.xml$)", base)
    end
    return m === nothing ? nothing : m.captures[1]
end

function filename_numeric_suffix(path::AbstractString)
    base, _ = splitext(basename(path))
    m = match(r"(\d+)(?!.*\d)", base)
    return m === nothing ? nothing : m.captures[1]
end

function save_plot(path::AbstractString, xvals::AbstractVector{T}, ins::Vector{Int}, outs::Vector{Int},
    xlabel::AbstractString=DEFAULT_XLABEL, id::Union{Nothing,AbstractString}=nothing) where {T<:Real}
    isempty(xvals) && return
    cin = cumsum(ins)
    cout = cumsum(outs)
    net = cout .- cin
    window = 500
    isempty(net) && return
    idxs = 1:window:length(xvals)
    length(idxs) <= 1 && return
    sampled_x = collect(xvals[idxs])
    sampled_net = net[idxs]
    length(sampled_net) <= 1 && return
    delta_vals = diff(sampled_net)
    plot_x = sampled_x[2:end]
    isempty(delta_vals) && return

    p1 = plot(xvals, cin; label="Cumulative In", c=:blue, lw=2, dpi=200)
    plot!(p1, xvals, cout; label="Cumulative Out", c=:red, lw=2)
    savefig(p1, "plot_cumulative.pdf")

    # p = plot(plot_x, delta_vals; label="Δ(cout-cin)", c=:purple, lw=2, xlabel=xlabel,
    #     ylabel="per-frame change", title="Net crossing change over $(window) frames",
    #     size=(PLOT_WIDTH, PLOT_HEIGHT), margin=5Plots.mm)
    # savefig(p, path)
    inferred_id = id === nothing ? filename_numeric_suffix(path) : id
    save_id = inferred_id === nothing ? "" : "_$(inferred_id)"
    # @save "tracks_low_tif/res_track/track_crossings_at_time$(save_id)_500.jld2" time = plot_x flux = delta_vals
end

function write_svg_plot(path::AbstractString, xvals::AbstractVector{T}, ins::Vector{Int}, outs::Vector{Int}) where {T<:Real}
    isempty(xvals) && return
    minf, maxf = minimum(xvals), maximum(xvals)
    maxy = maximum(vcat(ins, outs, [1]))  # avoid zero scale
    w, h, m = PLOT_WIDTH, PLOT_HEIGHT, MARGIN
    xrange = max(1, maxf - minf)

    function build_path(vals::Vector{Int})
        buf = IOBuffer()
        for (idx, f) in enumerate(xvals)
            x = m + (f - minf) / xrange * (w - 2 * m)
            y = h - m - vals[idx] / maxy * (h - 2 * m)
            if idx == 1
                print(buf, "M$x,$y")
            else
                print(buf, " L$x,$y")
            end
        end
        return String(take!(buf))
    end

    pin = build_path(ins)
    pout = build_path(outs)

    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$w" height="$h" viewBox="0 0 $w $h">""")
        println(io, """<rect x="0" y="0" width="$w" height="$h" fill="white" stroke="none"/>""")
        # Axes
        println(io, """<line x1="$m" y1="$(h-m)" x2="$(w-m)" y2="$(h-m)" stroke="#333" stroke-width="1"/>""")
        println(io, """<line x1="$m" y1="$m" x2="$m" y2="$(h-m)" stroke="#333" stroke-width="1"/>""")
        # Paths
        println(io, """<path d="$pin" fill="none" stroke="#1f77b4" stroke-width="2"/>""")
        println(io, """<path d="$pout" fill="none" stroke="#d62728" stroke-width="2"/>""")
        # Legend
        println(io, """<rect x="$(w-180)" y="$m" width="12" height="12" fill="#1f77b4"/><text x="$(w-160)" y="$(m+11)" font-size="12" fill="#000">in</text>""")
        println(io, """<rect x="$(w-180)" y="$(m+20)" width="12" height="12" fill="#d62728"/><text x="$(w-160)" y="$(m+31)" font-size="12" fill="#000">out</text>""")
        println(io, """</svg>""")
    end
end

"""
    load_trackmate(path)

Parse a TrackMate XML file and return:
  - `spots` : Dict{Int,(Int,Float64,Float64)} mapping spot ID → (frame, x, y)
  - `tracks`: Dict{Int,Vector{Int}} mapping track ID → list of spot IDs
Only base functionality is used to avoid external dependencies.
"""
function load_trackmate(path::AbstractString)
    spots = Dict{Int,Tuple{Int,Float64,Float64}}()
    tracks = Dict{Int,Vector{Int}}()

    in_spots = false
    in_tracks = false
    current_track = nothing

    open(path, "r") do io
        for line in eachline(io)
            # Section switches
            if occursin("<AllSpots", line)
                in_spots = true
                continue
            elseif occursin("</AllSpots>", line)
                in_spots = false
                continue
            elseif occursin("<AllTracks", line)
                in_tracks = true
                continue
            elseif occursin("</AllTracks>", line)
                in_tracks = false
                current_track = nothing
                continue
            end

            # Spot collection
            if in_spots && occursin("<Spot ", line)
                sid = parse(Int, require_attr(line, "ID"))
                frame = parse(Int, require_attr(line, "FRAME"))
                x = parse(Float64, require_attr(line, "POSITION_X"))
                y = parse(Float64, require_attr(line, "POSITION_Y"))
                spots[sid] = (frame, x, y)
                continue
            end

            # Track collection
            if in_tracks && occursin("<Track ", line)
                tid = parse(Int, require_attr(line, "TRACK_ID"))
                tracks[tid] = Int[]
                current_track = tid
                continue
            elseif in_tracks && occursin("</Track>", line)
                current_track = nothing
                continue
            elseif in_tracks && occursin("<Edge ", line)
                current_track === nothing && continue
                sid = parse(Int, require_attr(line, "SPOT_SOURCE_ID"))
                tid = parse(Int, require_attr(line, "SPOT_TARGET_ID"))
                push!(tracks[current_track], sid)
                push!(tracks[current_track], tid)
                continue
            end
        end
    end

    # Keep unique spot IDs per track
    tracks = Dict(k => unique(v) for (k, v) in tracks)
    return spots, tracks
end

"""
    count_crossings(spots, spot_ids; threshold, return_events=false)

Compute signed crossings for one track given its spot IDs.
Inside→outside contributes +1, outside→inside contributes -1.
If `return_events` is true, also return a vector of `(frame, dir)` with dir ∈ {+1,-1}.
"""
function count_crossings(spots::Dict{Int,Tuple{Int,Float64,Float64}}, spot_ids::Vector{Int};
    threshold::Real=THRESHOLD, return_events::Bool=false)
    filtered = Tuple{Int,Float64}[]
    for sid in spot_ids
        val = get(spots, sid, nothing)
        val === nothing && continue
        frame, x, y = val
        x = x - CENTER_X
        y = y - CENTER_Y
        r = sqrt(x^2 + y^2)
        push!(filtered, (frame, r))
    end

    if isempty(filtered)
        return return_events ? (Tuple{Int,Int}[], 0) : 0
    end

    sort!(filtered, by=first)  # order by frame
    inside = filtered[1][2] <= threshold
    total = 0
    events = Tuple{Int,Int}[]  # (frame, dir), dir: +1 out, -1 in
    @inbounds for i in 2:length(filtered)
        frame = filtered[i][1]
        r = filtered[i][2]
        new_inside = r <= threshold
        if new_inside != inside
            dir = inside ? 1 : -1  # inside→outside : outside→inside
            total += dir
            push!(events, (frame, dir))
        end
        inside = new_inside
    end
    return return_events ? (events, total) : total
end

const FramePoint = NamedTuple{(:x, :y, :tid),Tuple{Float64,Float64,Int}}

function circle_coords(radius::Real, n::Int=256)
    θ = range(0, 2π; length=n)
    return radius .* cos.(θ), radius .* sin.(θ)
end

function collect_frame_points(spots::Dict{Int,Tuple{Int,Float64,Float64}},
    tracks::Dict{Int,Vector{Int}}, track_ids::Vector{Int})
    frame_points = Dict{Int,Vector{FramePoint}}()
    track_series = Dict{Int,NamedTuple{(:frames, :xs, :ys),Tuple{Vector{Int},Vector{Float64},Vector{Float64}}}}()
    min_frame = typemax(Int)
    max_frame = typemin(Int)
    max_radius = 0.0
    for tid in track_ids
        spot_ids = get(tracks, tid, nothing)
        spot_ids === nothing && continue
        series = Tuple{Int,Float64,Float64}[]
        for sid in spot_ids
            val = get(spots, sid, nothing)
            val === nothing && continue
            frame, x, y = val
            x -= CENTER_X
            y -= CENTER_Y
            push!(series, (frame, x, y))
        end
        isempty(series) && continue
        sort!(series, by=first)
        frames_vec = Int[]
        xs = Float64[]
        ys = Float64[]
        for (frame, x, y) in series
            push!(frames_vec, frame)
            push!(xs, x)
            push!(ys, y)
            frame_points_frame = get!(frame_points, frame, FramePoint[])
            push!(frame_points_frame, (x=x, y=y, tid=tid))
            min_frame = min(min_frame, frame)
            max_frame = max(max_frame, frame)
            max_radius = max(max_radius, sqrt(x^2 + y^2))
        end
        track_series[tid] = (frames=frames_vec, xs=xs, ys=ys)
    end
    if isempty(frame_points)
        return frame_points, nothing, nothing, 0.0, track_series
    end
    return frame_points, min_frame, max_frame, max_radius, track_series
end

function save_crossing_animation(path::AbstractString,
    spots::Dict{Int,Tuple{Int,Float64,Float64}},
    tracks::Dict{Int,Vector{Int}}, tracked_ids::Vector{Int},
    track_roles::Dict{Int,Symbol})
    isempty(tracked_ids) && return nothing
    frame_points, min_frame, max_frame, max_radius, track_series = collect_frame_points(spots, tracks, tracked_ids)
    min_frame === nothing && return nothing

    frames = collect(min_frame:OUT_ANIM_FRAME_SKIP:max_frame)
    if isempty(frames) || frames[end] != max_frame
        push!(frames, max_frame)
    end
    lim = max(THRESHOLD, max_radius) * 1.2
    circle_x, circle_y = circle_coords(THRESHOLD)
    role_colors = Dict(:out => OUT_COLOR, :in => IN_COLOR, :both => BOTH_COLOR)

    anim = @animate for frame in frames
        pts = get(frame_points, frame, FramePoint[])
        p = plot(; xlims=(-lim, lim), ylims=(-lim, lim), size=OUT_ANIM_SIZE,
            aspect_ratio=:equal, legend=:topright, xlabel="x (centered px)",
            ylabel="y (centered px)", title="Crossing tracks — frame $(frame)")
        plot!(p, circle_x, circle_y; color=:gray, lw=1.5, label="threshold")
        for tid in tracked_ids
            ts = get(track_series, tid, nothing)
            ts === nothing && continue
            idx = searchsortedlast(ts.frames, frame)
            idx <= 1 && continue
            role = track_roles[tid]
            plot!(p, ts.xs[1:idx], ts.ys[1:idx]; color=role_colors[role], alpha=PATH_ALPHA,
                lw=1.5, label="")
        end
        if !isempty(pts)
            out_pts = [pt for pt in pts if track_roles[pt.tid] != :in]
            in_pts = [pt for pt in pts if track_roles[pt.tid] == :in]
            if !isempty(out_pts)
                xs = [pt.x for pt in out_pts]
                ys = [pt.y for pt in out_pts]
                role = any(track_roles[pt.tid] == :both for pt in out_pts) ? :both : :out
                scatter!(p, xs, ys; color=role_colors[role], markersize=5, marker=:circle,
                    markerstrokecolor=:white, label=role == :out ? "out" : "out & in")
            end
            if !isempty(in_pts)
                xs_in = [pt.x for pt in in_pts]
                ys_in = [pt.y for pt in in_pts]
                scatter!(p, xs_in, ys_in; color=role_colors[:in], markersize=5, marker=:square,
                    markerstrokecolor=:white, label="in")
            end
        end
    end
    out_path = replace(path, r"\.xml$" => "_crossings.mp4")
    out_path == path && (out_path = path * "_crossings.mp4")
    mp4(anim, out_path; fps=OUT_ANIM_FPS)
    return out_path
end

"""
    process_file(path) -> (per_track::Dict, file_total::Int, axis_vals::Vector{Int}, ins::Vector{Int}, outs::Vector{Int}, anim_path::Union{String,Nothing}, axis_label::String)
"""
function process_file(path::AbstractString)
    spots, tracks = load_trackmate(path)
    per_track = Dict{Int,Int}()
    total = 0
    in_counts = Dict{Int,Int}()   # frame -> entries
    out_counts = Dict{Int,Int}()  # frame -> exits
    outgoing_tracks = Int[]
    incoming_tracks = Int[]

    for (tid, spot_ids) in tracks
        events, c = count_crossings(spots, spot_ids; return_events=true)
        per_track[tid] = c
        total += c
        for (frame, dir) in events
            if dir == 1
                out_counts[frame] = get(out_counts, frame, 0) + 1
            else
                in_counts[frame] = get(in_counts, frame, 0) + 1
            end
        end
        if any(dir == 1 for (_, dir) in events)
            push!(outgoing_tracks, tid)
        end
        if any(dir == -1 for (_, dir) in events)
            push!(incoming_tracks, tid)
        end
    end

    if isempty(spots)
        return per_track, total, Int[], Int[], Int[], nothing, DEFAULT_XLABEL
    end
    min_frame = minimum(first.(values(spots)))
    max_frame = maximum(first.(values(spots)))
    frames = max_frame < min_frame ? Int[] : collect(min_frame:max_frame)
    ins = [get(in_counts, f, 0) for f in frames]
    outs = [get(out_counts, f, 0) for f in frames]
    axis_vals = frames
    axis_label = DEFAULT_XLABEL
    if !isempty(frames)
        custom_axis = figure_time_axis(path, length(frames))
        if custom_axis !== nothing
            axis_vals = custom_axis
            axis_label = FIGURE_TIME_XLABEL
        end
    end
    tracked_ids = sort(unique(vcat(outgoing_tracks, incoming_tracks)))
    track_roles = Dict{Int,Symbol}()
    for tid in tracked_ids
        in_flag = tid in incoming_tracks
        out_flag = tid in outgoing_tracks
        track_roles[tid] = out_flag && in_flag ? :both : (out_flag ? :out : :in)
    end
    anim_path = ENABLE_ANIMATION ? save_crossing_animation(path, spots, tracks, tracked_ids, track_roles) : nothing
    return per_track, total, axis_vals, ins, outs, anim_path, axis_label
end

function main()
    files = sort(filter(f -> endswith(f, ".xml"), readdir(TRACKS_DIR; join=true)))
    isempty(files) && error("No XML files found in $(TRACKS_DIR)")

    target_track = cli_track_filter()
    if target_track !== nothing
        println("Filtering to track ID suffix: ", target_track)
    end
    grand_total = 0
    processed_any = false
    println("Threshold radius: ", THRESHOLD)
    println()
    for f in files
        track_id = target_track === nothing ? nothing : track_suffix_id(f)
        if target_track !== nothing
            track_id === nothing && continue
            track_id == target_track || continue
        end
        per_track, file_total, axis_vals, ins, outs, anim_path, axis_label = process_file(f)
        if isempty(per_track) && isempty(axis_vals)
            println("File: $(basename(f)) (skipped: no TrackMate data?)")
            println("  plot: no data")
            if ENABLE_ANIMATION
                println("  animation: no data")
            else
                println("  animation: disabled")
            end
            continue
        end
        processed_any = true
        grand_total += file_total
        println("File: $(basename(f))")
        println("  tracks: ", length(per_track))
        println("  file crossing sum: ", file_total)
        if !isempty(axis_vals)
            png_path = replace(f, r"\.xml$" => "_in_out.png")
            png_path == f && (png_path = f * "_in_out.png")
            save_plot(png_path, axis_vals, ins, outs, axis_label, track_id)
            println("  plot: ", basename(png_path))
            axis_label == DEFAULT_XLABEL || println("  time axis: figure timestamps")
        else
            println("  plot: no data")
        end
        if !ENABLE_ANIMATION
            println("  animation: disabled")
        elseif anim_path === nothing
            println("  animation: no threshold crossings")
        else
            println("  animation: ", basename(anim_path))
        end
    end
    if target_track !== nothing && !processed_any
        println("No XML files matched track suffix ", target_track)
    end
    println("\nGrand total across all files: ", grand_total)
end

main()
