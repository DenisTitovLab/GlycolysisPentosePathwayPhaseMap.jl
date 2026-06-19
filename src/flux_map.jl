# Node layout (data units) mirroring tmp/glycolysis-PPP_figure.png and the approved mock
# (docs/superpowers/specs/2026-06-11-bibi-reference-flux-map-mock.html): glycolysis spine
# (left), oxidative PPP (top row), non-oxidative PPP as a 3x2 sugar grid. GAP and F6P are
# each ONE physical pool but drawn twice — on the glycolysis spine (:GAP/:F6P) and in the
# non-ox grid (:GAP_p/:F6P_p) — so cycle closure reads without long crossing arrows. The
# duplication is a DISPLAY device (matched halo + label), NOT two model pools.
const _PPP_NODES = Dict(
    # glycolysis spine
    :Glucose => (0.85, 5.25), :G6P => (0.85, 4.40), :F6P => (0.85, 3.30),
    :F16BP => (0.85, 2.20), :DHAP => (0.85, 1.15), :GAP => (2.10, 1.70),
    :LOWER => (2.10, 0.80),
    # oxidative PPP (top row)
    :PGLn => (2.35, 4.40), :PGA => (3.75, 4.40), :Ru5P => (5.05, 4.40),
    # non-oxidative PPP grid: col-1 X5P/GAP_p/F6P_p, col-2 R5P/S7P/E4P
    :X5P => (3.80, 3.20), :R5P => (5.60, 3.20),
    :GAP_p => (3.86, 2.05), :S7P => (5.60, 2.05),
    :F6P_p => (3.80, 0.90), :E4P => (5.60, 0.90),
    :NUCLEO => (6.65, 3.20),
)
const _PPP_LABELS = Dict(
    :Glucose => "Glucose", :G6P => "G6P", :F6P => "F6P", :F16BP => "F16BP",
    :DHAP => "DHAP", :GAP => "GAP", :LOWER => "lower glyc.",
    :PGLn => "PGLn", :PGA => "PGA", :Ru5P => "Ru5P",
    :X5P => "X5P", :R5P => "R5P", :GAP_p => "GAP", :S7P => "S7P",
    :F6P_p => "F6P", :E4P => "E4P", :NUCLEO => "R5P sink",
)
# (dx, dy, align) placing each label beside its dot.
const _PPP_LABEL_OFF = Dict(
    :Glucose => (0.0, 0.28, (:center, :bottom)),
    :G6P => (-0.30, 0.0, (:right, :center)), :F6P => (-0.30, 0.0, (:right, :center)),
    :F16BP => (-0.34, 0.0, (:right, :center)), :DHAP => (-0.30, 0.0, (:right, :center)),
    :GAP => (0.0, 0.24, (:center, :bottom)), :LOWER => (0.0, -0.30, (:center, :top)),
    :PGLn => (0.0, 0.26, (:center, :bottom)), :PGA => (0.0, 0.26, (:center, :bottom)),
    :Ru5P => (0.0, 0.26, (:center, :bottom)),
    :X5P => (-0.30, 0.0, (:right, :center)), :R5P => (0.0, -0.24, (:center, :top)),
    :GAP_p => (-0.20, -0.22, (:center, :top)), :S7P => (0.30, 0.0, (:left, :center)),
    :F6P_p => (-0.30, 0.0, (:right, :center)), :E4P => (0.30, 0.0, (:left, :center)),
    :NUCLEO => (0.0, 0.24, (:center, :bottom)),
)
# Background region tints (xlo, xhi, ylo, yhi, colour) echoing the template's boxes.
const _PPP_GROUPS = (
    (0.22, 1.54, 0.27, 5.52, (:gold,      0.05)),   # glycolysis
    (1.58, 5.60, 4.08, 4.70, (:tomato,    0.07)),   # oxidative PPP
    (2.55, 6.17, 0.20, 3.72, (:gray50,    0.12)),   # non-oxidative PPP
)
# Shared-pool nodes (drawn with a pale halo + bold label).
const _PPP_SHARED = (:F6P, :GAP, :F6P_p, :GAP_p)

# Enzyme-square style per region (glycolysis squares smaller, de-emphasised). The non-ox
# squares use the same gray scheme as their background box (`gray50` tint).
const _PPP_SQUARE_STYLE = Dict(
    :glyc  => (stroke = :goldenrod4, txt = :gray25, fill = :white,  fs = 14, hw = 0.30, hh = 0.16),
    :ox    => (stroke = :firebrick,  txt = :firebrick, fill = :white, fs = 15, hw = 0.40, hh = 0.19),
    :nonox => (stroke = :gray50,     txt = :gray30, fill = :white,  fs = 15, hw = 0.40, hh = 0.19),
)

# Reactions. Each carries ONE signed flux (its legs colour/scale/flip together), an enzyme
# `square` (centre, in front of the legs), a `region` (square style), and 1–2 `legs`:
#   (from, to, ctrl)  ctrl === nothing -> straight arrow;
#                     ctrl::Tuple      -> quadratic curve whose MIDPOINT sits at ctrl
#                                         (paired bi-bi legs share the square centre, so they
#                                          meet behind the square). Backbone pairings come
#                                          from PPP.jl isotope_tracing.jl atom maps.
const _PPP_REACTIONS = (
    (key = :V_HK1,   label = "HK1",   region = :glyc, square = (0.85, 4.82),
        legs = ((:Glucose, :G6P, nothing),)),
    (key = :V_GPI,   label = "GPI",   region = :glyc, square = (0.85, 3.85),
        legs = ((:G6P, :F6P, nothing),)),
    (key = :V_PFKP,  label = "PFKP",  region = :glyc, square = (0.85, 2.75),
        legs = ((:F6P, :F16BP, nothing),)),
    (key = :V_ALDO,  label = "ALDO",  region = :glyc, square = (0.85, 1.68),
        legs = ((:F16BP, :DHAP, nothing), (:F16BP, :GAP, nothing))),
    (key = :V_TPI,   label = "TPI",   region = :glyc, square = (1.50, 1.43),
        legs = ((:DHAP, :GAP, nothing),)),
    (key = :V_GAPDH, label = "GAPDH", region = :glyc, square = (2.10, 1.25),
        legs = ((:GAP, :LOWER, nothing),)),
    (key = :V_G6PD,  label = "G6PD",  region = :ox, square = (1.60, 4.40),
        legs = ((:G6P, :PGLn, nothing),)),
    (key = :V_PGLS,   label = "PGLS",   region = :ox, square = (3.05, 4.40),
        legs = ((:PGLn, :PGA, nothing),)),
    (key = :V_PGD,   label = "PGD",   region = :ox, square = (4.40, 4.40),
        legs = ((:PGA, :Ru5P, nothing),)),
    (key = :V_RPE,   label = "RPE",   region = :nonox, square = (4.42, 3.80),
        legs = ((:Ru5P, :X5P, nothing),)),
    (key = :V_RPI,   label = "RPI",   region = :nonox, square = (5.32, 3.80),
        legs = ((:Ru5P, :R5P, nothing),)),
    (key = :V_R5Pase, label = "R5Pase", region = :nonox, square = (6.12, 3.00),
        legs = ((:R5P, :NUCLEO, nothing),)),                       # square below the arrow
    (key = :V_TKT_Rxn1, label = "TKT", region = :nonox, square = (4.70, 2.62),
        legs = ((:X5P, :GAP_p, (4.70, 2.62)), (:R5P, :S7P, (4.70, 2.62)))),
    (key = :V_TA, label = "TA", region = :nonox, square = (4.70, 1.48),
        legs = ((:GAP_p, :F6P_p, (4.70, 1.48)), (:S7P, :E4P, (4.70, 1.48)))),
    (key = :V_TKT_Rxn2, label = "TKT R2", region = :nonox, square = (4.70, 0.68),
        legs = ((:E4P, :F6P_p, (4.70, 0.68)), (:X5P, :GAP_p, (3.54, 2.62))),
        extra = ((3.54, 2.62), "TKT R2")),   # label the separated X5P<->GAP arm of R2
)
const _PPP_FLUX_KEYS = Tuple(unique(r.key for r in _PPP_REACTIONS))
const _LEG_GAP = 0.14   # data-unit gap between an arrow end and its node dot

# Draw an enzyme square (opaque, in front of its legs) with the label centred inside.
function _draw_enzyme_square!(ax, pos, label, region)
    s = _PPP_SQUARE_STYLE[region]
    poly!(ax, Point2f[(pos[1] - s.hw, pos[2] - s.hh), (pos[1] + s.hw, pos[2] - s.hh),
                      (pos[1] + s.hw, pos[2] + s.hh), (pos[1] - s.hw, pos[2] + s.hh)];
          color = s.fill, strokecolor = s.stroke, strokewidth = 1)
    text!(ax, pos[1], pos[2]; text = label, align = (:center, :center),
          fontsize = s.fs, font = :italic, color = s.txt)
    ax
end

# Trimmed cubic Bézier polyline from node centre A (tail) to B (head) bowing through M (the
# enzyme square centre), with `gA`/`gB` data-unit gaps trimmed off each end so the shaft starts/
# stops shy of the node dots. Returns the point list (tuples). The handle length s = 4/3 makes
# the curve pass through M at t = 0.5 so paired bi-bi legs meet behind the square.
function _curved_pts(A, B, M, gA, gB)
    s = 4 / 3
    P1 = (A[1] + s * (M[1] - A[1]), A[2] + s * (M[2] - A[2]))
    P2 = (B[1] + s * (M[1] - B[1]), B[2] + s * (M[2] - B[2]))
    bez(t) = (u = 1 - t; (u^3*A[1] + 3u^2*t*P1[1] + 3u*t^2*P2[1] + t^3*B[1],
                          u^3*A[2] + 3u^2*t*P1[2] + 3u*t^2*P2[2] + t^3*B[2]))
    pts = [bez(t) for t in range(0, 1; length = 64)]
    i0 = findfirst(p -> hypot(p[1] - A[1], p[2] - A[2]) >= gA, pts)
    i1 = findlast(p  -> hypot(p[1] - B[1], p[2] - B[2]) >= gB, pts)
    i0 === nothing && (i0 = 1)
    (i1 === nothing || i1 <= i0) && (i1 = length(pts))
    return pts[i0:i1]
end

# Pure geometry for ONE flux leg. Returns the shaft polyline (as Point2f), colour (reverse =
# steel blue, forward = tomato), width/tip pixel sizes, and the arrowhead position + rotation.
# A negative flux flips tail<->head so the arrowhead points the reverse way. LOWER is a text-only
# label (no dot), so its end takes no gap.
function _leg_drawspec(from, to, ctrl, v, lw, tip)
    col = v < 0 ? :steelblue : :tomato
    nf, nt = from, to
    A = _PPP_NODES[from]; B = _PPP_NODES[to]
    if v < 0
        (A, B) = (B, A); (nf, nt) = (nt, nf)
    end
    gA = nf === :LOWER ? 0.0 : _LEG_GAP
    gB = nt === :LOWER ? 0.0 : _LEG_GAP
    if ctrl === nothing
        L = hypot(B[1] - A[1], B[2] - A[2])
        if L > gA + gB
            ux, uy = (B[1] - A[1]) / L, (B[2] - A[2]) / L
            A = (A[1] + ux * gA, A[2] + uy * gA)
            B = (B[1] - ux * gB, B[2] - uy * gB)
        end
        pts = [A, B]
    else
        pts = _curved_pts(A, B, ctrl, gA, gB)
    end
    he = pts[end]; ref = pts[max(1, lastindex(pts) - 1)]
    headrot = atan(he[2] - ref[2], he[1] - ref[1]) - pi/2   # utriangle points +y at rotation 0
    shaft = Point2f[Point2f(p[1], p[2]) for p in pts]
    return (shaft = shaft, col = col, lw = lw, tip = tip,
            head = Point2f(he[1], he[2]), headrot = headrot)
end

# Fixed flat list of legs (built once at load): the leg COUNT never changes, so the explorer can
# build exactly one shaft + one arrowhead plot per entry and only update their data on reselect.
const _PPP_LEGS = [(key = r.key, region = r.region, from = leg[1], to = leg[2], ctrl = leg[3])
                   for r in _PPP_REACTIONS for leg in r.legs]

"""
    draw_ppp_flux_map!(ax, fluxes_obs::Observable; title = "", maxflux = nothing)
    draw_ppp_flux_map!(ax, fluxes::NamedTuple; title = "", maxflux = nothing)

Draw the glycolysis+PPP wiring. The static structure (region tints, enzyme squares, shared-pool
halos, node dots, labels) is drawn once. The reaction legs are driven by `fluxes_obs`: a single
`lift` recomputes per-leg geometry (colour = sign, width proportional to |flux|, arrowhead
direction = sign) and each leg's shaft (`lines!`) + arrowhead (`scatter!`) updates IN PLACE — no
plot is added or removed on reselect. Glycolysis and PPP self-normalize on separate width scales
(PPP max width = half glycolysis) unless `maxflux` overrides the glycolysis reference. `title` may
be a String or an `Observable{String}`. The `NamedTuple` method wraps the value in a constant
`Observable`, so the static CairoMakie path evaluates each `lift` exactly once.
"""
function draw_ppp_flux_map!(ax, fluxes_obs::Observable; title = "", maxflux = nothing)
    hidedecorations!(ax); hidespines!(ax)
    if title isa Observable
        on(title; update = true) do t
            ax.title = t
        end
    else
        ax.title = title
    end
    xlims!(ax, 0.0, 7.1); ylims!(ax, 0.1, 5.7)

    # static background region tints
    for (xlo, xhi, ylo, yhi, col) in _PPP_GROUPS
        poly!(ax, Point2f[(xlo, ylo), (xhi, ylo), (xhi, yhi), (xlo, yhi)]; color = col)
    end
    # Dynamic per-leg draw specs: ONE lift, so gmax/pmax are computed once per update. Two
    # independent width scales keep the PPP branch legible even when glycolytic flux dwarfs it:
    # glycolysis normalizes to its OWN max (max width 14px), the PPP (oxidative + non-oxidative)
    # to its OWN max but at HALF the range (max width 7px). So the widest PPP arrow is half the
    # widest glycolytic one — a visual cue that the PPP carries less flux while still showing
    # within-PPP contrast. `maxflux`, if given, overrides the glycolysis reference.
    legdata = lift(fluxes_obs) do f
        _gm(pred) = begin
            m = maximum((abs(getfield(f, r.key)) for r in _PPP_REACTIONS if pred(r.region));
                        init = 0.0) do v
                isnan(v) ? 0.0 : v
            end
            (m == 0 || isnan(m)) ? 1.0 : m
        end
        gmax = maxflux === nothing ? _gm(==(:glyc)) : maxflux
        gmax = (gmax == 0 || isnan(gmax)) ? 1.0 : gmax
        pmax = _gm(!=(:glyc))
        map(_PPP_LEGS) do L
            v = getfield(f, L.key)
            av = isnan(v) ? 0.0 : abs(v)
            lw = L.region === :glyc ? 1.0 + 13.0 * av / gmax : 1.0 + 6.0 * av / pmax
            tip = 8.0 + 2.2 * lw
            _leg_drawspec(L.from, L.to, L.ctrl, v, lw, tip)
        end
    end

    # one shaft (lines!) + one arrowhead (scatter!) per leg, built ONCE, lifted by index. Drawn
    # AFTER the region tints and BEFORE the squares so legs stay behind the enzyme squares.
    for j in eachindex(_PPP_LEGS)
        col = lift(d -> d[j].col, legdata)
        lines!(ax, lift(d -> d[j].shaft, legdata); color = col,
               linewidth = lift(d -> d[j].lw, legdata))
        scatter!(ax, lift(d -> d[j].head, legdata); marker = :utriangle,
                 markersize = lift(d -> d[j].tip, legdata),
                 rotation = lift(d -> d[j].headrot, legdata), color = col)
    end

    # static enzyme squares on top (+ any `extra` label, e.g. a reaction's separated arm)
    for r in _PPP_REACTIONS
        _draw_enzyme_square!(ax, r.square, r.label, r.region)
        ex = get(r, :extra, nothing)
        ex === nothing || _draw_enzyme_square!(ax, ex[1], ex[2], r.region)
    end
    # static shared-pool halos, then nodes (LOWER is text-only), then labels
    for name in _PPP_SHARED
        (x, y) = _PPP_NODES[name]
        scatter!(ax, [x], [y]; markersize = 26, color = (:gold, 0.25))
    end
    for (name, (x, y)) in _PPP_NODES
        name === :LOWER && continue
        scatter!(ax, [x], [y]; markersize = 12, color = :white,
                 strokecolor = :black, strokewidth = 1.8)
    end
    for (name, (x, y)) in _PPP_NODES
        off = get(_PPP_LABEL_OFF, name, (0.0, 0.42, (:center, :bottom)))
        text!(ax, x + off[1], y + off[2]; text = get(_PPP_LABELS, name, string(name)),
              align = off[3], fontsize = 15,
              font = (name in _PPP_SHARED ? :bold : :regular))
    end
    return ax
end

draw_ppp_flux_map!(ax, fluxes::NamedTuple; title = "", maxflux = nothing) =
    draw_ppp_flux_map!(ax, Observable(fluxes); title = title, maxflux = maxflux)
