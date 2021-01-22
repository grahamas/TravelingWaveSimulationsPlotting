using WilsonCowanModel: wcm_du_defn, wcm_dv_defn, WCMParams, AbstractNullclineParams
using Makie: @lift

function lifted_wcm_param(;
    Aee=1., Aie=1., Aei=1.5, Aii=0.25,
    θef=0.125, θif=0.4, θib=7., τ=0.4, β=50.,
    decaye=1., decayi=1.)
    @lift WCMParams(;
            Aee=$(Node(Aee)), Aie=$(Node(Aie)), Aei=$(Node(Aei)), Aii=$(Node(Aii)),
            θef=$(Node(θef)), θif=$(Node(θif)), θib=$(Node(θib)), τ=$(Node(τ)), β=$(Node(β)),
            decaye=$(Node(decaye)), decayi=$(Node(decayi))
        )
end

# plotting
function nullclines(args...)
    scene, layout = layoutscene()
    layout[1,1] =  nullclines!(scene, args...)
    return (scene, layout)
end

function nullclines!(scene, p::AbstractNullclineParams, dx=0.01)
    nullclines!(scene, p, 0.:dx:1., 0.:dx:1.)
end
function nullclines!(scene, p::AbstractNullclineParams, us::AbstractVector, vs::AbstractVector)

    layout = GridLayout()
    layout[1,1] = ax = LAxis(scene)

    dus = [wcm_du_defn(u, v, p) for u in us, v in vs]
    dvs = [wcm_dv_defn(u, v, p) for u in us, v in vs]
    u_nullclines = lines(contour(us, vs, dus, 0.))
    v_nullclines = lines(contour(us, vs, dvs, 0.))

    for line in u_nullclines
        xs, ys = coordinates(line)
        Makie.lines!(ax, xs, ys, color=:blue)
    end

    for line in v_nullclines
        xs, ys = coordinates(line)
        Makie.lines!(ax, xs, ys, color=:red, linestyle=:dash, linesize=5)
    end
    layout[1,0] = LText(scene, "v", tellheight=false)
    layout[end+1,2] = LText(scene, "u", tellwidth=false)
    return layout
end
