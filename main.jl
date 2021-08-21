include("Solvers.jl")
include("utils.jl")
include("potential_fields2D.jl")
using GLMakie

Base.@kwdef mutable struct Params
    Nx::Int64
    Ny::Int64 = Nx
    Nt::Int64
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
    t_max::Float64
    
    dt::Float64             = t_max/(Nt-1)
    dx::Float64             = (abs(x_min)+abs(x_max))/(Nx-1)
    dy::Float64             = (abs(y_min)+abs(y_max))/(Ny-1)
    t::Array{Float64}       = collect(0:dt:t_max)
    x::Array{Float64}       = collect(x_min:dx:x_max)
    y::Array{Float64}       = collect(y_min:dy:y_max)

    kx_0::Float64
    ky_0::Float64
    σx::Float64
    σy::Float64
    x_0::Float64
    y_0::Float64
    
    #Huzzahhh!
    ψ::Array{ComplexF64,2} = ( 1 / sqrt(2π*σx*σy) ).*[exp(-((x-x_0)^2) / (2*σx^2))*exp(-((y-y_0)^2) / (2*σy^2))*exp(im*kx_0*x + im*ky_0*y) for x in x, y in y]
    potential_field::Array{Float64,2} = zeros(Int64, Nx, Ny)
end

fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),resolution = (1200, 900))
ax1 = Axis(fig[1:10, 1:10], title = "ℝ(ψ)",tellwidth=true)
s1, interval_layout = layoutscene(padding=0)
s2, barrier_layout  = layoutscene(padding=0)

N_slider       = labelslider!(fig, "N:", 49:2:501; format = x -> "$(x)")
x_0_slider     = labelslider!(fig, "x_0:", -50:1:50; format = x -> "$(x)")
y_0_slider     = labelslider!(fig, "y_0:", -50:1:50; format = x -> "$(x)")
kx_slider      = labelslider!(fig, "kx:", -15:0.2:15; format = x -> "$(x)")
ky_slider      = labelslider!(fig, "ky:", -15:0.2:15; format = x -> "$(x)")
σx_slider      = labelslider!(fig, "σx:", 0.1:0.1:5; format = x -> "$(x)")
σy_slider      = labelslider!(fig, "σy:", 0.1:0.1:5; format = x -> "$(x)")
Nt_slider      = labelslider!(fig, "Nt:", 10:1:512; format = x -> "$(x)")
t_max_slider   = labelslider!(fig, "t:", 1:0.1:10; format = x -> "$(x)")

barrier_layout[1,1]  = Label(fig, "----- Barrier -----")
pos_slider           = labelslider!(fig, "pos:", -50:0.1:50; format = x -> "$(x)")
slit_barrier_slider  = labelslider!(fig, "width:", 0.1:0.1:10; format = x -> "$(x)")
slit_size_slider     = labelslider!(fig, "gap:", 0.1:0.05:10; format = x -> "$(x)")

interval_layout[2,1] = Label(fig, "x,x")
interval_layout[2,2] = x_range = IntervalSlider(fig, range = -50:1:50, startvalues = (-10, 10), tellwidth=true)
interval_layout[2,3] = Label(fig, @lift string($(x_range.interval)))
interval_layout[3,1] = Label(fig, "y,y")
interval_layout[3,2] = y_range = IntervalSlider(fig, range = -50:1:50, startvalues = (-10, 10), tellwidth=true)
interval_layout[3,3] = Label(fig, @lift string($(y_range.interval)))

fig[1:8,11:14] = vgrid!(interval_layout,N_slider.layout,
                    x_0_slider.layout, y_0_slider.layout, 
                    kx_slider.layout, ky_slider.layout, 
                    σx_slider.layout, σy_slider.layout, 
                    Nt_slider.layout, t_max_slider.layout, 
                    barrier_layout,
                    pos_slider.layout, slit_barrier_slider.layout,
                    slit_size_slider.layout)

toggle = Toggle(fig[10,11:14])
fig[10,11:14] = grid!(hcat(Label(fig, "Initialise") , toggle ,  Label(fig, "Animate") ))

init_params = @lift Params(Nx=($(N_slider.slider.value)), 
                            Nt = ($(Nt_slider.slider.value)),
                            x_min=$(x_range.interval)[1], 
                            x_max=$(x_range.interval)[2],
                            y_min=$(y_range.interval)[1], 
                            y_max=$(y_range.interval)[2], 
                            t_max=($(t_max_slider.slider.value)),
                            kx_0=($(kx_slider.slider.value)), 
                            ky_0=($(ky_slider.slider.value)),
                            σx=($(σx_slider.slider.value)), 
                            σy=($(σy_slider.slider.value)),
                            x_0=($(x_0_slider.slider.value)),
                            y_0=($(y_0_slider.slider.value)))

A = @lift Utils.draw_potential_slits_barrier(($(pos_slider.slider.value)), ($(slit_barrier_slider.slider.value)), 
                                                ($(slit_size_slider.slider.value)), ($(init_params)))  
                                                
Δx = @lift round((($(x_range.interval)[2]-$(x_range.interval)[1]))/$(N_slider.slider.value), digits=4)
Δy = @lift round((($(y_range.interval)[2]-$(y_range.interval)[1]))/$(N_slider.slider.value), digits=4)
Δt = @lift round($(t_max_slider.slider.value)/$(Nt_slider.slider.value), digits=4)
Δt_label = Label(fig,@lift string("Δt=", $(Δt)))
Δx_label = Label(fig,@lift string("Δx=", $(Δx)))
Δy_label = Label(fig,@lift string("Δy=", $(Δy)))

interval_layout[1,1] = Δx_label
interval_layout[1,2] = Δy_label
interval_layout[1,3] = Δt_label
                            
heatmap!(ax1,@lift($(init_params).x), @lift($(init_params).y),  @lift(real.($(init_params).ψ)), colormap=:corkO)
poly!(ax1, (@lift ($(A)[1])), color=:black)
poly!(ax1, (@lift ($(A)[2])), color=:black)
poly!(ax1, (@lift ($(A)[3])), color=:black)

@lift limits!(ax1, $(x_range.interval)[1], $(x_range.interval)[2], $(y_range.interval)[1], $(y_range.interval)[2])

##################################################################################

function animate_plot(A)

    init_params[].potential_field = Potentials2D.make_U_doubleslit(init_params[], pos_slider.slider.value[],
                                                                    slit_barrier_slider.slider.value[], slit_size_slider.slider.value[])
    Solver_par = Solver2D.Solver_Params(par = init_params[])

    ax_abs_2d  = Axis(fig[1:5, 1:5], title = "ℝ(ψ)",tellwidth=true)
    ax_real_2d = Axis(fig[6:10, 1:5], title = "ℝ(ψ)",tellwidth=true)
    ax_real_3d = Axis3(fig[1:5, 6:10], title = "𝕀(Ψ)",tellwidth=true)
    ax_abs_3d  = Axis3(fig[6:10, 6:10], title = "|Ψ|²",tellwidth=true)
    
    ψ_real = Node(real.(init_params[].ψ))
    ψ_imag = Node(imag.(init_params[].ψ))
    ψ_abs = Node(abs.(init_params[].ψ))
    #Utils.set_limits(init_params[], ax3, ax4)
    poly!(ax_abs_2d, A[][1], color=:black)
    poly!(ax_abs_2d, A[][2], color=:black)
    poly!(ax_abs_2d, A[][3], color=:black)
    heatmap!(ax_abs_2d, init_params[].x, init_params[].y, ψ_abs, colormap=:corkO)
    heatmap!(ax_real_2d, init_params[].x, init_params[].y, ψ_real, colormap=:corkO)
    surface!(ax_real_3d, init_params[].x, init_params[].y, ψ_imag,axis=(type=Axis3,), colormap=:corkO)
    surface!(ax_abs_3d, init_params[].x, init_params[].y, ψ_abs,axis=(type=Axis3,), colormap=:corkO)

    record(fig, "test.gif") do io
        for i in init_params[].t

            init_params[].ψ = Solver2D.ADI_solver(init_params[].ψ, Solver_par)

            ψ_real[] = real.(init_params[].ψ) 
            ψ_imag[] = imag.(init_params[].ψ) 
            ψ_abs[]  = abs.(init_params[].ψ) 

            ax_abs_3d.title = "|ψ|²     t = $(round(i, digits=4))"
            recordframe!(io)
        end
    end
end

@lift begin
    if $(toggle.active) == true
        empty!(ax1); delete!(ax1)
        @async animate_plot(A)
    end
end
