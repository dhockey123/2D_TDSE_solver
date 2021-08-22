include("Solvers_2D.jl")
include("utils.jl")
include("potential_fields2D.jl")
using GLMakie

Base.@kwdef mutable struct Params
    Nx::Int64
    Ny::Int64 = Nx
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
    t_max::Float64
    
    dx::Float64             = (abs(x_min)+abs(x_max))/(Nx-1)
    dy::Float64             = (abs(y_min)+abs(y_max))/(Ny-1)
    dt::Float64             = 0.5*dx^2
    t::Array{Float64}       = collect(0:dt:t_max)
    x::Array{Float64}       = collect(x_min:dx:x_max)
    y::Array{Float64}       = collect(y_min:dy:y_max)

    kx_0::Float64
    ky_0::Float64
    σx::Float64
    σy::Float64
    x_0::Float64
    y_0::Float64
    
    # Init Ψ at t=0
    ψ::Array{ComplexF64,2} = ( 1 / sqrt(2π*σx*σy) ).*[exp(-((x-x_0)^2) / (2*σx^2))*exp(-((y-y_0)^2) / (2*σy^2))*exp(im*kx_0*x + im*ky_0*y) for x in x, y in y]
    ψ_copy::Array{ComplexF64,2} = ( 1 / sqrt(2π*σx*σy) ).*[exp(-((x-x_0)^2) / (2*σx^2))*exp(-((y-y_0)^2) / (2*σy^2))*exp(im*kx_0*x + im*ky_0*y) for x in x, y in y]
    
    potential_field::Array{Float64,2} = zeros(Int64, Nx, Ny)
end

fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),resolution = (1500, 1000))
ax1 = Axis(fig[1:10, 1:10], title = "ℝ(ψ)",tellwidth=true)

########################################################################################
################################ INTERACTION ###########################################

# Link sliders with initial params.

N_slider       = labelslider!(fig, "N:", 11:2:501; format = x -> "$(x)")
x_0_slider     = labelslider!(fig, "xᵢ:", -50:1:50; format = x -> "$(x)")
y_0_slider     = labelslider!(fig, "yᵢ :", -50:1:50; format = x -> "$(x)")
kx_slider      = labelslider!(fig, "kx :", -15:0.2:15; format = x -> "$(x)")
ky_slider      = labelslider!(fig, "ky :", -15:0.2:15; format = x -> "$(x)")
σx_slider      = labelslider!(fig, "σx:", 0.1:0.1:5; format = x -> "$(x)")
σy_slider      = labelslider!(fig, "σy:", 0.1:0.1:5; format = x -> "$(x)")
t_max_slider   = labelslider!(fig, "t:", 1:0.1:10; format = x -> "$(x)")
x_range        = IntervalSlider(fig, range = -50:1:50, startvalues = (-10, 10), tellwidth=true)
y_range        = IntervalSlider(fig, range = -50:1:50, startvalues = (-10, 10), tellwidth=false)

barrier_layout       = Label(fig, "----- Barrier -----")
pos_slider           = labelslider!(fig, "pos:", -50:0.1:50; format = x -> "$(x)")
barrier_width_slider = labelslider!(fig, "width:", 0:0.1:10; format = x -> "$(x)")
slit_size_slider     = labelslider!(fig, "gap:", 0:0.05:10; format = x -> "$(x)")
cylinder_layout      = Label(fig, "----- Cylinder -----")
cylinder_x           = labelslider!(fig, "x:", -50:0.1:50; format = x -> "$(x)")
cylinder_y           = labelslider!(fig, "y:", -50:0.1:50; format = x -> "$(x)")
cylinder_radius      = labelslider!(fig, "r:", 0:0.1:50; format = x -> "$(x)")
cylinder_magnitude   = labelslider!(fig, "|U|:", -500:10:500; format = x -> "$(x)")

# Toggles for changing numerical solver and potential field type.

toggle = Toggle(fig[10,11:15])
fig[10,11:15] = grid!(hcat(Label(fig, "Initialise") , toggle ,  Label(fig, "Animate") ))
potentials_menu = Menu(fig[9,11:15], options = ["Oscillator", "Barrier", "Cylinder"])
solver_menu = Menu(fig[9,11:15], options = ["Pseudospec_FFT", "ADI_solver"])
fig[9,11:15]=hgrid!(potentials_menu, solver_menu)

# Initialise params with slider values.

init_params = @lift Params(Nx=($(N_slider.slider.value)), 
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

Δx = @lift round((($(x_range.interval)[2]-$(x_range.interval)[1]))/$(N_slider.slider.value), digits=4)
Δx_label = Label(fig,@lift string("Δx=", $(Δx)))
Δy_label = Label(fig,@lift string("Δy=", $(Δx)))
Δt_label = Label(fig,@lift string("Δt=", round(0.5*$(Δx)^2, digits = 4)))

# Constructs Layout for sliders and stuff. 

fig[1:7,11:15] = vgrid!(hgrid!(Δx_label, Δy_label, Δt_label),
                        hgrid!(Label(fig, "x,x"), x_range, Label(fig, @lift string($(x_range.interval)))),
                        hgrid!(Label(fig, "y,y"), y_range, Label(fig, @lift string($(y_range.interval)))),
                        N_slider.layout,
                        x_0_slider.layout, y_0_slider.layout, 
                        kx_slider.layout, ky_slider.layout, 
                        σx_slider.layout, σy_slider.layout, 
                        t_max_slider.layout, 
                        barrier_layout,
                        pos_slider.layout,
                        hgrid!(slit_size_slider.layout, barrier_width_slider.layout), 
                        cylinder_layout,
                        cylinder_x.layout, cylinder_y.layout, 
                        hgrid!(cylinder_radius.layout, cylinder_magnitude.layout) )

#################################################################################################################
######################################### INITIAL PLOTS #####################################################
# 

draw_barrier = @lift Utils.draw_U_barrier( ($(init_params)), ($(pos_slider.slider.value)), ($(barrier_width_slider.slider.value)), 
                                                ($(slit_size_slider.slider.value)) )                                                                  

heatmap!(ax1,@lift($(init_params).x), @lift($(init_params).y),  @lift(real.($(init_params).ψ)), colormap=:corkO)
lines!(ax1, (@lift ($(draw_barrier)[1])), color=:black)
lines!(ax1, (@lift ($(draw_barrier)[2])), color=:black)
lines!(ax1, (@lift ($(draw_barrier)[3])), color=:black)
lines!(ax1, @lift Circle(Point2f0(($(cylinder_x.slider.value)), ($(cylinder_y.slider.value))), ($(cylinder_radius.slider.value)) ))
@lift limits!(ax1, $(x_range.interval)[1], $(x_range.interval)[2], $(y_range.interval)[1], $(y_range.interval)[2])

##################################################################################

function animate_plot(draw_barrier)

    ax_abs_2d  = Axis(fig[1:5, 1:5], title = "ℝ(ψ)",tellwidth=true)
    ax_real_2d = Axis(fig[6:10, 1:5], title = "ℝ(ψ)",tellwidth=true)
    ax_real_3d = Axis3(fig[1:5, 6:10], title = "𝕀(Ψ)",tellwidth=true)
    ax_abs_3d  = Axis3(fig[6:10, 6:10], title = "|Ψ|²",tellwidth=true)
    
    # Determines which potential field will be constructed 

    if potentials_menu.selection[] == "Barrier"
        init_params[].potential_field = Potentials2D.doubleslit_potential(init_params[], pos_slider.slider.value[],
                                                                          barrier_width_slider.slider.value[], slit_size_slider.slider.value[])
        lines!(ax_abs_2d, draw_barrier[][1], color=:black)
        lines!(ax_abs_2d, draw_barrier[][2], color=:black)
        lines!(ax_abs_2d, draw_barrier[][3], color=:black)

    elseif potentials_menu.selection[] == "Cylinder"
        init_params[].potential_field = Potentials2D.cylinder_potential(init_params[], cylinder_x.slider.value[], cylinder_y.slider.value[],
                                                                        cylinder_radius.slider.value[], cylinder_magnitude.slider.value[])

        lines!(ax_abs_2d, Circle(Point2f0(cylinder_x.slider.value[], cylinder_y.slider.value[]), cylinder_radius.slider.value[]), color=:black)       

    elseif potentials_menu.selection[] == "Oscillator"
        init_params[].potential_field = Potentials2D.oscillator_2D(init_params[], 2)
    end

    # Initialise constructor in potential_fields.jl

    Solver_par = Solver2D.Solver_Params(par = init_params[])

    # Initialize observable nodes for updating plot coordinates

    ψ_real = Node(real.(init_params[].ψ))
    ψ_imag = Node(imag.(init_params[].ψ))
    ψ_abs  = Node(abs.(init_params[].ψ))

    heatmap!(ax_abs_2d, init_params[].x, init_params[].y, ψ_abs, colormap=:corkO)
    heatmap!(ax_real_2d, init_params[].x, init_params[].y, ψ_real, colormap=:corkO)
    surface!(ax_real_3d, init_params[].x, init_params[].y, ψ_imag,axis=(type=Axis3,), colormap=:corkO)
    surface!(ax_abs_3d, init_params[].x, init_params[].y, ψ_abs,axis=(type=Axis3,), colormap=:corkO)

    record(fig, "test.gif") do io
        for i in init_params[].t

            if solver_menu.selection[] == "Pseudospec_FFT" 
                init_params[].ψ = Solver2D.Pseudospec_FFT(init_params[].ψ, Solver_par)
            elseif solver_menu.selection[] == "ADI_solver"
                init_params[].ψ = Solver2D.ADI_solver(init_params[].ψ, Solver_par)
            end
            
            init_params[].ψ = Solver2D.ADI_solver(init_params[].ψ, Solver_par)

            ψ_real[] = real.(init_params[].ψ) 
            ψ_imag[] = imag.(init_params[].ψ) 
            ψ_abs[]  = abs.(init_params[].ψ) 

            # ax_abs_3d.title = "|ψ|²     t = $(round(i, digits=4))"
            recordframe!(io)
        end
    end
end

# Toggle to start animation

@lift begin
    if $(toggle.active) == true
        empty!(ax1); delete!(ax1)
        @async animate_plot(draw_barrier)
    end
end
