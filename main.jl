include("Solvers.jl")
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
end

fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),resolution = (1400, 1000))
ax1 = Axis3(fig[1:5, 1:5], title = "ℝ(ψ) ",tellwidth=true)
ax2 = Axis3(fig[6:10, 1:5], title = "ℝ(ψ) ",tellwidth=true)
ax3 = Axis3(fig[1:5, 6:10], title = "ℝ(ψ) ",tellwidth=true)
ax4 = Axis3(fig[6:10, 6:10], title = "ℝ(ψ) ",tellwidth=true)
s1, interval_layout = layoutscene(padding=0)

N_slider     = labelslider!(fig, "N:", 50:2:500; format = x -> "$(x)")
x_0_slider   = labelslider!(fig, "x_0:", -30:1:30; format = x -> "$(x)")
y_0_slider   = labelslider!(fig, "x_0:", -30:1:30; format = x -> "$(x)")
kx_slider    = labelslider!(fig, "kx:", -15:0.2:15; format = x -> "$(x)")
ky_slider    = labelslider!(fig, "ky:", -15:0.2:15; format = x -> "$(x)")
σx_slider    = labelslider!(fig, "σx:", 0.1:0.1:5; format = x -> "$(x)")
σy_slider    = labelslider!(fig, "σy:", 0.1:0.1:5; format = x -> "$(x)")
Nt_slider    = labelslider!(fig, "Nt:", 10:1:512; format = x -> "$(x)")
t_max_slider = labelslider!(fig, "t:", 1:0.1:10; format = x -> "$(x)")

interval_layout[1,1] = Label(fig, "x,x")
interval_layout[1,2] = x_range = IntervalSlider(fig, range = -30:1:30, startvalues = (-10, 10), tellwidth=true)
interval_layout[1,3] = Label(fig, @lift string($(x_range.interval)))
interval_layout[2,1] = Label(fig, "y,y")
interval_layout[2,2] = y_range = IntervalSlider(fig, range = -30:1:30, startvalues = (-10, 10), tellwidth=true)
interval_layout[2,3] = Label(fig, @lift string($(y_range.interval)))

fig[1:8,11:14] = vgrid!(interval_layout,N_slider.layout,
                    x_0_slider.layout, y_0_slider.layout, 
                    kx_slider.layout, ky_slider.layout, 
                    σx_slider.layout, σy_slider.layout, 
                    Nt_slider.layout, t_max_slider.layout)

toggle = Toggle(fig[9:10,11:14])
fig[9:10,11:14] = grid!(hcat(Label(fig, "Initialise") , toggle ,  Label(fig, "Animate") ))

init_params = @lift Params( Nx=200, 
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
                            y_0=($(y_0_slider.slider.value)) )

surface!(ax1, @lift($(init_params).x), @lift($(init_params).y), @lift(real.($(init_params).ψ)),axis=(type=Axis3,))
surface!(ax2, @lift($(init_params).x), @lift($(init_params).y), @lift(imag.($(init_params).ψ)),axis=(type=Axis3,))
surface!(ax3, @lift($(init_params).x), @lift($(init_params).y), @lift(real.($(init_params).ψ)),axis=(type=Axis3,))
surface!(ax4, @lift($(init_params).x), @lift($(init_params).y), @lift(abs2.($(init_params).ψ)),axis=(type=Axis3,))
@lift limits!(ax1, $(x_range.interval)[1], $(x_range.interval)[2], $(y_range.interval)[1], $(y_range.interval)[2], minimum(real.($(init_params).ψ)), maximum(real.($(init_params).ψ)))
@lift limits!(ax2, $(x_range.interval)[1], $(x_range.interval)[2], $(y_range.interval)[1], $(y_range.interval)[2], minimum(imag.($(init_params).ψ)), maximum(imag.($(init_params).ψ)))
@lift limits!(ax3, $(x_range.interval)[1], $(x_range.interval)[2], $(y_range.interval)[1], $(y_range.interval)[2], minimum(abs.($(init_params).ψ)), maximum(abs.($(init_params).ψ)))
@lift limits!(ax4, $(x_range.interval)[1], $(x_range.interval)[2], $(y_range.interval)[1], $(y_range.interval)[2], minimum(real.($(init_params).ψ)), maximum(real.($(init_params).ψ)))

# # #################################################################################

function animate_plot()
    #delete!(ax1); delete!(ax2); delete!(ax3); delete!(ax4) 
    ax1 = Axis3(fig[1:5, 1:5], title = "ℝ(ψ) ",tellwidth=true)
    ax2 = Axis3(fig[6:10, 1:5], title = "ℝ(ψ) ",tellwidth=true)
    ax3 = Axis3(fig[1:5, 6:10], title = "ℝ(ψ) ",tellwidth=true)
    ax4 = Axis3(fig[6:10, 6:10], title = "ℝ(ψ) ",tellwidth=true)

    Solver_par = Solver2D.Solver_Params(par = init_params[])
    real_node = Node(real.(init_params[].ψ))
    imag_node = Node(imag.(init_params[].ψ))
    abs_node = Node(abs.(init_params[].ψ))
    
    surface!(ax1, init_params[].x, init_params[].y, real_node,axis=(type=Axis3,))
    surface!(ax2, init_params[].x, init_params[].y, imag_node,axis=(type=Axis3,))
    surface!(ax3, init_params[].x, init_params[].y, abs_node,axis=(type=Axis3,))
    surface!(ax4, init_params[].x, init_params[].y, real_node,axis=(type=Axis3,))
    
    record(fig, "test.gif") do io
        for i in init_params[].t
            init_params[].ψ = Solver2D.ADI_solver(init_params[].ψ, Solver_par)
            real_node[] = real.(init_params[].ψ) 
            imag_node[] = imag.(init_params[].ψ) 
            abs_node[] = abs.(init_params[].ψ) 
            ax1.title = "|ψ|²     t = $(round(i, digits=4))"
            recordframe!(io)
        end
    end
end


@lift begin
    if $(toggle.active) == true
        delete!(ax1);delete!(ax2);delete!(ax3);delete!(ax4)
        @async animate_plot()
    end
end