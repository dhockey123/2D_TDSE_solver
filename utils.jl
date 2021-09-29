module Utils
using GLMakie; limits!, poly

function set_limits(par, ax1, ax2)

    lims=[minimum(par.x), maximum(par.x), minimum(par.y), maximum(par.y), minimum(real.(par.ψ)), maximum(real.(par.ψ))]
    
    A = limits!(ax1, lims[1],lims[2],lims[3],lims[4],lims[5],lims[6])
    B = limits!(ax2, lims[1],lims[2],lims[3],lims[4],lims[5],lims[6])
    # C = limits!(ax3, lims[1],lims[2],lims[3],lims[4],lims[5],lims[6])
    return A, B
end

# Returns coordinate points needed to draw 3 squares on 2d-axis for a double slit barrier parallel to y-axis.
function draw_U_barrier(par, pos, width, slit_size)
    p, w = pos, width   
    y_mid = (par.y_max + par.y_min)/2

    # draw_slit() is used to get coordinates of each of the 3 squares needed for the double slit barrier
    draw_slit(p, w, y1, y2)     = [(p, y1),(p+w, y1),(p+w, y2),(p, y2),(p, y1)]
    draw_top                    = draw_slit(p, w, par.y_max, y_mid+1.5slit_size)
    draw_mid                    = draw_slit(p, w, y_mid-.5slit_size, y_mid+.5slit_size)
    draw_bottom                 = draw_slit(p, w, y_mid-1.5slit_size, par.y_min)
    A = [draw_top, draw_mid, draw_bottom]
    return A
end

# Used to draw circle for the cylindrical potential field
function draw_U_cylinder(pos_x, pos_y, radius)
    return Circle(Point2f0(pos_x, pos_y), radius)
end

# Returns grid coordinate for plotting intensity pattern
function draw_interference_fringes(par, pos_x)
    x_total_length   = abs(par.x_max-par.x_min)
    x_length2pattern = abs(pos_x-par.x_min)
    return Int32(round(x_length2pattern/x_total_length * par.Nx))
end

# Algorithm for getting the position of the peaks in a given intensity pattern.
# if x_3 < x_2 and x_2 > x_1; x_2 = peak.
# fringe : 1D slice of init_params.ψ for where intensity pattern is detected ( e.g init_params.ψ[ : , 100] ) 

# function get_peaks(fringe)
#     peaks = Float64[]
#     for i in 2:1:length(fringe)
#         if fringe[i] > 0.01 && fringe[i] < fringe[i-1] && fringe[i-1] > fringe[i-2]
#              push!(peaks, x[i-1])
#         end
#     end
#     for i in 1:1:length(peaks)
#         if peaks[i] == 7.5
#             peaks[i] = 0
#         elseif peaks[i] > 7.5
#             peaks[i] = round(peaks[i]-7.5, digits=6)
#         else
#             peaks[i] = round(-(7.5-peaks[i]), digits=6)
#         end
#     end
#     return peaks
# end

end
