module Potentials2D
    export make_U_doubleslit

    # Constructs the double slit potential field. 
    # par : initial params in main.jl, required for domain values (Nx, Ny, x/y_min, x/y_max)
    # barrier_pos : user defined position barrier is centered at (in graphs coordinates x, y)
    # barrier_width : width of barrier (in graphs coordinates x, y)
    # slit_gap : size of slits (in graphs coordinates x, y)
    
    function doubleslit_potential(par, barrier_pos, barrier_width , slit_gap)
        V = par.potential_field
        x_total_length = abs(par.x_max-par.x_min)
        y_total_length = abs(par.y_max-par.y_min)

        # Gets the y_midpoint and converts to domain coordinates
        # Used to define origin point for constructing slit placement
        y_mid = abs(par.y_max-par.y_min)/2
        y_mid = Int32(round( (par.Ny+1)*(y_mid/y_total_length) ))
        
        # Converts barrier graph coordinates (x, y) to domain coordinates (Nx, Ny)
        barrier_pos   = Int32(round( par.Nx*(abs(barrier_pos-par.x_min)/x_total_length ) ))
        barrier_width = Int32(round( par.Nx*barrier_width/x_total_length ))
        slit_gap      = Int32(round( ((par.Ny-2)*( slit_gap/y_total_length )) /2 ))
        
        # Fills in potential field at given domain coordinates
        V[ : ,                             barrier_pos:barrier_pos+barrier_width] .= 1000
        V[y_mid-3slit_gap:y_mid-slit_gap, barrier_pos:barrier_pos+barrier_width] .= 0
        V[y_mid+slit_gap:y_mid+3slit_gap, barrier_pos:barrier_pos+barrier_width] .= 0

        return transpose(V)
    end

    function oscillator_2D(par, mag)
        V  = par.potential_field
        V .= [0.5*mag*(x^2 + y^2 ) for x in par.x, y in par.y] 
    end

    # Constructs the cylindrical potential. 
    # par : initial params in main.jl, required for domain values (Nx, Ny, x/y_min, x/y_max)
    # pos_x, pos_y : origin point of circle (in graphs coordinates x, y)
    # radius : radius of circle (in graphs coordinates x, y)
    # U_magnitude : magnitude of potential field

    function cylinder_potential(par, pos_y, pos_x, radius, U_magnitude)
        V = par.potential_field
        x_total_length = abs(par.x_max-par.x_min)
        y_total_length = abs(par.y_max-par.y_min)

        # Gets center of circle points (a,b) and radius (r) as domain coordinates
        a = par.Nx * abs(pos_x - par.x_min) / x_total_length
        b = par.Ny * abs(pos_y - par.y_min) / y_total_length
        r = radius * (par.Nx / x_total_length)

        # Goes through each coordinate in 2D array and 
        # fills in the potential magnitude when the following circle formula is satisfied.. 
        for i in 1:par.Nx
            for j in 1:par.Ny
                if abs.((i-a)^2 + (j-b)^2) < r^2
                    V[j,i] = U_magnitude
                end
            end
        end
        return V
    end
end


