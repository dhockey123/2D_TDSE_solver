module Potentials2D
    export make_U_doubleslit

    function make_U_doubleslit(par, x_pos, slit_size, slit_width)

        V = zeros(par.Nx,par.Ny)
        x_pos = x_pos*par.Nx / par.x_max
        slit_size = (par.Nx-1)*slit_size/maximum(par.x)
        width = slit_width*(par.Nx-1)/par.y_max
        pos_main = Int32(round(x_pos))
        pos = Int32(floor(par.Nx*.5))+1
        println((slit_size)/2, "\t", Int32(round((slit_size)/2)) , "\t Values must be the same for slit_width accuracy- Ensure 1/slit_width = integer")
        gap = Int32(round((slit_size)/2))
        println(width, "\t", Int32(round(width)), "\t Values must be the same for slit_size accuracy- Ensure 1/slit_size = integer")
        width = Int32(round(width))
        for i in 1:width
            V[ : , pos_main+i] .= 1000
            V[pos+gap:pos+gap, pos_main+i] .= 0
            V[pos-gap:pos-gap, pos_main+i] .= 0
            # V[pos_main+i, : ] .= 1000
            # V[pos_main+i, pos+gap:pos+3gap] .= 0
            # V[pos_main+i, pos-3gap:pos-gap] .= 0
        end
        return V
    end
end


function make_U_doubleslit(par, x_pos, slit_size, slit_width)
    N=11
    slit_size=2
    width=1
    x_pos=3
    x_max = 22
    y_max=20
    slit_width=3
    V = zeros(N, N)
    x_pos = x_pos*N / x_max
    slit_size = (N-1)*slit_size/x_max
    width = slit_width*(N-1)/y_max
    pos_main = Int32(round(x_pos))
    pos = Int32(floor(N*.5))+1
    #println(pos, "\tMust be odd integer")
    println((slit_size)/2, "\t", Int32(round((slit_size)/2)) , "\t Values must be the same for slit_width accuracy- Ensure 1/slit_width = integer")
    gap = Int32(round((slit_size)/2))
    println(width, "\t", Int32(round(width)), "\t Values must be the same for slit_size accuracy- Ensure 1/slit_size = integer")
    width = Int32(round(width))
    #println(V[pos+gap:pos+3gap, pos_main+width])
    for i in 0:width
        V[ : , pos_main+i] .= 1000
        V[pos+gap:pos+3gap, pos_main+i] .= 0
        V[pos-3gap:pos-gap, pos_main+i] .= 0
        # V[pos_main+i, : ] .= 1000
        # V[pos_main+i, pos+gap:pos+3gap] .= 0
        # V[pos_main+i, pos-3gap:pos-gap] .= 0
    end
    V