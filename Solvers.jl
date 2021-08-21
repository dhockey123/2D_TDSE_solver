module Solver2D
export ADI_solver

Threads.nthreads()
    using LinearAlgebra
    Base.@kwdef mutable struct Solver_Params
        par = init_params
        Nx::Int32 = par.Nx
        Ny::Int32 = par.Ny
        λ::Float64 = (2par.dx^2) / par.dt
        A::Array{ComplexF64,2} = diff_matrix(1+2im/λ, -im/λ, Nx, Ny)
        B::Array{ComplexF64,2} = diff_matrix(1-2im/λ, im/λ, Nx, Ny)
        C::Array{ComplexF64,2} = B / A
    end

    function diff_matrix(diag,offdiag, Nx, Ny)
        A = zeros(ComplexF64, Nx, Ny)
        A[diagind(A,0)] .= diag
        A[diagind(A,1)] .= A[diagind(A,-1)] .= offdiag
        return A
    end

    function ADI_solver(ψ, par::Solver_Params)
        #Threads.@threads for i in 2:par.Nx-1
        for i in 2:par.Nx-1
            tmp_ψ = ψ[ : , i ]
            ψ[ : , i ] = par.C * tmp_ψ
        end
        #Threads.@threads for j in 2:par.Ny-1
        for j in 2:par.Ny-1
            tmp_ψ = ψ[ j , : ]
            ψ[ j , : ] = par.C * tmp_ψ
        end
        return ψ
    end

end



# function ADI_solver()

# end