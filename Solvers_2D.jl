module Solver2D
export ADI_solver, Pseudospec_FFT
# const m, ħ = 1, 1
Threads.nthreads()
    using LinearAlgebra, FFTW

    Base.@kwdef mutable struct Solver_Params
        par = init_params
        Nx::Int32 = par.Nx
        Ny::Int32 = par.Ny
        V::Array{Float64,2} = par.potential_field

        ############### FDM ADI METHOD ##################################
        λ::Float64 = (2par.dx^2) / par.dt
        A::Array{ComplexF64,2} = diff_matrix(1+2im/λ, -im/λ, Nx, Ny)
        B::Array{ComplexF64,2} = diff_matrix(1-2im/λ, im/λ, Nx, Ny)
        C::Array{ComplexF64,2} = B / A
        U_a::Array{ComplexF64,2}  = exp.(im*par.dt*V/2)
        U_b::Array{ComplexF64,2}  = exp.(-im*par.dt*V/2)
        ################ PS FFT METHOD ###################################
        k = [ collect(0:((par.Nx-1)/2)) ; collect(-(par.Nx)/2:-1) ] *2π/(par.Nx*par.dx)
        FFT_a = exp.(-im.*k.^2*par.dt/4) 
        FFT_U = exp.(-im*V*par.dt/2)
    end

    # Constructs finite difference matrix ≡ 1D discrete laplacian matrix
    function diff_matrix(diag,offdiag, Nx, Ny)
        A = zeros(ComplexF64, Nx, Ny)
        A[diagind(A,0)] .= diag
        A[diagind(A,1)] .= A[diagind(A,-1)] .= offdiag
        return A
    end

    function ADI_solver(ψ, Solv_par::Solver_Params)
        #Threads.@threads for i in 2:par.Nx-1
        for i in 2:Solv_par.Nx-1
            tmp_ψ = ψ[ : , i ]
            ψ[ : , i ] = Solv_par.C * tmp_ψ
        end
        #Threads.@threads for j in 2:par.Ny-1
        for j in 2:Solv_par.Ny-1
            tmp_ψ = ψ[ j , : ]
            ψ[ j , : ] = Solv_par.C * tmp_ψ
            ψ[ j , : ] = (ψ[ j , : ]) .* (Solv_par.U_b[ j , : ] ./ Solv_par.U_a[ j , : ])
        end
        return ψ
    end

    function Pseudospec_FFT(ψ, Solv_par::Solver_Params)
        
        for i in 2:Solv_par.Nx-1
            tmp_ψ = fft(ψ[ : , i ])
            ψ[ : , i ] =Solv_par.FFT_a.*tmp_ψ
        end

        for j in 2:Solv_par.Ny-1
            tmp_ψ = fft(ψ[ j , : ])
            ψ[ j , : ] = Solv_par.FFT_a .* tmp_ψ

        end
        ψ = Solv_par.FFT_U .* ifft(ψ)
        ψ = Solv_par.FFT_a .* fft(ψ)
        return ifft(ψ)

        # for i in 2:Solv_par.Nx-1
        #     tmp_ψ = fft(ψ[ : , i ])
        #     ψ[ : , i ] =Solv_par.FFT_a.*tmp_ψ
        # end

        # for j in 2:Solv_par.Ny-1
        #     tmp_ψ = fft(ψ[ j , : ])
        #     ψ[ j , : ] = Solv_par.FFT_a .* tmp_ψ
        #     ψ[ j , : ] = (ψ[ j , : ]) .* (Solv_par.U_b[ j , : ] ./ Solv_par.U_a[ j , : ])
        # end
        # return ifft(ψ)

    end

end
