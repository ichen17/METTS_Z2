using ITensors
using PyCall
using NPZ

let
    function avg_err(v::Vector)
        N = length(v)
        avg = v[1] / N
        avg2 = v[1]^2 / N
        for j in 2:N
          avg += v[j] / N
          avg2 += v[j]^2 / N
        end
        return avg, √((avg2 - avg^2) / N)
    end
    
    
    #function main(; N=10, cutoff=1E-8, δτ=0.1, beta=2.0, NMETTS=3000, Nwarm=10)
    
    #mu_step=parse(Int64, ARGS[1])
    N=13
    cutoff=1E-9
    δτ=0.05
    beta=5.0
    NMETTS=20
    Nwarm=10

    Num_metts=100
    
    hz=0.1
    hz_edge=hz
    m=0
    mu=-0.4 #-0.88#+mu_step*0.05
    a=1
    # Make an array of 'site' indices
    s = siteinds("S=1/2", N)
    gates = ITensor[]
    
    for j in 1:(N)
        s1 = s[j]
        hi = op("Sz", s1)
        if j ==1 || j==N
            Gi = exp(hz_edge*δτ* hi)
        else
            Gi = exp(hz*δτ* hi)
        end
        push!(gates, Gi)
    end
    
    for j in 1:(N - 2)
        s1 = s[j]
        s2 = s[j + 1]
        s3 = s[j + 2]
        hi = op("Id", s1)*op("Sx", s2)* op("Id", s3)+ -4*op("Sz", s1) * op("Sx", s2) * op("Sz", s3)
    
        Gi = exp(-δτ/a/4 * hi)
    
        push!(gates, Gi)
    end
    
    for j in 1:(N - 1)
        s1 = s[j]
        s2 = s[j + 1]
    
        hj = op("Sz", s1) * op("Sz", s2)
        Gj = exp(-1*(mu-m*(-1)^j)*δτ* hj)
        push!(gates, Gj)
    end
    
    
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))
    
    # Make y-rotation gates to use in METTS collapses
    Rz_gates = ops([("Rz", n, (θ=-π / 2,)) for n in 1:N], s)
    #Ry_gates = ops([("Ry", n, (θ=-π / 2,)) for n in 1:N], s)
    #Ry_gates = ops([("Ry", n, (θ=-π / 2,)) for n in 1:N], s)
    #Rx_gates = ops([("Rx", n, (θ=π,)) for n in 1:N], s)
    #H_gates = ops([("H", n) for n in 1:N], s)
    H_gates = ITensor[op("H", s[j]) for j=1:N]
    
    # Arbitrary initial state
    psi = randomMPS(s)
    
    # Make H for measuring the energy
    os = OpSum()
    os += -2*hz_edge,"Sz",1
    for j=2:N-1
        os += -2*hz,"Sz",j
    end
    os += -2*hz_edge,"Sz",N
    
    for j=1:N-1
        os += 2*(mu-m*((-1)^j)),"Sz",j,"Sz",j+1
    end
    
    for j=1:N-2
        os += -2/a,"Sz",j,"Sx",j+1,"Sz",j+2
        os += 1/2/a,"Sx",j+1
    end
    
    H = MPO(os,s)
    #H = MPO(terms, s)
    
    os_N = OpSum()
    for j=1:N-1
        os_N += -2,"Sz",j,"Sz",j+1
    end
    NN = MPO(os_N,s)
    
    NNN=apply(NN, NN; cutoff)
    
    os_E = OpSum()
    
    for j=1:N
        os_E += -2*hz,"Sz",j
    end
    
    for j=1:N-1
        os_E += 2*(-m*((-1)^j)),"Sz",j,"Sz",j+1
    #os += 4*mu,"Sz",j,"Sz",j+1
    end
    
    for j=1:N-2
        os_E += -2/a,"Sz",j,"Sx",j+1,"Sz",j+2
        os_E += 1/2/a,"Sx",j+1
    end
    
    OE = MPO(os_E,s)
    
    os_C = OpSum()
    
    for j=1:N-1
        os_C += 2*(-1*((-1)^j)),"Sz",j,"Sz",j+1
    end
    
    OC = MPO(os_C,s)
    
    # Make τ_range and check δτ is commensurate
    τ_range = δτ:δτ:(beta / 2)
    if norm(length(τ_range) * δτ - beta / 2) > 1E-10
        error("Time step δτ=$δτ not commensurate with beta/2=$(beta/2)")
    end

    N_i_dict=zeros( N-1, Num_metts, NMETTS)

    #energies = Dict(string(1)=>Float64[])#Float64[]

    S_set = Set()
    for times in 1:Num_metts
        for step in 1:(Nwarm + NMETTS)
            if step <= Nwarm
                println("Making warmup METTS number $step")
            else
                println("Making METTS number $(step-Nwarm)")
            end
        
        # Do the time evolution by applying the gates
            for τ in τ_range
                psi = apply(gates, psi; cutoff)
                normalize!(psi)
            end
        
        # Measure properties after >= Nwarm 
        # METTS have been made
            if step > Nwarm

                energy = inner(psi', OE, psi)
                Nums = inner(psi', NN, psi)
                
                sum_sz=0
                for j in 1:N-1
                    os_N = OpSum()
                    os_N += -2,"Sz",j,"Sz",j+1
                    n_i=MPO(os_N,s)
                    N_i_val=real(inner(psi',n_i,psi))
                    sum_sz+=N_i_val

                    N_i_dict[j,times,step-Nwarm]=N_i_val
                end

            end
        
            if step % 2 == 1
                
                psi = apply(Rz_gates, psi)
                psi = apply(H_gates, psi)

                samp = sample!(psi)
                new_state = [samp[j] == 1 ? "Y+" : "Y-" for j in 1:N]
            else
                samp = sample!(psi)
                new_state = [samp[j] == 1 ? "Z+" : "Z-" for j in 1:N]
            end
            psi = productMPS(ComplexF64, s, new_state)
            push!(S_set,samp)
        end
    end
    mu=round(mu, digits=3)
    file_name="Data_Ni_h_"*string(hz)*"_a_"*string(a)*"_L"*string(N)*"_mu_"*string(mu)*"_500_T_"*string(beta)*"_S_"*string(NMETTS)*"_Y_Z_1.npz"
    npzwrite(file_name, N_i_dict)

end

#Aver_h_corr=h_data_point 



