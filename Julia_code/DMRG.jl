using ITensors
using PyCall
using NPZ

let
    num=121
    N_point=Float64[]
    E_point=Float64[]
    C_point=Float64[]

    N = 101
    sites = siteinds("S=1/2",N)
    hz_edg=5
    hz=5
    m=0.0
    a=1

    os_N = OpSum() # Number operator defined as -\sum ZZ/2
    for j=1:N-1
        os_N += -2,"Sz",j,"Sz",j+1
    end
    NN = MPO(os_N,sites)

    os_E = OpSum() # Energy 

    os_E += -2*hz,"Sz",1
    os_E += -2*hz,"Sz",N
    for j=2:N-1
        os_E += -2*hz,"Sz",j
    end
    for j=1:N-1
        os_E += 2*(-m*((-1)^j)),"Sz",j,"Sz",j+1
    end
    for j=1:N-2
        os_E += -2/a,"Sz",j,"Sx",j+1,"Sz",j+2
        os_E += 1/2/a,"Sx",j+1
    end
    OE = MPO(os_E,sites)

    os_C = OpSum() # Chiral operator

    for j=1:N-1
        os_C += 2*(-1*((-1)^j)),"Sz",j,"Sz",j+1
    end
    OC = MPO(os_C,sites)

    for k=1:num
        mu=-0.05+0.05*k # DMRG varying with chemical potential

        os = OpSum() # The hamiltonian including chemical potential
        os += -2*hz_edg,"Sz",1
        os += -2*hz_edg,"Sz",N
        for j=2:N-1
            os_E += -2*hz,"Sz",j
        end

        for j=1:N-1
            os += 2*(mu-m*((-1)^j)),"Sz",j,"Sz",j+1
        end
        for j=1:N-2
            os += -2/a,"Sz",j,"Sx",j+1,"Sz",j+2
            os += 1/2/a,"Sx",j+1
        end

        H = MPO(os,sites)

        nsweeps = 100 # number of sweeps is 5
        maxdim = [10,20,100,100,200,200,500,1000,1500] # gradually increase states kept
        cutoff = [1E-10]

        psi = randomMPS(sites,256)

        energy,psi = dmrg(H,psi; nsweeps, maxdim, cutoff)

        N_val=real(inner(psi',NN,psi))
        push!(N_point,N_val)

        E_val=real(inner(psi',OE,psi))
        push!(E_point,E_val)

        C_val=real(inner(psi',OC,psi))
        push!(C_point,C_val)
    end

    file_name="Data_E_h_"*string(hz)*"_m_"*string(m)*"_a_"*string(a)*"_L"*string(N)*"_no"
    npzwrite(file_name, E_point)
    file_name="Data_N_h_"*string(hz)*"_m_"*string(m)*"_a_"*string(a)*"_L"*string(N)*"_no"
    npzwrite(file_name, N_point)
    file_name="Data_C_h_"*string(hz)*"_m_"*string(m)*"_a_"*string(a)*"_L"*string(N)*"_no"
    npzwrite(file_name, C_point)
end

#Aver_h_corr=h_data_point 



