module SetGlobals

using CSV
using DataFrames
using Interpolations
using KernelDensity
using Statistics 


export G, M_bh, int_steps, coulomb_log, m_test, grid_size
export find_psi, find_Ftot, find_Ntot
export read_in, Ntot_grid, Ftot_grid, Ntot, Ftot, psi_itp
export m_tab, r_tab, vr_tab, vt_tab
#export id_tab, startype_tab, binflag
export vmin, vmax, rmin, rmax


#Constants
const G=1
const M_bh=8.500250825e-03
const int_steps=1000
const grid_size=0.001


# Psi via Henon
function find_psi(r_, m_)
    N = length(r_)
    
    log_r = log.(r_)
    r_inv = [1/r_[i] for i in 1:N]
    push!(r_inv, 0.0)
    
    # Starting arrays
    psi_ = zeros(N + 1)
    M_ = zeros(N + 1)
    psi_[N + 1] = 0.0 
    M_[N] = sum(m_)

    # Iteration
    for k in N:-1:1
        psi_[k] = psi_[k + 1] - G * M_[k] * (r_inv[k] - r_inv[k + 1]) - G * M_bh * (r_inv[k]) # Added potential from central BH
        if k > 1
            M_[k - 1] = M_[k] - m_[k]
        end
    end
    psi_array = psi_[1:N]

    return interpolate((log_r,), psi_array, Gridded(Linear())) # Interpolate in log space 
end

function find_Ftot(tab_r,tab_v,tab_m,rmin,rmax,vmin,vmax,dr,dv)

    Ni = Int(floor((rmax - rmin) / dr))
    Nj = Int(floor((vmax - vmin) / dv))

    m2_sum_grid = zeros(Float64, Ni, Nj)

    # Bin data into cells by (i,j)
    for (r, v, m) in zip(tab_r, tab_v, tab_m)
        i = Int(floor((r - rmin) / dr)) + 1
        j = Int(floor((v - vmin) / dv)) + 1
        if 1 <= i <= Ni && 1 <= j <= Nj
            m2_sum_grid[i, j] += m ^ 2
        end
    end

    # Compute Ftot
    Ftot_grid = m2_sum_grid ./ (dr * dv)

    # Using lower left cell corners
    r_grid = rmin .+ dr .* (0:(Ni - 1))
    v_grid = vmin .+ dv .* (0:(Nj - 1))

    # Bilinear interpolation
    itp_Ftot = interpolate((log.(r_grid),log.(v_grid)), Ftot_grid, Gridded(Linear()))

    # Double check bounds 
    function Ftot(r, v)
        if r < rmin || r > r_grid[end] || v < vmin || v > v_grid[end]
            return 0.0
        else
            return itp_Ftot(log(r), log(v))
        end
    end

    return Ftot_grid, Ftot
end

function find_Ntot(tab_r,tab_v,tab_m,rmin,rmax,vmin,vmax,dr,dv)

    Ni = Int(floor((rmax - rmin) / dr))
    Nj = Int(floor((vmax - vmin) / dv))

    m_sum_grid = zeros(Float64, Ni, Nj)

    # Bin data into cells by (i,j)
    for (r, v, m) in zip(tab_r, tab_v, tab_m)
        i = Int(floor((r - rmin) / dr)) + 1
        j = Int(floor((v - vmin) / dv)) + 1
        if 1 <= i <= Ni && 1 <= j <= Nj
            m_sum_grid[i, j] += m
        end
    end

    # Compute Ntot
    Ntot_grid = m_sum_grid ./ (dr * dv)

    # Using lower left cell corners
    r_grid = rmin .+ dr .* (0:(Ni - 1))
    v_grid = vmin .+ dv .* (0:(Nj - 1))

    # Bilinear interpolation
    itp_Ntot = interpolate((log.(r_grid), log.(v_grid)), Ntot_grid, Gridded(Linear()))

    # Double check bounds 
    function Ntot(r, v)
        if r < rmin || r > r_grid[end] || v < vmin || v > v_grid[end]
            return 0.0
        else
            return itp_Ntot(log(r), log(v))
        end
    end

    return Ntot_grid, Ntot
end

function read_in(filename)

    # read in data
    data=CSV.File(filename)
    omega_cen=DataFrame(data)

    # set up arrays
    global id_tab, m_tab, r_tab, vr_tab, vt_tab, startype_tab, binflag = eachcol(omega_cen[:,1:7])
    global v_tab = sqrt.(vr_tab .^ 2 + vt_tab .^2)
    

    # find bounds
    global vmin=minimum(v_tab)
    global vmax=maximum(v_tab)

    global rmin=minimum(r_tab)
    global rmax=maximum(r_tab)

    global coulomb_log = M_bh / mean(m_tab)
    global m_test = mean(m_tab) # Test mass

    #find df, potential 
    global psi_itp=find_psi(r_tab,m_tab)
    global Ntot_grid,Ntot=find_Ntot(r_tab,v_tab,m_tab,rmin,rmax,vmin,vmax,grid_size,grid_size)
    global Ftot_grid,Ftot=find_Ftot(r_tab,v_tab,m_tab,rmin,rmax,vmin,vmax,grid_size,grid_size)
end

end

