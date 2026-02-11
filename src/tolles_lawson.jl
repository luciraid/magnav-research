using HDF5
using LinearAlgebra
using Statistics
using Plots
using MagNav

# ---------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------

function load_flight_data(filepath)

    h5 = h5open(filepath, "r")

    Bx = read(h5["flux_b_x"])
    By = read(h5["flux_b_y"])
    Bz = read(h5["flux_b_z"])

    B_igrf = read(h5["mag_1_igrf"])
    t  = read(h5["tt"])

    close(h5)

    return Bx, By, Bz, B_igrf, t
end

# ---------------------------------------------------------
# 2. NUMERICAL DERIVATIVE
# ---------------------------------------------------------

function derivative(x, t)

    dx = similar(x)
    dx[1] = 0.0

    for i in 2:length(x)
        dt = t[i] - t[i-1]
        dx[i] = (x[i] - x[i-1]) / dt
    end

    return dx
end

# ---------------------------------------------------------
# 3. BUILD TL MATRIX
# ---------------------------------------------------------

function build_tl_matrix(Bx, By, Bz, dBx, dBy, dBz)

    N = length(Bx)
    A = zeros(N, 15)

    for i in 1:N

        A[i,1:3] = [1, 1, 1]

        A[i,4:6] = [Bx[i], By[i], Bz[i]]

        A[i,7:9] = [Bx[i]*By[i],
                    Bx[i]*Bz[i],
                    By[i]*Bz[i]]

        A[i,10:12] = [dBx[i], dBy[i], dBz[i]]

        A[i,13:15] = [dBx[i]*By[i],
                      dBx[i]*Bz[i],
                      dBy[i]*Bz[i]]
    end

    return A
end

# ---------------------------------------------------------
# 4. MAIN PIPELINE
# ---------------------------------------------------------

function run_tl(filepath)

    println("Loading data...")
    Bx, By, Bz, B_igrf, t = load_flight_data(filepath)

    println("Computing total field...")
    B_total = sqrt.(Bx.^2 .+ By.^2 .+ Bz.^2)

    println("Computing derivatives...")
    dBx = derivative(Bx, t)
    dBy = derivative(By, t)
    dBz = derivative(Bz, t)

    println("Building TL matrix...")
    A = build_tl_matrix(Bx, By, Bz, dBx, dBy, dBz)

    println("Forming disturbance...")
    disturbance = B_total .- B_igrf

    println("Removing invalid rows...")

    valid = .!(isnan.(disturbance) .|
               isnan.(sum(A, dims=2)[:]) .|
               isinf.(disturbance))

    A_clean = A[valid, :]
    disturbance_clean = disturbance[valid]
    B_total_clean = B_total[valid]
    t_clean = t[valid]

    println("Solving TL regression...")
    c = qr(A_clean) \ disturbance_clean

    println("Applying compensation...")
    correction = A_clean * c
    B_total_comp = B_total_clean .- correction

    println("Variance Reduction:")
    println("Total field: ", var(B_total_clean), " → ", var(B_total_comp))

    plot(t_clean, B_total_clean, label="Raw Total Field")
    plot!(t_clean, B_total_comp, label="TL Compensated Total Field")
    savefig("figures/Flt1002_TL_compensation.png")

    println("Saved plot to figures/Flt1002_TL_compensation.png")
end

# ---------------------------------------------------------
# RUN
# ---------------------------------------------------------

dataset_path = sgl_2020_train()
flight_path = joinpath(dataset_path, "Flt1002_train.h5")

println("Using file: ", flight_path)

run_tl(flight_path)
