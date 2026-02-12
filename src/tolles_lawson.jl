using HDF5
using LinearAlgebra
using Statistics
using Plots
using MagNav

# ---------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------

function load_flight_data(filepath)

    h5 = h5open(filepath, "r")

    Bx = read(h5["flux_b_x"])
    By = read(h5["flux_b_y"])
    Bz = read(h5["flux_b_z"])
    yaw = read(h5["ins_yaw"])
    yaw_rate = read(h5["yaw_rate"])
    t  = read(h5["tt"])

    close(h5)

    return Bx, By, Bz, yaw, yaw_rate, t
end

# ---------------------------------------------------------
# 9-TERM TL MATRIX
# ---------------------------------------------------------

function build_9term_matrix(l, m, n)

    N = length(l)
    A = zeros(N, 9)

    for i in 1:N
        A[i,1:3] = [1, 1, 1]
        A[i,4:6] = [l[i], m[i], n[i]]
        A[i,7:9] = [l[i]*m[i], l[i]*n[i], m[i]*n[i]]
    end

    return A
end

# ---------------------------------------------------------
# SEGMENTED TL
# ---------------------------------------------------------

function run_segmented_tl(filepath)

    println("Loading data...")
    Bx, By, Bz, yaw, yaw_rate, t = load_flight_data(filepath)

    valid = .!(isnan.(Bx) .| isnan.(By) .| isnan.(Bz) .|
               isinf.(Bx) .| isinf.(By) .| isinf.(Bz))

    Bx = Bx[valid]
    By = By[valid]
    Bz = Bz[valid]
    yaw = yaw[valid]
    yaw_rate = yaw_rate[valid]
    t  = t[valid]

    B_total = sqrt.(Bx.^2 .+ By.^2 .+ Bz.^2)

    mag_valid = B_total .> 1.0

    Bx = Bx[mag_valid]
    By = By[mag_valid]
    Bz = Bz[mag_valid]
    yaw = yaw[mag_valid]
    yaw_rate = yaw_rate[mag_valid]
    t  = t[mag_valid]
    B_total = B_total[mag_valid]

    l = Bx ./ B_total
    m = By ./ B_total
    n = Bz ./ B_total

    A_full = build_9term_matrix(l, m, n)

    # -----------------------------------------------------
    # TRAIN ONLY ON TURNS
    # -----------------------------------------------------

    turn_idx = abs.(yaw_rate) .> 1.0

    println("Training samples (turning only): ", sum(turn_idx))
    println("Total samples: ", length(turn_idx))

    A_train = A_full[turn_idx, :]

    Bx_train = (Bx .- mean(Bx))[turn_idx]
    By_train = (By .- mean(By))[turn_idx]
    Bz_train = (Bz .- mean(Bz))[turn_idx]

    println("Solving TL on turning segments...")

    cx = A_train \ Bx_train
    cy = A_train \ By_train
    cz = A_train \ Bz_train

    println("Applying compensation to full flight...")

    Bx_comp = Bx .- A_full * cx
    By_comp = By .- A_full * cy
    Bz_comp = Bz .- A_full * cz

    B_total_comp = sqrt.(Bx_comp.^2 .+ By_comp.^2 .+ Bz_comp.^2)

    # -----------------------------------------------------
    # VALIDATION
    # -----------------------------------------------------

    corr_before = cor(B_total, yaw)
    corr_after  = cor(B_total_comp, yaw)

    println("Heading Correlation:")
    println("Before TL: ", corr_before)
    println("After TL : ", corr_after)

    println("Variance Reduction:")
    println(var(B_total), " → ", var(B_total_comp))

    plot(t, B_total, label="Raw")
    plot!(t, B_total_comp, label="Segmented TL")
    savefig("figures/Flt1002_segmented_TL.png")

    println("Saved segmented TL plot.")
end

# ---------------------------------------------------------
# RUN
# ---------------------------------------------------------

dataset_path = sgl_2020_train()
flight_path = joinpath(dataset_path, "Flt1002_train.h5")

println("Using file: ", flight_path)

run_segmented_tl(flight_path)
