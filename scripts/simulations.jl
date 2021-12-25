# This is a Julia script to
# plot the simulated orbits

"Check if the value of the option is the last argument"
function check_last(i)
    if i + 1 == length(ARGS)
        println("The last argument should be the input directory.")
        exit(1)
    end
end

# Define default values for optional arguments
POSTFIX = ""

# Parse the options
for i in eachindex(ARGS)
    # A number of integration iterations used
    if ARGS[i] == "-n"
        check_last(i)
        try
            global N = parse(UInt, ARGS[i+1])
        catch
            println("Couldn't parse the value of the `-n` argument.")
            exit(1)
        end
    end
    # A number of simulation iterations
    if ARGS[i] == "-s"
        check_last(i)
        try
            global S = parse(UInt, ARGS[i+1])
        catch
            println("Couldn't parse the value of the `-h` argument.")
            exit(1)
        end
    end
    # A postfix for the names of output files
    if ARGS[i] == "--postfix"
        check_last(i)
        try
            global POSTFIX = " ($(ARGS[i+1]))"
        catch
            println("Couldn't parse the value of the `--postfix` argument.")
            exit(1)
        end
    end
end

# Check for required arguments (exc)
if !@isdefined(N) || !@isdefined(S) || length(ARGS) == 0
    println("""
        Usage:
        julia --project=. scripts/simulations.jl -n <N> -s <S> <input> [--postfix <POSTFIX>]"""
    )
    exit(1)
end

# Define the input directory
INPUT_DIR = ARGS[end]

println('\n', " "^4, "> Loading the packages...")

using LaTeXStrings
using Plots
using Statistics

# Use the GR backend for plots
gr()

# Change some of the default parameters for plots
default(fontfamily = "Computer Modern", dpi = 300, legend = :topright)

# Define the number of minor ticks for heatmaps
minorticks = 8

# Define the paths to output directories
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
PLOTS_DIR = joinpath(ROOT_DIR, "plots")
SIMULATIONS_DIR = joinpath(PLOTS_DIR, "simulations")
OUTPUT_DIR = joinpath(SIMULATIONS_DIR, basename(INPUT_DIR))

# Make sure the needed directories exist
mkpath(OUTPUT_DIR)

"Plot a heatmap from the coordinates data"
function plot_heatmap(name, object_dir, x, y, xlabel, ylabel, plane)
    # Plot the orbit
    p = plot(
        x[1:100:(N+1)],
        y[1:100:(N+1)];
        label = "",
        title = name,
        xlabel,
        ylabel,
        size = (400, 400),
        minorticks,
        minorgrid = true
    )

    # Define the ticks step
    h₁ = abs(xticks(p)[1][1][2] - xticks(p)[1][1][1]) / minorticks
    h₂ = abs(yticks(p)[1][1][2] - yticks(p)[1][1][1]) / minorticks

    # Define the ranges of the heatmap matrix
    r₁ = range(
        (xticks(p)[1][1][1] - h₁ * floor((xticks(p)[1][1][1] - xlims(p)[1]) / h₁)) - 5 * h₁,
        (xticks(p)[1][1][end] + h₁ * floor((xlims(p)[2] - xticks(p)[1][1][end]) / h₁)) + 5 * h₁;
        step = h₁
    )
    r₂ = range(
        (yticks(p)[1][1][1] - h₂ * floor((yticks(p)[1][1][1] - ylims(p)[1]) / h₂)) - 5 * h₂,
        (yticks(p)[1][1][end] + h₂ * floor((ylims(p)[2] - yticks(p)[1][1][end]) / h₂)) + 5 * h₂;
        step = h₂
    )

    # Get lengths of the ranges
    l₁ = length(r₁)
    l₂ = length(r₂)

    # Define the heatmap matrix
    m = Matrix{Int}(undef, l₂ - 1, l₁ - 1)
    mₚ = Matrix{Int}(undef, l₂ - 1, l₁ - 1)
    m .= 0

    # For each simulation iteration
    for i = 1:(S+1)
        # Get a copy of matrix on the previous step
        mₚ .= m
        # For each pair of the coordinates
        for j = ((i-1)*N+i):(i*N+i)
            # Check if the orbit is in the limits of the grid
            if (r₁[begin] ≤ x[j] ≤ r₁[end]) && (r₂[begin] ≤ y[j] ≤ r₂[end])
                # Get the indices of the square
                k = 1 + floor(Int, abs(x[j] - r₁[begin]) / h₁)
                l = 1 + floor(Int, abs(y[j] - r₂[begin]) / h₂)

                # Mark the crossed squares, once
                if m[l, k] == mₚ[l, k]
                    m[l, k] += 1
                end
            end
        end
    end

    # Plot the heatmap
    heatmap!(
        p,
        r₁,
        r₂,
        m / (S + 1);
        color = cgrad([RGBA{Float64}(0.0, 0.0, 0.0, 0.0); cgrad(:inferno)[2:end]]),
        xticks = [xticks(p)[1][1][begin] - h₁ * minorticks; xticks(p)[1][1]; xticks(p)[1][1][end] + h₁ * 8],
        yticks = [yticks(p)[1][1][begin] - h₂ * minorticks; yticks(p)[1][1]; yticks(p)[1][1][end] + h₂ * 8]
    )

    # Redraw the orbit
    plot!(
        p,
        x[1:100:(N+1)],
        y[1:100:(N+1)];
        label = "",
        color = palette(:default)[1]
    )

    # Point out the starting position
    scatter!(p, [x[1],], [y[1],]; label = "")

    # Save the figure as PDF and PNG
    savefig(p, joinpath(object_dir, "$(name) (Simulated orbits, $(plane))$(POSTFIX).pdf"))
    savefig(p, joinpath(object_dir, "$(name) (Simulated orbits, $(plane))$(POSTFIX).png"))
end

# For each path in the output directory
for path in readdir(INPUT_DIR; join = true)
    # Check if the path is a directory
    if isdir(path)
        # Define the paths to the binary files
        r_path = joinpath(path, "r.bin")
        z_path = joinpath(path, "z.bin")
        x_path = joinpath(path, "x.bin")
        y_path = joinpath(path, "y.bin")
        apo_path = joinpath(path, "apo.bin")
        peri_path = joinpath(path, "peri.bin")

        # Get the ID of the object
        name = basename(path)

        # Define the object directory
        object_dir = joinpath(OUTPUT_DIR, name)
        mkpath(object_dir)

        # Plot the histogram of apocentric distances if the corresponding data file exists
        if isfile(apo_path) && isfile(peri_path)
            # Prepare the array
            apo = Vector{Float64}(undef, (S + 1))

            # Read the data
            read!(apo_path, apo)

            println(
                " "^4, "> Plotting the histogram of apocentric distances from \"$(name)\"...", '\n',
                " "^4, "  > max: ", maximum(apo), '\n',
                " "^4, "  > min: ", minimum(apo), '\n',
                " "^4, "  > no var.: ", apo[1]
            )

            # Calculate the number of bins
            bins = ceil(Int, (maximum(apo) - minimum(apo)) / (3.49 * std(apo) / S^(1 / 3)) * 1.5)

            # Plot the histogram of apocentric distances
            p = histogram(
                apo;
                bins,
                label = "",
                title = name,
                xlabel = L"r_a \; [\mathrm{kpc}]",
                xminorticks = 10,
                color = "#80cdfd"
            )
            plot!(p, yticks = range(0, ceil(Int, ylims(p)[2]); step = 5))
            vline!(p, [apo[1],]; label = "", lw = 1.5)

            # Save the figure as PDF and PNG
            savefig(p, joinpath(object_dir, "$(name) (Apocentric distances)$(POSTFIX).pdf"))
            savefig(p, joinpath(object_dir, "$(name) (Apocentric distances)$(POSTFIX).png"))


        end

        # Plot the histogram of pericentric distances if the corresponding data file exists
        if isfile(peri_path)
            # Prepare the array
            peri = Vector{Float64}(undef, (S + 1))

            # Read the data
            read!(peri_path, peri)

            # Calculate the number of bins
            bins = ceil(Int, (maximum(peri) - minimum(peri)) / (3.49 * std(peri) / S^(1 / 3)) * 1.5)

            println(
                " "^4, "> Plotting the histogram of pericentric distances from \"$(name)\"...", '\n',
                " "^4, "  > max: ", maximum(peri), '\n',
                " "^4, "  > min: ", minimum(peri), '\n',
                " "^4, "  > no var.: ", peri[1]
            )

            # Plot the histogram of pericentric distances
            p = histogram(
                peri;
                bins,
                label = "",
                title = name,
                xlabel = L"r_p \; [\mathrm{kpc}]",
                xminorticks = 10,
                color = "#80cdfd"
            )
            plot!(p, yticks = range(0, ceil(Int, ylims(p)[2]); step = 5))
            vline!(p, [peri[1],]; label = "", lw = 1.5)

            # Save the figure as PDF and PNG
            savefig(p, joinpath(object_dir, "$(name) (Pericentric distances)$(POSTFIX).pdf"))
            savefig(p, joinpath(object_dir, "$(name) (Pericentric distances)$(POSTFIX).png"))
        end

        # Plot the heatmap of simulated orbits in the RZ plane if the corresponding data files exist
        if isfile(r_path) && isfile(z_path)
            # Prepare the arrays
            r = Vector{Float64}(undef, (N + 1) * (S + 1))
            z = Vector{Float64}(undef, (N + 1) * (S + 1))

            try
                # Read the data
                read!(r_path, r)
                read!(z_path, z)

                println(" "^4, "> Plotting the heatmap of simulated orbits in the RZ plane from \"$(name)\"...")
                plot_heatmap(name, object_dir, r, z, L"R \; [\mathrm{kpc}]", L"Z \; [\mathrm{kpc}]", "RZ")
            catch
            end
        end

        # Plot the heatmap of simulated orbits in the XY plane if the corresponding data files exist
        if isfile(x_path) && isfile(y_path)
            # Prepare the arrays
            x = Vector{Float64}(undef, (N + 1) * (S + 1))
            y = Vector{Float64}(undef, (N + 1) * (S + 1))

            try
                # Read the data
                read!(x_path, x)
                read!(y_path, y)

                println(" "^4, "> Plotting the heatmap of simulated orbits in the XY plane from \"$(name)\"...")
                plot_heatmap(name, object_dir, x, y, L"X \; [\mathrm{kpc}]", L"Y \;\, [\mathrm{kpc}]", "XY")
            catch
            end
        end
    end
end

println()
