# This is a Julia script to
# plot the integrated orbits

println('\n', " "^4, "> Loading the packages...")

using LaTeXStrings
using Plots
using Statistics

# Use the GR backend for plots
gr()

# Change some of the default parameters for plots
default(fontfamily = "Computer Modern", dpi = 300, legend = :topright)

# Define the number of integration iterations used
n = 100000

# Define the time step used
h = -0.01

# Define the number of simulation iterations
s = 200

# Define the number of minor ticks for heatmaps
minorticks = 8

# Define postfix for names of files
postfix = if length(ARGS) > 0
    " ($(ARGS[1]))"
else
    ""
end

# Define the paths to output directories
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
DATA_OUTPUT = joinpath(ROOT_DIR, "data", "output")
PLOTS_OUTPUT = joinpath(ROOT_DIR, "plots")

# Plot a heatmap from the coordinates data
function plot_heatmap(name, output_dir, x, y, xlabel, ylabel, plane)
    # Plot the orbit
    p = plot(
        x[1:100:(n+1)],
        y[1:100:(n+1)];
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
    for i = 1:(s+1)
        # Get a copy of matrix on the previous step
        mₚ .= m
        # For each pair of the coordinates
        for j = ((i-1)*n+i):(i*n+i)
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
        m / (s + 1);
        color = cgrad([RGBA{Float64}(0.0, 0.0, 0.0, 0.0); cgrad(:inferno)[2:end]]),
        xticks = [xticks(p)[1][1][begin] - h₁ * minorticks; xticks(p)[1][1]; xticks(p)[1][1][end] + h₁ * 8],
        yticks = [yticks(p)[1][1][begin] - h₂ * minorticks; yticks(p)[1][1]; yticks(p)[1][1][end] + h₂ * 8]
    )

    # Redraw the orbit
    plot!(
        p,
        x[1:100:(n+1)],
        y[1:100:(n+1)];
        label = "",
        color = palette(:default)[1]
    )

    # Point out the starting position
    scatter!(p, [x[1],], [y[1],]; label = "")

    # Save the figure as PDF and PNG
    savefig(p, joinpath(output_dir, "$(name) (Simulated orbits, $(plane))$(postfix).pdf"))
    savefig(p, joinpath(output_dir, "$(name) (Simulated orbits, $(plane))$(postfix).png"))
end

# For each path in the output directory
for path in readdir(DATA_OUTPUT; join = true)
    # Check if the path is a directory
    if isdir(path)
        # Define the paths to the binary files
        r_path = joinpath(path, "r.bin")
        z_path = joinpath(path, "z.bin")
        x_path = joinpath(path, "x.bin")
        y_path = joinpath(path, "y.bin")
        e_path = joinpath(path, "e.bin")
        apo_path = joinpath(path, "apo.bin")
        peri_path = joinpath(path, "peri.bin")
        # Check if the data files exist
        if isfile(r_path) && isfile(z_path) && isfile(apo_path) && isfile(peri_path)
            # Get the ID of the object
            name = basename(path)

            # Define the output directories for plots
            output_dir = joinpath(PLOTS_OUTPUT, name)
            simulations_dir = joinpath(output_dir, "simulations")

            # Create directories for the object and
            # its simulations in the plots folder
            # if they don't exist
            if !isdir(output_dir)
                mkdir(output_dir)
            end
            if !isdir(simulations_dir)
                mkdir(simulations_dir)
            end

            # Prepare the arrays
            r = Vector{Float64}(undef, (n + 1) * (s + 1))
            z = Vector{Float64}(undef, (n + 1) * (s + 1))
            x = Vector{Float64}(undef, (n + 1) * (s + 1))
            y = Vector{Float64}(undef, (n + 1) * (s + 1))
            e = Vector{Float64}(undef, (n + 1) * (s + 1))
            apo = Vector{Float64}(undef, (s + 1))
            peri = Vector{Float64}(undef, (s + 1))

            # Read the data
            read!(r_path, r)
            read!(z_path, z)
            read!(x_path, x)
            read!(y_path, y)
            read!(e_path, e)
            read!(apo_path, apo)
            read!(peri_path, peri)

            println(" "^4, "> Plotting the data from \"$(name)\"...")

            # Calculate the number of bins
            bins = ceil(Int, (maximum(apo) - minimum(apo)) / (3.49 * std(apo) / s^(1 / 3)) * 1.5)

            # Plot the histogram of the apocentric distances
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
            savefig(p, joinpath(output_dir, "$(name) (Apocentric distances)$(postfix).pdf"))
            savefig(p, joinpath(output_dir, "$(name) (Apocentric distances)$(postfix).png"))

            # Calculate the number of bins
            bins = ceil(Int, (maximum(peri) - minimum(peri)) / (3.49 * std(peri) / s^(1 / 3)) * 1.5)

            # Plot the histogram of the pericentric distances
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
            savefig(p, joinpath(output_dir, "$(name) (Pericentric distances)$(postfix).pdf"))
            savefig(p, joinpath(output_dir, "$(name) (Pericentric distances)$(postfix).png"))

            # Plot the heatmaps of simulated orbits in the RZ- and XY-planes
            plot_heatmap(name, output_dir, r, z, L"R \; [\mathrm{kpc}]", L"Z \; [\mathrm{kpc}]", "RZ")
            plot_heatmap(name, output_dir, x, y, L"X \; [\mathrm{kpc}]", L"Y \;\, [\mathrm{kpc}]", "XY")
        end
    end
end

println()
