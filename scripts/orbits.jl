# This is a Julia script to
# plot the integrated orbits

println('\n', " "^4, "> Loading the packages...")

using LaTeXStrings
using Plots

# Use the GR backend for plots
gr()

# Change some of the default parameters for plots
default(fontfamily = "Computer Modern", dpi = 300, legend = :topright)

# Define the number of integration iterations used
n = 100000

# Define the time step used
h = -0.01

# Define if reverse mode was used
rev = false

# Define postfix for names of files
postfix = if length(ARGS) > 0
    " ($(ARGS[1]))"
else
    ""
end

# Define the paths to output directories
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
DATA_OUTPUT_DIR = joinpath(ROOT_DIR, "data", "output")
PLOTS_DIR = joinpath(ROOT_DIR, "plots")
ORBITS_DIR = joinpath(PLOTS_DIR, "orbits")

# Make sure the needed directories exist
mkpath(ORBITS_DIR)

# For each path in the output directory
for path in readdir(DATA_OUTPUT_DIR; join = true)
    # Check if the path is a directory
    if isdir(path)
        # Define the paths to the binary files
        r_path = joinpath(path, "r.bin")
        z_path = joinpath(path, "z.bin")
        x_path = joinpath(path, "x.bin")
        y_path = joinpath(path, "y.bin")
        e_path = joinpath(path, "e.bin")
        # Check if the data files exist
        if isfile(r_path) && isfile(z_path) && isfile(x_path) && isfile(y_path) && isfile(e_path)
            # Get the ID of the object
            name = basename(path)

            # Define the output directory for plots
            output_dir = joinpath(ORBITS_DIR, name)
            mkpath(output_dir)

            # Prepare the arrays
            r = Vector{Float64}(undef, n + 1)
            z = Vector{Float64}(undef, n + 1)
            x = Vector{Float64}(undef, n + 1)
            y = Vector{Float64}(undef, n + 1)
            e = Vector{Float64}(undef, n + 1)

            # Read the data
            read!(r_path, r)
            read!(z_path, z)
            read!(x_path, x)
            read!(y_path, y)
            read!(e_path, e)

            println(" "^4, "> Plotting the data from \"$(name)\"...")

            # Plot the orbit in the RZ-plane
            p = plot(
                r[1:100:end],
                z[1:100:end];
                label = "",
                title = name,
                xlabel = L"R \; [\mathrm{kpc}]",
                ylabel = L"Z \; [\mathrm{kpc}]"
            )

            # Point out the starting position
            scatter!(p, [r[1],], [z[1],]; label = "")

            # Save the figure as PDF and PNG
            savefig(p, joinpath(output_dir, "$(name) (Orbit, RZ)$(postfix).pdf"))
            savefig(p, joinpath(output_dir, "$(name) (Orbit, RZ)$(postfix).png"))

            # Plot the orbit in the XY-plane
            p = plot(
                x[1:100:end],
                y[1:100:end];
                label = "",
                title = name,
                xlabel = L"X \; [\mathrm{kpc}]",
                ylabel = L"Y \;\, [\mathrm{kpc}]"
            )

            # Point out the starting position
            scatter!(p, [x[1],], [y[1],]; label = "")

            # Save the figure as PDF and PNG
            savefig(p, joinpath(output_dir, "$(name) (Orbit, XY)$(postfix).pdf"))
            savefig(p, joinpath(output_dir, "$(name) (Orbit, XY)$(postfix).png"))

            # Plot the total energy over time
            p = plot(
                0.0:(h*100/1000):(h*(length(e)-1)/1000),
                e[1:100:end];
                label = "",
                title = name,
                xlabel = L"T \;\, [\mathrm{Gyr}]",
                ylabel = L"E \; [\mathrm{km^2 \, s^{-2}}]"
            )

            # Point out the initial value
            scatter!(p, [0.0,], [e[1],]; label = "")

            # Save the figure as PDF and PNG
            savefig(p, joinpath(output_dir, "$(name) (Total energy)$(postfix).pdf"))
            savefig(p, joinpath(output_dir, "$(name) (Total energy)$(postfix).png"))
        end
    end
end

println()
