# This is a Julia script to
# plot the integrated orbits

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
    # A time step
    if ARGS[i] == "-h"
        check_last(i)
        try
            global H = parse(Float64, ARGS[i+1])
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
if !@isdefined(N) || !@isdefined(H) || length(ARGS) == 0
    println("""
        Usage:
        julia --project=. scripts/orbits.jl -n <N> -h <H> <input> [--postfix <POSTFIX]>"""
    )
    exit(1)
end

# Define the input directory
INPUT_DIR = ARGS[end]

println('\n', " "^4, "> Loading the packages...")

using LaTeXStrings
using Plots

# Use the GR backend for plots
gr()

# Change some of the default parameters for plots
default(fontfamily = "Computer Modern", dpi = 300, legend = :topright)

# Define the paths to output directories
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
PLOTS_DIR = joinpath(ROOT_DIR, "plots")
ORBITS_DIR = joinpath(PLOTS_DIR, "orbits")
OUTPUT_DIR = joinpath(ORBITS_DIR, basename(INPUT_DIR))

# Make sure the needed directories exist
mkpath(OUTPUT_DIR)

# For each path in the input directory
for path in readdir(INPUT_DIR; join = true)
    # Check if the path is a directory
    if isdir(path)
        # Define the paths to the binary files
        r_path = joinpath(path, "r.bin")
        z_path = joinpath(path, "z.bin")
        x_path = joinpath(path, "x.bin")
        y_path = joinpath(path, "y.bin")
        e_path = joinpath(path, "e.bin")

        # Get the ID of the object
        name = basename(path)

        # Define the object directory
        object_dir = joinpath(OUTPUT_DIR, name)
        mkpath(object_dir)

        # Plot the orbit in the RZ plane if the corresponding data files exist
        if isfile(r_path) && isfile(z_path)
            # Prepare the arrays
            r = Vector{Float64}(undef, N + 1)
            z = Vector{Float64}(undef, N + 1)

            # Read the data
            read!(r_path, r)
            read!(z_path, z)

            println(" "^4, "> Plotting the orbit in the RZ plane from \"$(name)\"...")

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
            savefig(p, joinpath(object_dir, "$(name) (Orbit, RZ)$(POSTFIX).pdf"))
            savefig(p, joinpath(object_dir, "$(name) (Orbit, RZ)$(POSTFIX).png"))
        end

        # Plot the orbit in the XY plane if the corresponding data files exist
        if isfile(x_path) && isfile(y_path)
            # Prepare the arrays
            x = Vector{Float64}(undef, N + 1)
            y = Vector{Float64}(undef, N + 1)

            # Read the data
            read!(x_path, x)
            read!(y_path, y)

            println(" "^4, "> Plotting the orbit in the XY plane from \"$(name)\"...")

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
            savefig(p, joinpath(object_dir, "$(name) (Orbit, XY)$(POSTFIX).pdf"))
            savefig(p, joinpath(object_dir, "$(name) (Orbit, XY)$(POSTFIX).png"))
        end

        # Plot the total energy variation if the corresponding data file exists
        if isfile(e_path)
            # Prepare the arrays
            e = Vector{Float64}(undef, N + 1)

            # Read the data
            read!(e_path, e)

            println(" "^4, "> Plotting the total energy variation from \"$(name)\"...")

            # Plot the total energy variation
            p = plot(
                0.0:(H*100/1000):(H*(length(e)-1)/1000),
                e[1:100:end];
                label = "",
                title = name,
                xlabel = L"T \;\, [\mathrm{Gyr}]",
                ylabel = L"E \; [\mathrm{km^2 \, s^{-2}}]"
            )

            # Point out the initial value
            scatter!(p, [0.0,], [e[1],]; label = "")

            # Save the figure as PDF and PNG
            savefig(p, joinpath(object_dir, "$(name) (Total energy)$(POSTFIX).pdf"))
            savefig(p, joinpath(object_dir, "$(name) (Total energy)$(POSTFIX).png"))
        end
    end
end

println()
