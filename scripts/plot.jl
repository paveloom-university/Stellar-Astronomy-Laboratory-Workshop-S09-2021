# This is a Julia script to
# plot the integrated orbits

println('\n', " "^4, "> Loading the packages...")

using LaTeXStrings
using Plots

# Use the GR backend for plots
gr()

# Change some of the default parameters for plots
default(fontfamily="Computer Modern", dpi=300, legend=:topright)

# Define the number of iterations used in the input data
n = 5000000

# Define the step used in the input data
h = -0.001

# For each path in the output directory
for path in readdir("$(@__DIR__)/../data/output"; join=true)
    # Check if the path is a directory
    if isdir(path)
        # Define the paths to the binary files
        r_path = joinpath(path, "r.bin")
        z_path = joinpath(path, "z.bin")
        e_path = joinpath(path, "e.bin")
        # Check if the data files exist
        if isfile(r_path) && isfile(z_path) && isfile(e_path)
            # Get the ID of the object
            name = basename(path)

            # Define the output directory for plots
            output_dir = joinpath("$(@__DIR__)/../plots/$(name)")

            # Create a directory for the object in the plots folder
            if !isdir(output_dir)
                mkdir(output_dir)
            end

            # Prepare the arrays
            r = Vector{Float64}(undef, n + 1)
            z = Vector{Float64}(undef, n + 1)
            e = Vector{Float64}(undef, n + 1)

            # Read the data
            read!(r_path, r)
            read!(z_path, z)
            read!(e_path, e)

            println(" "^4, "> Plotting the data from \"$(name)\"...")

            # Plot the orbit on the R x Z plane
            p = plot(r[1:100:end], z[1:100:end]; label="", title=name, xlabel="R", ylabel="Z");

            # Point out the starting position
            scatter!(p, [r[1],], [z[1],]; label="");

            # Save the figure as PDF and PNG
            savefig(p, joinpath(output_dir, "$(name) (Orbit).pdf"))
            savefig(p, joinpath(output_dir, "$(name) (Orbit).png"))

            # Plot the total energy over time
            p = plot(
                0.0:(h * 100 / 1000):(h * (length(e) - 1) / 1000),
                e[1:100:end];
                label="",
                title=name,
                xlabel=L"T, Gyr",
                ylabel=L"E, km^2/s^2"
            );

            # Point out the initial value
            scatter!(p, [0.0,], [e[1],]; label="");

            # Save the figure as PDF and PNG
            savefig(p, joinpath(output_dir, "$(name) (Total energy).pdf"))
            savefig(p, joinpath(output_dir, "$(name) (Total energy).png"))
        end
    end
end

println()
