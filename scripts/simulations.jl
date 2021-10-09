# This is a Julia script to
# plot the integrated orbits

println('\n', " "^4, "> Loading the packages...")

using LaTeXStrings
using Plots
using Statistics

# Use the GR backend for plots
gr()

# Change some of the default parameters for plots
default(fontfamily="Computer Modern", dpi=300, legend=:topright)

# Define the number of integration iterations used
n = 100000

# Define the time step used
h = -0.01

# Define the number of simulation iterations
s = 200

# Define postfix for names of files
postfix = if length(ARGS) > 0
    " ($(ARGS[1]))"
else
    ""
end

# For each path in the output directory
for path in readdir("$(@__DIR__)/../data/output"; join=true)
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
            output_dir = joinpath("$(@__DIR__)/../plots/$(name)")
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
                label="",
                title=name,
                xlabel=L"r_a \; [\mathrm{kpc}]",
                xminorticks=10,
                color="#80cdfd",
            );
            plot!(p, yticks=range(0, ceil(Int, ylims(p)[2]); step=5));
            vline!(p, [apo[1],]; label="", lw=1.5)

            # Save the figure as PDF and PNG
            savefig(p, joinpath(output_dir, "$(name) (Apocentric distances)$(postfix).pdf"))
            savefig(p, joinpath(output_dir, "$(name) (Apocentric distances)$(postfix).png"))

            # Calculate the number of bins
            bins = ceil(Int, (maximum(peri) - minimum(peri)) / (3.49 * std(peri) / s^(1 / 3)) * 1.5)

            # Plot the histogram of the pericentric distances
            p = histogram(
                peri;
                bins,
                label="",
                title=name,
                xlabel=L"r_p \; [\mathrm{kpc}]",
                xminorticks=10,
                color="#80cdfd",
            );
            plot!(p, yticks=range(0, ceil(Int, ylims(p)[2]); step=5));
            vline!(p, [peri[1],]; label="", lw=1.5)

            # Save the figure as PDF and PNG
            savefig(p, joinpath(output_dir, "$(name) (Pericentric distances)$(postfix).pdf"))
            savefig(p, joinpath(output_dir, "$(name) (Pericentric distances)$(postfix).png"))

            # For each simulation iteration
            for i in 1:(s + 1)
                # Define the output directory for this iteration
                s_output_dir = i == 1 ? output_dir : joinpath("$(simulations_dir)/$(i - 1)")
                    
                # Create a directory for this iteration if it doesn't exist
                if !isdir(s_output_dir)
                    mkdir(s_output_dir)
                end

                # Plot the orbit in the RZ-plane
                p = plot(
                    r[((i - 1) * n + i):100:(i * n + i)],
                    z[((i - 1) * n + i):100:(i * n + i)];
                    label="",
                    title=name,
                    xlabel=L"R \; [\mathrm{kpc}]",
                    ylabel=L"Z \; [\mathrm{kpc}]"
                );

                # Point out the starting position
                scatter!(p, [r[((i - 1) * n + i)],], [z[((i - 1) * n + i)],]; label="");

                # Save the figure as PDF and PNG
                savefig(p, joinpath(s_output_dir, "$(name) (Orbit, RZ)$(postfix).pdf"))
                savefig(p, joinpath(s_output_dir, "$(name) (Orbit, RZ)$(postfix).png"))

                # Plot the orbit in the XY-plane
                p = plot(
                    x[((i - 1) * n + i):100:(i * n + i)],
                    y[((i - 1) * n + i):100:(i * n + i)];
                    label="",
                    title=name,
                    xlabel=L"X \; [\mathrm{kpc}]",
                    ylabel=L"Y \;\, [\mathrm{kpc}]"
                );

                # Point out the starting position
                scatter!(p, [x[((i - 1) * n + i)],], [y[((i - 1) * n + i)],]; label="");

                # Save the figure as PDF and PNG
                savefig(p, joinpath(s_output_dir, "$(name) (Orbit, XY)$(postfix).pdf"))
                savefig(p, joinpath(s_output_dir, "$(name) (Orbit, XY)$(postfix).png"))

                # Plot the total energy over time
                p = plot(
                    0.0:(h * 100 / 1000):(h * n / 1000),
                    e[((i - 1) * n + i):100:(i * n + i)];
                    label="",
                    title=name,
                    xlabel=L"T \;\, [\mathrm{Gyr}]",
                    ylabel=L"E \; [\mathrm{km^2 \, s^{-2}}]"
                );

                # Point out the initial value
                scatter!(p, [0.,], [e[((i - 1) * n + i)],]; label="");

                # Save the figure as PDF and PNG
                savefig(p, joinpath(s_output_dir, "$(name) (Total energy)$(postfix).pdf"))
                savefig(p, joinpath(s_output_dir, "$(name) (Total energy)$(postfix).png"))
            end
        end
    end
end

println()
