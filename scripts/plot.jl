# This is a Julia script to
# plot the integrated orbits

println('\n', " "^4, "> Loading the packages...")

using CSV
using DataFrames
using Plots

# Use the GR backend for plots
gr()

# Change the default font for plots
default(fontfamily="Computer Modern", dpi=300, legend=:topright)

for file in readdir("$(@__DIR__)/../data/output"; join=true)
    name_ext = basename(file)
    name = splitext(name_ext)[1]
    println(" "^4, "> Plotting the data from \"$(name_ext)\"...")
    df = DataFrame(CSV.File(file))
    p = plot(df.r, df.z; label=name);
    savefig(p, "$(@__DIR__)/../plots/$(name).pdf")
end

println()
