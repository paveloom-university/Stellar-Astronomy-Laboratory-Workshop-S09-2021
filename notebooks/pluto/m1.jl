### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 2f0ce4ef-4c45-4175-a9a7-df70cb268629
begin
    CURRENT_DIR = @__DIR__
    ROOT_DIR = if basename(CURRENT_DIR) == "pluto"
        dirname(dirname(CURRENT_DIR))
    else
        CURRENT_DIR
    end

    using Pkg
    Pkg.activate(ROOT_DIR; io = devnull)

    using Logging
    disable_logging(Logging.Warn)

    println('\n', " "^4, "> Loading the packages...")

    using LaTeXStrings
    using Plots
    using PlutoUI

    # Use the GR backend for plots
    gr()

    # Change some of the default parameters for plots
    default(fontfamily = "Computer Modern", dpi = 300, legend = nothing)

    # Define the paths to output directories
    PLOTS_DIR = joinpath(ROOT_DIR, "plots")
    NOTEBOOKS_PLOTS_DIR = joinpath(PLOTS_DIR, "notebooks")
    M1_DIR = joinpath(NOTEBOOKS_PLOTS_DIR, "m1")

    # Make sure the needed directories exist
    mkpath(M1_DIR)

    # An alias for a floating point type
    F = Float64

    md"This is a Pluto notebook to test the parameters of the first model."
end

# ╔═╡ 52342993-f23a-4a50-b3e0-e9e4b032f1cf
md"""
Masses are in units of ``M_0 = 2.325 \times 10^7 M_\odot``, other parameters are in ``\text{kpc}``.

Initial values are taken from Bajkova, Bobylev (2020), ranges are taken from Granados et al. (2021).
"""

# ╔═╡ 1fa86bf7-e1dc-4cd1-9ff0-73c14225ef2c
begin
    # Bulge
    M_B_S = @bind M_B NumberField(
        43.0:0.001:2580.0,
        default = 443.0,
    )
    A_B_S = @bind A_B NumberField(0.0:0.001:0.05, default = 0.0)
    B_B_S = @bind B_B NumberField(0.0:0.001:2.0, default = 0.2672)

    # Disk
    M_D_S = @bind M_D NumberField(
        1500.0:0.001:12903.0,
        default = 2798.0,
    )
    A_D_S = @bind A_D NumberField(1.0:0.001:20.0, default = 4.40)
    B_D_S = @bind B_D NumberField(0.1:0.001:16.0, default = 0.3084)

    # Halo
    M_H_S = @bind M_H NumberField(
        43.0:0.001:12903.0,
        default = 12474.0,
    )
    A_H_S = @bind A_H NumberField(0.1:0.001:30.0, default = 7.7)

    md"""
    ``M_b``: $(M_B_S) ``a_b``: $(A_B_S) ``b_b``: $(B_B_S) \
    ``M_d``: $(M_D_S) ``a_d``: $(A_D_S) ``b_d``: $(B_D_S) \
    ``M_h``: $(M_H_S) ``a_h``: $(A_H_S)
    """
end

# ╔═╡ a5d896ff-570f-49a9-adf0-19f11cb899e0
begin
    "Calculate the value of the R derivative of the Miyamoto & Nagai potential
    ``[ 100 \\, \\text{km} \\, \\text{s}^{-2} ]``"
    function miyamoto_nagai_phi_dr(r, z, m, a, b)::F
        return m * r / (r^2 + (a + √(z^2 + b^2))^2)^(1.5)
    end

    "Calculate the value of the R derivative of the Navarro-Frenk-White potential
    ``[ 100 \\, \\text{km} \\, \\text{s}^{-2} ]``"
    function navarro_frenk_white_phi_dr(r, z, m, a)::F
        sq_sum = r^2 + z^2
        k_1 = √(sq_sum) / a + 1.0
        return m * r * log(k_1) / sq_sum^(1.5) - m * r / a / sq_sum / k_1
    end

    "Calculate the value of circular velocity according to
    ``v_c^2(r) = R \\left. \\frac{\\partial \\Phi(R, Z)}{\\partial R} \\right|_{Z=0}``
    ``[ \\text{km}^2 \\, \\text{s}^{-2} ]``"
    circ_vel(r) = √(r * phi_dr(r, 0) * 100)

    "Calculate the value of ``\\partial \\Phi(R, Z) / \\partial R``
    ``[ 100 \\, \\text{km} \\, \\text{s}^{-2} ]``"
    function phi_dr(r, z)::F
        return miyamoto_nagai_phi_dr(r, z, M_B, A_B, B_B) +
               miyamoto_nagai_phi_dr(r, z, M_D, A_D, B_D) +
               navarro_frenk_white_phi_dr(r, z, M_H, A_H)
    end

    println(" "^4, "> Plotting the rotation curve...")

    # Plot the rotation curve
    p1 = plot(
        circ_vel, 0.2, 20;
        title = "Rotation curve",
        xlabel = L"R \;\, [\mathrm{kpc}]",
        ylabel = L"v_c \; [\mathrm{km^2 \, s^{-2}}]",
        xminorticks = 8,
        yminorticks = 8
    )

    # Save the figure
    savefig(p1, joinpath(M1_DIR, "Rotation curve.pdf"))
    savefig(p1, joinpath(M1_DIR, "Rotation curve.png"))

    p1
end

# ╔═╡ d4ba276c-a578-4a64-9f54-c33c81225a0b
begin
    # Plot the rotation curve
    p2 = plot(
        circ_vel, 0.2, 200;
        title = "Rotation curve",
        xlabel = L"R \;\, [\mathrm{kpc}]",
        ylabel = L"v_c \; [\mathrm{km^2 \, s^{-2}}]",
        xticks = [1e-1, 1, 1e1, 1e2],
        xlims = (0.1, 390),
        xscale = :log10,
        xminorticks = 9,
        yminorticks = 8
    )

    # Save the figure
    savefig(p2, joinpath(M1_DIR, "Rotation curve (log).pdf"))
    savefig(p2, joinpath(M1_DIR, "Rotation curve (log).png"))

    println()

    p2
end

# ╔═╡ bb727eda-9d7c-48e2-bac6-23e90aa6aab1
md"""
From Bajkova, Bobylev (2020):

$(PlutoUI.LocalResource(joinpath(ROOT_DIR, "assets", "Rotation curve of the NFWBB model (Bajkova, Bobylev, 2020).png")))
"""

# ╔═╡ Cell order:
# ╟─2f0ce4ef-4c45-4175-a9a7-df70cb268629
# ╟─52342993-f23a-4a50-b3e0-e9e4b032f1cf
# ╟─1fa86bf7-e1dc-4cd1-9ff0-73c14225ef2c
# ╟─a5d896ff-570f-49a9-adf0-19f11cb899e0
# ╟─d4ba276c-a578-4a64-9f54-c33c81225a0b
# ╟─bb727eda-9d7c-48e2-bac6-23e90aa6aab1
