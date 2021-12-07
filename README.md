## How to reproduce the results

1. Check out this repository:

- [Codeberg](https://codeberg.org/paveloom-university/Stellar-Astronomy-Laboratory-Workshop-S09-2021)
- [GitHub](https://github.com/paveloom-university/Stellar-Astronomy-Laboratory-Workshop-S09-2021)
- [GitLab](https://gitlab.com/paveloom-g/university/s09-2021/stellar-astronomy-laboratory-workshop)
- [SourceHut](https://sr.ht/~paveloom/Stellar-Astronomy-Laboratory-Workshop-S09-2021/)

2. Integrate the orbits for 5 Gyr backward using both models:

    ```bash
    cargo run --release -- --model 1 -n 500000 -h -0.01 -f=x,y,r,z -o data/output/M1 data/input/initial.dat
    cargo run --release -- --model 2 -n 500000 -h -0.01 -f=x,y,r,z -o data/output/M2 data/input/initial.dat
    ```

3. Generate plots of the orbits in XY and RZ planes:

    3.1. Instantiate the Julia project:

    ```bash
    julia --project=. -e "using Pkg; Pkg.instantiate()"
    ```

    3.2. Run the script:

    ```bash
    julia --project=. scripts/orbits.jl -n 500000 -h -0.01 data/output/M1
    julia --project=. scripts/orbits.jl -n 500000 -h -0.01 data/output/M2
    ```

    *or*

    ```bash
    ./julia.bash scripts/orbits.jl -n 500000 -h -0.01 data/output/M1
    ./julia.bash scripts/orbits.jl -n 500000 -h -0.01 data/output/M2
    ```

    The latter will start a Julia [daemon](https://github.com/dmolina/DaemonMode.jl) in the background. To kill it, run

    ```bash
    ./julia.bash kill
    ```

    > ***HINT:*** See the results in `plots/orbits`.

## Notices

Reference mirrors:
- [GitHub Pages](https://paveloom-university.github.io/Stellar-Astronomy-Laboratory-Workshop-S09-2021)
- [GitLab Pages](https://paveloom-g.gitlab.io/university/s09-2021/stellar-astronomy-laboratory-workshop)

### KaTeX

To build this crate's documentation with [KaTeX](https://katex.org/) support, run:

```bash
cargo doc
RUSTDOCFLAGS="--html-in-header assets/katex-header.html" cargo doc --no-deps --open
```

### Julia

This project provides [Julia](https://julialang.org) scripts. Make sure to use the project files (`Project.toml`) when running them:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. scripts/script.jl
```

Alternatively, you can use the `julia.bash` script, which starts a [daemon](https://github.com/dmolina/DaemonMode.jl) and runs scripts through it:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
./julia.bash scripts/script.jl
```

To kill the daemon run

```bash
./julia.bash kill
```

### Pluto

This project provides Pluto notebooks. You can interact with them in the web interface:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=.
```

```julia
using Pluto
Pluto.run()
```

Alternatively, you can run them as scripts:

```
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. notebooks/pluto/notebook.jl
```

### wxMaxima

This project provides a [wxMaxima](https://wxmaxima-developers.github.io/wxmaxima/) notebook.
