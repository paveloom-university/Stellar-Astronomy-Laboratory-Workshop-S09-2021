## How to reproduce the results

1. Check out this repository:

- [Codeberg](https://codeberg.org/paveloom-university/Stellar-Astronomy-Laboratory-Workshop-S09-2021)
- [GitHub](https://github.com/paveloom-university/Stellar-Astronomy-Laboratory-Workshop-S09-2021)
- [GitLab](https://gitlab.com/paveloom-g/university/s09-2021/stellar-astronomy-laboratory-workshop)

2. Make sure you have the following installed:

- [Julia](https://julialang.org)
- [Rust](https://www.rust-lang.org)
- [Tectonic](https://tectonic-typesetting.github.io)

3. Prepare the input data for calculating orbits:

    ```bash
    sed -i '3,17s/# "/"/' data/input/initial.dat
    ```

    > ***HINT:*** Execute this and all of the following commands in the root directory of this repository.

4. Integrate the orbits for 5 Gyr backward using both models:

    4.1. Prepare the output directories:

    ```bash
    mkdir -p data/output/M1
    mkdir -p data/output/M2
    ```

    4.2. Build and run the program

    ```bash
    cargo run --release -- --model 1 -n 500000 -h -0.01 -f=x,y,r,z,e -o data/output/M1 data/input/initial.dat
    cargo run --release -- --model 2 -n 500000 -h -0.01 -f=x,y,r,z,e -o data/output/M2 data/input/initial.dat
    ```

5. Generate plots of the orbits in XY and RZ planes, plus total energy variations:

    5.1. Instantiate the Julia project:

    ```bash
    julia --project=. -e "using Pkg; Pkg.instantiate()"
    ```

    5.2. Run the script:

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

6. Prepare the input data for running simulations:

    ```bash
    sed -i '3,17s/^"/# "/' data/input/initial.dat
    sed -i 's/# "NGC 5927/"NGC 5927/g' data/input/initial.dat
    sed -i 's/# "Pyxis/"Pyxis/g' data/input/initial.dat
    ```

7. Run 200 simulations using the Monte Carlo method for 1 Gyr backward using both models:

    ```bash
    cargo run --release -- --model 1 -n 100000 -h -0.01 -s 200 --simulate -f=r,z,x,y,apo,peri -o data/output/M1 data/input/initial.dat
    cargo run --release -- --model 2 -n 100000 -h -0.01 -s 200 --simulate -f=r,z,x,y,apo,peri -o data/output/M2 data/input/initial.dat
    ```

    > ***NOTE:*** This will take some time and quite a bit of disk space.

8. Generate histograms of apocentric and pericentric distances and plots of the simulated orbits:

    ```bash
    julia --project=. scripts/simulations.jl -n 100000 -s 200 data/output/M1
    julia --project=. scripts/simulations.jl -n 100000 -s 200 data/output/M2
    ```

    *or*

    ```bash
    ./julia.bash scripts/simulations.jl -n 100000 -s 200 data/output/M1
    ./julia.bash scripts/simulations.jl -n 100000 -s 200 data/output/M2
    ```

    The latter will start a Julia [daemon](https://github.com/dmolina/DaemonMode.jl) in the background. To kill it, run

    ```bash
    ./julia.bash kill
    ```

    > ***HINT:*** See the results in `plots/simulations`.

9. Compile the report:

    ```bash
    tectonic -X compile report/report.tex
    ```

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
