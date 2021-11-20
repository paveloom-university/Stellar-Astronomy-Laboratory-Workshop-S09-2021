# Notices

Git mirrors:
- [Codeberg](https://codeberg.org/paveloom-university/Stellar-Astronomy-Laboratory-Workshop-S09-2021)
- [GitHub](https://github.com/paveloom-university/Stellar-Astronomy-Laboratory-Workshop-S09-2021)
- [GitLab](https://gitlab.com/paveloom-g/university/s09-2021/stellar-astronomy-laboratory-workshop)
- [SourceHut](https://sr.ht/~paveloom/Stellar-Astronomy-Laboratory-Workshop-S09-2021/)

Reference mirrors:
- [GitHub Pages](https://paveloom-university.github.io/Stellar-Astronomy-Laboratory-Workshop-S09-2021)
- [GitLab Pages](https://paveloom-g.gitlab.io/university/s09-2021/stellar-astronomy-laboratory-workshop)

To build this crate's documentation with [KaTeX](https://katex.org/) support, run:

```bash
cargo doc
RUSTDOCFLAGS="--html-in-header assets/katex-header.html" cargo doc --no-deps --open
```

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

This project also provides Pluto notebooks. You can interact with them in the web interface:

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
julia --project=. notebooks/notebook.jl
```
