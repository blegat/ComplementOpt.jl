# MathOptComplements.jl

[![Build Status](https://github.com/jump-dev/MathOptComplements.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/jump-dev/MathOptComplements.jl/actions?query=workflow%3ACI)
[![Codecov branch](https://codecov.io/gh/jump-dev/MathOptComplements.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jump-dev/MathOptComplements.jl/branch/main)

[MathOptComplements.jl](https://github.com/jump-dev/MathOptComplements.jl) is a
JuMP extension for reformulating complementarity constraints.

## License

`MathOptComplements.jl` is licensed under the [MIT License](https://github.com/jump-dev/MathOptComplements.jl/blob/main/LICENSE.md).

## Getting help

If you need help, please ask a question on the [JuMP community forum](https://jump.dev/forum).

If you have a reproducible example of a bug, please [open a GitHub issue](https://github.com/jump-dev/MathOptLazy.jl/issues/new).

## Installation

Install `MathOptComplements` using `Pkg.add`:

```julia
import Pkg
Pkg.add(; url = "https://github.com/jump-dev/MathOptComplements.jl")
```

## Use with JuMP

Use `MathOptComplements.jl` with JuMP as follows:
```julia
using JuMP
import Ipopt
import MathOptComplements
model = Model(() -> MathOptComplements.Optimizer(Ipopt.Optimizer()))
set_attribute(
    model,
    MathOptComplements.DefaultComplementarityReformulation(),
    MathOptComplements.ScholtesRelaxation(0.0),
)
@variable(model, z[1:2])
set_lower_bound(z[2], 0)
@objective(model, Min, (z[1] - 1)^2 + z[2]^2)
@constraint(model, z[2] - z[1] ⟂ z[2])
optimize!(model)
```

If you use Ipopt, we recommend setting the following options to improve the
performance:
```julia
set_attribute(model, "mu_strategy", "adaptive")
set_attribute(model, "bound_push", 1e-1)
set_attribute(model, "bound_relax_factor", 0.0)
```

## Supported reformulations

You can change the reformulation by using the optimizer attribute
`MathOptComplements.DefaultComplementarityReformulation`. The following values
are supported. Check their docstrings for details.

- `MathOptComplements.ScholtesRelaxation(tau)` (**default**)
- `MathOptComplements.FischerBurmeisterRelaxation(tau)`
- `MathOptComplements.LiuFukushimaRelaxation(tau)`
- `MathOptComplements.KanzowSchwarzRelaxation(tau)`

Most reformulations are not equivalent to the original problem, which is why
they are not activated by default. [This arXiv paper](https://arxiv.org/html/2312.11022v2)
has a recent benchmark comparing the different reformulations on [MacMPEC](https://www.mcs.anl.gov/~leyffer/macmpec/).

## Funding
We acknowledge support from the [Fondation Mathématiques Jacques Hadamard](https://www.fondation-hadamard.fr/fr/)
which has funded the PGMO-IROE project "A new optimization suite for large-scale
market equilibrium".
