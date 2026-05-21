# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module MathOptComplements

import MathOptInterface as MOI

"""
    AbstractComplementarityRelaxation

Abstract type to implement any complementarity function ``\\psi``.

"""
abstract type AbstractComplementarityRelaxation end

"""
    DefaultComplementarityReformulation <: MOI.AbstractOptimizerAttribute

Optimizer attribute that sets the default
[`AbstractComplementarityRelaxation`](@ref) used to reformulate all
complementarity constraints.

This default is used for any constraint that does not have a constraint-specific
reformulation set via [`ComplementarityReformulation`](@ref).

## Example

```julia
julia> using JuMP, MathOptComplements

julia> model = Model();

julia> set_attribute(
           model,
           MathOptComplements.DefaultComplementarityReformulation(),
           MathOptComplements.ScholtesRelaxation(0.0),
       )
```
"""
struct DefaultComplementarityReformulation <: MOI.AbstractOptimizerAttribute end

"""
    ComplementarityReformulation <: MOI.AbstractConstraintAttribute

Constraint attribute that overrides the
[`AbstractComplementarityRelaxation`](@ref) for a specific complementarity
    constraint.

When set, this takes precedence over the model-wide default set via
[`DefaultComplementarityReformulation`](@ref). When not set, [`MOI.get`](@ref)
returns the model-wide default.

## Example

```julia
julia> using JuMP, MathOptComplements

julia> model = Model();

julia> set_attribute(
           model,
           MathOptComplements.DefaultComplementarityReformulation(),
           MathOptComplements.ScholtesRelaxation(0.0),
       )

julia> c = @constraint(model, x ⟂ y);

julia> set_attribute(
           c,
           MathOptComplements.ComplementarityReformulation(),
           MathOptComplements.FischerBurmeisterRelaxation(1e-8),
       )
```
"""
struct ComplementarityReformulation <: MOI.AbstractConstraintAttribute end

"""
    ComplementsWithSetType{S<:MOI.AbstractSet} <: MOI.AbstractVectorSet

Complementarity constraint where each slack variable (second half of the vector)
is asserted to belong to set type `S`. For a constraint with `dimension = 2n`,
the `n` complementarity pairs are `(x[i], x[i+n])` for `i = 1, …, n`.

`S` can be a scalar set type (`MOI.GreaterThan{T}`, `MOI.LessThan{T}`,
`MOI.EqualTo{T}`, `MOI.Interval{T}`) or a vector set type (`MOI.Nonnegatives`,
`MOI.Nonpositives`, `MOI.Zeros`).
"""
struct ComplementsWithSetType{S<:MOI.AbstractSet} <: MOI.AbstractVectorSet
    dimension::Int
end

MOI.dimension(set::ComplementsWithSetType) = set.dimension

include("Bridges/Bridges.jl")

using .Bridges:
    ScholtesRelaxation,
    FischerBurmeisterRelaxation,
    LiuFukushimaRelaxation,
    KanzowSchwarzRelaxation

include("MOI_wrapper.jl")

end # module MathOptComplements
