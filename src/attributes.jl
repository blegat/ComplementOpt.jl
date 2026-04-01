"""
    AbstractComplementarityRelaxation

Abstract type to implement any complementarity function ``\\psi``.

"""
abstract type AbstractComplementarityRelaxation end

"""
    DefaultComplementarityReformulation <: MOI.AbstractOptimizerAttribute

Optimizer attribute that sets the default [`AbstractComplementarityRelaxation`](@ref)
used to reformulate all complementarity constraints.

This default is used for any constraint that does not have a constraint-specific
reformulation set via [`ComplementarityReformulation`](@ref).

## Example

```julia
MOI.set(model, ComplementOpt.DefaultComplementarityReformulation(), ComplementOpt.ScholtesRelaxation(0.0))
```
"""
struct DefaultComplementarityReformulation <: MOI.AbstractOptimizerAttribute end

"""
    ComplementarityReformulation <: MOI.AbstractConstraintAttribute

Constraint attribute that overrides the [`AbstractComplementarityRelaxation`](@ref)
for a specific complementarity constraint.

When set, this takes precedence over the model-wide default set via
[`DefaultComplementarityReformulation`](@ref). When not set, [`MOI.get`](@ref) returns the
model-wide default.

## Example

```julia
MOI.set(model, ComplementOpt.DefaultComplementarityReformulation(), ComplementOpt.ScholtesRelaxation(0.0))
c = @constraint(model, x ⟂ y)
MOI.set(model, ComplementOpt.ComplementarityReformulation(), c, ComplementOpt.FischerBurmeisterRelaxation(1e-8))
```
"""
struct ComplementarityReformulation <: MOI.AbstractConstraintAttribute end
