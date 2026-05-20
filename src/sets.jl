# Copyright (c) 2025 François Pacaud, Benoît Legat, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    ComplementsWithSetType{S<:MOI.AbstractSet} <: MOI.AbstractVectorSet

Complementarity constraint where each slack variable (second half of the vector)
is asserted to belong to set type `S`. For a constraint with `dimension = 2n`,
the `n` complementarity pairs are `(x[i], x[i+n])` for `i = 1, …, n`.

`S` can be a scalar set type (`MOI.GreaterThan{T}`, `MOI.LessThan{T}`,
`MOI.EqualTo{T}`, `MOI.Interval{T}`) or a vector set type
(`MOI.Nonnegatives`, `MOI.Nonpositives`, `MOI.Zeros`).

"""
struct ComplementsWithSetType{S<:MOI.AbstractSet} <: MOI.AbstractVectorSet
    dimension::Int
end

MOI.dimension(set::ComplementsWithSetType) = set.dimension
