
using JuMP

# Solution reported by Sven Leyffer in MacMPEC:
# https://wiki.mcs.anl.gov/leyffer/index.php/MacMPEC
const MACMPEC_SOLUTIONS = Dict{Symbol, Float64}(
    :dempe_model => 28.25,
    :design_centering_model => -1.86065,
    :desilva_model => -1.0,
    :gauvin_model => 20.0,
    :qpec2_model => 45.0,
    :ralph_2_model => 0.0,
    :scale_1_model => 1.0,
    :scholtes4_model => -2e-4, # Here, we found a different solution than in MacMPEC (3.07336e-7)
    :water_net_model => 927.264324 # Solution reported is better than in MacMPEC (929.169)
)

# Ex. (2.2) from "Local convergence of SQP methods for MPECs".
# Solution is (0.5, 0.5), basic multipliers is (0, 1, 0)
function fletcher_leyffer_ex1_model()
    model = Model()
    @variable(model, z[1:2])
    set_lower_bound(z[2], 0)
    @objective(model, Min, (z[1] - 1)^2 + z[2]^2)
    @constraint(model, [z[2] - z[1], z[2]] ∈ MOI.Complements(2))
    return model
end

# Ex. (2.2) from "Local convergence of SQP methods for MPECs".
# Solution is (0, 1), dual (0.5, 1, 0, 0)
function fletcher_leyffer_ex2_model()
    model = Model()
    @variable(model, z[1:2] >= 0)
    @objective(model, Min, z[1] + z[2])
    @constraint(model, z[2]^2 >= 1)
    @constraint(model, [z[1], z[2]] ∈ MOI.Complements(2))
    return model
end

# dempe.mod	QQR2-MN-4-3
# Original AMPL coding by Sven Leyffer

# An MPEC from S. Dempe, "A necessary and sufficient optimality
# condition for bilevel programming problems", Optimization 25,
# pp. 341-354, 1992. From book Math. Progr. with Equil. Constr,
# by Luo, Pang & Ralph, CUP, 1997, p. 354.

# Number of variables:   2 + 1 multipliers
# Number of constraints: 2
# Nonlinear complementarity constraints
function dempe_model()
    model = Model()
    @variable(model, x, start=0.183193)
    @variable(model, z, start=0.428106)
    @variable(model, 0 <= w, start=3.00379)
    @objective(model, Min, (x - 3.5)^2 + (z + 4.0)^2)
    @constraint(model, z - 3.0 + 2.0 * z * w == 0.0)
    # Complementary constraint
    @constraint(model, [(x - z^2), w] ∈ MOI.Complements(2))
    return model
end

# design-cent-2.mod    QOR-MY-NLP-13-13-3
#
# Design centering problem cast as an MPEC, from an idea by
# O. Stein and G. Still, "Solving semi-infinite optimization
# problems with Interior Point techniques", Lehrstuhl C fuer
# Mathematik, Rheinisch Westfaelische Technische Hochschule,
# Preprint No. 96, November 2001.
#
# Maximize the volume of the parameterized body B(x) contained
# in a second body G, described by a set of convex inequalities.
#
# Original AMPL coding by Sven Leyffer, University of Dundee, Jan. 2002
function design_centering_model()
    x0 = [0.0, 0.0, 1.0]
    y0 = [
        -1.000000000000073 0.2425356250359245 0;
        0 0.9701425001468117 -1.000000000026019
    ]
    l0 = [0.5, 0.5153882031999911, 0.4999999999777709]
    model = Model()
    @variable(model, x[i=1:3], start=x0[i])
    @variable(model, y[j=1:2, k=1:3], start=y0[j, k])
    @variable(model, 0 <= l[k=1:3], start=l0[k])
    # ... maximize the volume of the inscribed body
    @objective(model, Min, -pi*x[3]^2)
    # ... lower level solutions lie in body G
    @constraint(model, g1, -y[1, 1] - y[2, 1]^2 <= 0.0)
    @constraint(model, g2, y[1, 2] / 4.0 + y[2, 2] <= 0.75)
    @constraint(model, g3, -y[2, 3] <= 1.0)
    # ... first order conditions for 3 lower level problem
    @constraint(model, 1.0 + 2.0 * (y[1, 1] - x[1]) * l[1] == 0.0)
    @constraint(model, 2*y[2, 1] + 2.0 * (y[2, 1] - x[2]) * l[1] == 0.0)
    @constraint(model, -0.25 + 2.0 * (y[1, 2] - x[1]) * l[2] == 0.0)
    @constraint(model, -1.00 + 2.0 * (y[2, 2] - x[2]) * l[2] == 0.0)
    @constraint(model, 0.00 + 2.0 * (y[1, 3] - x[1]) * l[3] == 0.0)
    @constraint(model, 1.00 + 2.0 * (y[2, 3] - x[2]) * l[3] == 0.0)
    # complementarity
    for k = 1:3
        @constraint(
            model,
            [-(y[1, k] - x[1])^2 - (y[2, k] - x[2])^2 + x[3]^2, l[k]] ∈ MOI.Complements(2)
        )
    end
    return model
end

# scholtes4.mod	LQR2-MN-3-2
# Original AMPL coding by Sven Leyffer

# An LPEC from S. Scholtes, Judge Inst., University of Cambridge.
function scholtes4_model()
    z0 = [0, 1]
    model = Model()
    @variable(model, z[i=1:2] >= 0.0, start=z0[i])
    @variable(model, z3, start=0.0)
    @objective(model, Min, z[1] + z[2] - z3)
    @constraint(model, -4 * z[1] + z3 <= 0)
    @constraint(model, -4 * z[2] + z3 <= 0)
    @constraint(model, [z[1], z[2]] ∈ MOI.Complements(2))
    return model
end

# Number of variables:   3 slack
# Number of constraints: 2

# desilva.mod	QQR2-MN-8-5
# Original AMPL coding by Sven Leyffer, University of Dundee

# An MPEC from F. Facchinei, H. Jiang and L. Qi, A smoothing method for
# mathematical programs with equilibrium constraints, Universita di Roma
# Technical report, 03.96. Problem number 5

# Number of variables:   6
# Number of constraints: 4
function desilva_model()
    model = Model()
    @variable(model, 0.0 <= x[1:2] <= 2.0)
    @variable(model, y[1:2])
    @variable(model, 0.0 <= l[1:2])
    @objective(model, Min, x[1]^2 - 2*x[1] + x[2]^2 - 2*x[2] + y[1]^2 + y[2]^2)
    for i = 1:2
        @constraint(model, 2.0*y[i] - 2.0*x[i] + 2.0 * (y[i] - 1.0) * l[i] == 0.0)
        @constraint(model, [0.25 - (y[i] - 1.0)^2, l[i]] ∈ MOI.Complements(2))
    end
    return model
end

function ralph_2_model()
    model = Model()
    @variable(model, x >= 0.0, start=1.0)
    @variable(model, y >= 0.0, start=1.0)
    @objective(model, Min, x^2 + y^2 - 4*x*y)
    @constraint(model, [x, y] ∈ MOI.Complements(2))
    return model
end

function scale_1_model()
    model = Model()
    @variable(model, x[1:2] >= 0.0)
    @objective(model, Min, (100.0 * x[1] - 1.0)^2 + (x[2] - 1.0)^2)
    @constraint(model, [x[1], x[2]] ∈ MOI.Complements(2))
    return model
end

#=
  MPEC example taken from Gauvin and Savard, "Working Paper G9037",
  GERAD, Ecole Polytechnique de Montreal (1992) (first version).

   Let Y(x) := { y | 0 <= y <= 20 - x }.
   min  f(x,y) := x**2 + (y-10)**2
   s.t.      0 <= x <= 15
             y solves MCP((F(y) := 4(x + 2y - 30)), Y(x))

  We have modified the problem by adding a dual variable and
  incorporating the constraint y <= 20 - x into the MCP function F.

  From a GAMS model by S.P. Dirkse & M.C. Ferris (MPECLIB),
  (see http://www.gams.com/mpec/).
=#
function gauvin_model()
    model = Model()
    @variable(model, 0 <= x <= 15, start=7.5)
    @variable(model, 0 <= y)
    @variable(model, 0 <= u, start=1.0)
    @objective(model, Min, x^2 + (y-10.0)^2)
    @constraint(model, [4 * (x + 2*y - 30) + u, y] ∈ MOI.Complements(2))
    @constraint(model, [(20 - x - y), u] ∈ MOI.Complements(2))
    return model
end

# MPEC of taxation model taken from (6)-(10) of
# Light, M., "Optimal taxation: An application of mathematical
# programming with equilibrium constraints in economics",
# Department of Economics, University of Colorado, Boulder, 1999.
# attributed to Hakonsen, L. "Essays on Taxation, Efficiency and the
# Environment", PhD thesis, Norwegian School of Economics & Business
# Administration, April 1998.
# TODO: require to add support for MOI.VectorNonlinearFunction
function hakonsen_model_broken()
    L = 100.0
    G = 25.0
    pL = 1.0
    model = Model()
    @variable(model, x[1:2] >= 0.0, start=1.0)
    @variable(model, l >= 0.0, start=1.0)
    @variable(model, 0.0 <= p[1:2])
    @variable(model, t[1:2] >= 0.0)
    @objective(model, Min, (x[1] * x[2] * l)^(1/3))
    @constraint(model, prices[i=1:2], [(pL - p[i]), x[i]] ∈ MOI.Complements(2))
    @constraint(
        model,
        consum2[i=1:2],
        [(x[i] * (3.0 * p[i] * (1+t[i])) - 100.0 * pL), p[i]] ∈ MOI.Complements(2)
    )
    @constraint(model, L*pL == sum(x[i] * p[i] for i = 1:2) + l*pL + G)
    @constraint(model, revenue, sum(p[i] * t[i] * x[i] for i = 1:2) >= G)
    return model
end

# qpec2.mod	QQR2-MN-v-v
# Original AMPL coding by Sven Leyffer, University of Dundee

# A QPEC from H. Jiang and D. Ralph, Smooth SQP methods for mathematical
# programs with nonlinear complementarity constraints, University of
# Melbourne, December 1997.

# Number of variables:  v
# Number of constraints: v
function qpec2_model()
    n = 10
    m = 20
    rr = -1.0
    ss = -2.0
    model = Model()
    @variable(model, x[1:n], start=1.0)
    @variable(model, y[1:m] >= 0.0, start=1.0)
    @variable(model, s[1:n] >= 0.0)
    @objective(model, Min, sum((x[i] + rr)^2 for i = 1:n) + sum((y[j] + ss)^2 for j = 1:m))

    @constraint(model, [i=1:n], y[i] - x[i] >= 0.0)
    @constraint(model, [i=1:n], [(y[i] - x[i]), y[i]] ∈ MOI.Complements(2))
    # This constraint is degenerate and should be replaced by y = 0
    @constraint(model, [i=(n+1):m], [y[i], y[i]] ∈ MOI.Complements(2))
    return model
end

#
# water-net.mod   OOR2-AN-x-x
# AMPL coding: Sven Leyffer, University of Dundee, April 2000.
#
# This version expresses the integer restrictions as complementarity
# i.e. the MINLP part,
#
#      var z{arcs} >= 0, <= 1, binary;
#
#      qpup{(i,j) in arcs}: qp[i,j] <= maxq * z[i,j];
#      qnup{(i,j) in arcs}: qn[i,j] <= maxq * (1 - z[i,j]);
#
# is modelled as complementarity betwen qn[i,j] and qp[i,j], the negative
# or positive flows/pressures in each arc. Thus discard z[i,j] and write
#
#      0 <= qp[i,j] complements qn[i,j] >= 0
#
# Here model complementarity with a nonlinear equation sum qn[i,j]*qp[i,j] = 0
#
# From an original MINLP GAMS model.
#
# THIS EXAMPLE ILLUSTRATES THE USE OF NONLINEAR PROGRAMMING IN THE DESIGN OF
# WATER DISTRIBUTION SYSTEMS. THE MODEL CAPTURES THE MAIN FEATURES OF AN
# ACTUAL APPLICATION FOR A CITY IN INDONESIA.
#
# REFERENCES: BROOKE A, DRUD A AND MEERAUS A, MODELING SYSTEMS AND
#             NONLINEAR PROGRAMMING IN A RESEARCH ENVIRONMENT,
#             IN R RAGAVAN AND S M ROHDE (EDS), COMPUTERS IN
#             ENGINEERING 1985, 1985.
#
#             DRUD D AND ROSENBORG A, DIMENSIONING WATER DISTRIBUTION
#                  NETWORKS, M. SC. DISSERTATION, (IN DANISH), INSTITUTE OF
#                  MATHEMATICAL STATISTICS AND OPERATIONS RESEARCH, TECHNICAL
#                  UNIVERSITY OF DENMARK, 1973.
function _from(arcs, i)
    return [(k, l) for (k, l) in arcs if k == i]
end
function _to(arcs, i)
    return [(k, l) for (k, l) in arcs if l == i]
end

function water_net_model()
    nodes = [:NW, :E, :CC, :W, :SW, :S, :SE, :N]
    reservoirs = [:NW, :E]
    consumers = [:CC, :W, :SW, :S, :SE, :N]
    arcs = [
        (:NW, :W),
        (:NW, :CC),
        (:NW, :N),
        (:E, :N),
        (:E, :CC),
        (:E, :S),
        (:E, :SE),
        (:CC, :W),
        (:CC, :SW),
        (:CC, :S),
        (:CC, :N),
        (:S, :SE),
        (:S, :SW),
        (:SW, :W),
    ]
    supply = Dict{Symbol,Float64}(:NW => 2.5, :E => 6.0)
    wcost = Dict{Symbol,Float64}(:NW => 0.2, :E => 0.17)
    pcost = Dict{Symbol,Float64}(:NW => 1.02, :E => 1.02)
    demand = Dict{Symbol,Float64}(
        :NW => 0.0,
        :E => 0.0,
        :CC => 1.212,
        :W => 0.452,
        :SW => 0.245,
        :S => 0.652,
        :SE => 0.252,
        :N => 0.456,
    )
    height = Dict{Symbol,Float64}(
        :NW => 6.50,
        :E => 3.25,
        :CC => 3.02,
        :W => 5.16,
        :SW => 4.20,
        :S => 1.50,
        :SE => 0.00,
        :N => 6.30,
    )
    x = Dict{Symbol,Float64}(
        :NW => 1200,
        :E => 4000,
        :CC => 2000,
        :W => 750,
        :SW => 900,
        :S => 2000,
        :SE => 4000,
        :N => 3700,
    )
    y = Dict{Symbol,Float64}(
        :NW => 3600,
        :E => 2200,
        :CC => 2300,
        :W => 2400,
        :SW => 1200,
        :S => 1000,
        :SE => 900,
        :N => 3500,
    )

    dist = Dict((i, j) => sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2) for (i, j) in arcs)

    dpow = 5.33      # ... power on diameter in pressure loss equation
    dmin = 0.15      # ... minimum diameter of pipe
    dmax = 2.00      # ... maximum diameter of pipe
    hloss = 1.03e-3  # ... constant in the pressure loss equation
    dprc = 6.90e-2   # ... scale factor in the investment cost equation
    cpow = 1.29      # ... power on diameter in the cost equation
    r = 0.10         # ... interest rate
    maxq = 2.00      # ... bound on qp and qn
    # ... average diameter (geometric mean)
    davg = sqrt(dmin * dmax)
    # ... ratio of demand to supply
    rr = sum(demand[i] for i in consumers) / sum(supply[i] for i in reservoirs)
    hl = Dict(i => height[i] for i in nodes)
    for i in nodes
        if i in consumers
            hl[i] += 7.5 + 5.0 * demand[i]
        end
    end

    model = Model()

    # N.B. : solution found is highly sensitive to the initial values of (qp, qn)
    @variable(model, 0 <= qp[arcs], start=10.0)
    @variable(model, 0 <= qn[arcs], start=10.0)
    @variable(model, dmin <= d[arcs] <= dmax, start=davg)
    @variable(model, h[i in nodes], start=hl[i] + 5.0, lower_bound=hl[i], upper_bound=100)
    @variable(model, 0 <= s[i in reservoirs] <= supply[i], start=rr*supply[i])

    @objective(
        model,
        Min,
        sum(s[i] * pcost[i] * (h[i] - height[i]) + s[i] * wcost[i] for i in reservoirs) /
        r +
        dprc * sum(dist[(i, j)] * d[(i, j)]^cpow for (i, j) in arcs) +
        sum(qp[(i, j)] + qn[(i, j)] for (i, j) in arcs),
    )
    # ... max flow through arcs
    @constraint(model, [(i,j) in arcs], qp[(i, j)] + qn[(i, j)] <= maxq)

    # ... flow conservation equation at each node
    @constraint(
        model,
        flow_cons[i in nodes],
        sum(qp[(j, i)] - qn[(j, i)] for (j, i) in _to(arcs, i)) -
        sum(qp[(i, j)] - qn[(i, j)] for (i, j) in _from(arcs, i)) +
        ((i in reservoirs) ? s[i] : 0.0) == demand[i]
    )
    # ... pressure loss on each arc (assumes qpow = 2)
    @constraint(
        model,
        loss[(i, j) in arcs],
        h[i] - h[j] ==
        hloss * dist[(i, j)] * (qp[(i, j)]^2 - qn[(i, j)]^2) / (d[(i, j)]^dpow),
    )
    # ... complementarity
    @constraint(model, [(i, j) in arcs], [qp[(i, j)], qn[(i, j)]] ∈ MOI.Complements(2))

    return model
end
