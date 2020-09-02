include("../src/Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using Printf
using Combinatorics
using .BAN

τ = 4;

#P = [-10, -10, -20, -20, -1, -10];
#P = [-10η, -10η, -1, -1, -η, -10η];
P = [-10, -10, -α, -α, -1, -10];

#    1 2 3 4 5 6
D = [0 1 0 0 0 1;
     1 0 1 0 1 1;
	 0 0 0 1 1 0;
	 0 0 1 0 1 0;
	 0 1 1 1 0 1;
	 1 1 0 0 1 0];

D = convert(Matrix{Bool}, D);



x = solve_patrolling(P, D, τ, verbose = true, genLatex = false);

println(x);