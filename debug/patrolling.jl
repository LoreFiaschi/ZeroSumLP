include("../src/Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")

using Printf
using .BAN

τ = 4;

P = [10, 10, 20, 20, 1, 10];

#    1 2 3 4 5 6
D = [1 1 0 0 0 1;
     1 1 1 0 1 1;
	 0 0 1 1 1 0;
	 0 0 1 1 1 0;
	 0 1 1 1 1 1;
	 1 1 0 0 1 1];

D = convert(Matrix{Bool}, D);

tol = Ban(0, ones(SIZE).*1e-5);

u, f = solve_patrolling(P, D, τ, verbose = true, eps=tol, genLatex = false);
print("\tu: "); println(u);
println("");
print("\tf: "); println(f);