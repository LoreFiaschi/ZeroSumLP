include("../src/Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
using .BAN

τ = 4;

#P = [-10, -10, -50, -50, -1, -10];
#P = [-10, -10, -20, -20, -1, -10];
#P = [-10, -10, -α, -α, -1, -10];
P = [-10η, -10η, -1, -1, -η, -10η];

#    1 2 3 4 5 6
D = [0 1 0 0 0 1;
     1 0 1 0 1 1;
	 0 0 0 1 1 0;
	 0 0 1 0 1 0;
	 0 1 1 1 0 1;
	 1 1 0 0 1 0];
	 
D = convert(Matrix{Bool}, D);

tol = Ban(0, ones(SIZE).*1e-5)

x, s = solve_patrolling(P, D, τ, verbose = true, eps=tol, genLatex = false);

x = denoise(x,1e-5);
for (xi, si) in zip(x, s)
	print("\t"); print(xi); print("\t"); println(si);
end