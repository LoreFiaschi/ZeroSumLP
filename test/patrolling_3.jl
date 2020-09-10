include("../src/Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
using .BAN

τ = 3;


#P = [-10, -10, -α, -α, -1, -10];
#P = [-10, -10, -α, -α, -η, -10];
#P = [-30η, -10η, -1, -1, -η, -10η];
#P = [-10η, -10η, -1, -1, -η*η, -10η];
#P = [-10, -10, -100, -100, -0.01, -10];
P = [-30, -10, -100, -100, -1, -10];

#    1 2 3 4 5 6
D = [0 1 0 0 0 1;
     1 0 1 0 1 1;
	 0 0 0 1 1 0;
	 0 0 1 0 1 0;
	 0 1 1 1 0 1;
	 1 1 0 0 1 0];



#=
#P = [-3, -10η, -1, -5η];
P = [-300, -10, -100, -5];


D = [0 1 0 1;
     1 0 1 0;
	 0 1 0 1;
	 1 0 1 0];
=#





D = convert(Matrix{Bool}, D);

tol = Ban(0, ones(SIZE).*1e-5)

x, s = solve_patrolling(P, D, τ, verbose = true, eps=tol, genLatex = false);

print("P: "); println(P); print("τ: "); println(τ); println("");

x = denoise(x,1e-5);
for (xi, si) in zip(x, s)
	print("\t"); print(xi); print("\t"); println(si);
end