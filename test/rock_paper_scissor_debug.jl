include("../Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/BAN.jl")
using .BAN

A = ones(Ban, 3, 3);

A[1,2] *= 2; A[1,3] +=3;
A[2,1] *= 4; A[2,2] *= 5; A[2,3] *= 6;
A[3,1] *= 7; A[3,2] *= 8; A[3,3] *= 9;

x = solve(-A');
print("\tx: "); println(x);
println("");

y = solve(A, verbose=false, genLatex = false);
print("\ty: "); println(y);
