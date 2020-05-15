include("../src/Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
using .BAN


A = zeros(Ban, 3, 3);

A[1,2] =  one(Ban); A[1,3] = -one(Ban);
A[2,1] = -one(Ban); A[2,3] =  one(Ban);
A[3,1] =  one(Ban); A[3,2] = -one(Ban);

x = solve(-A');
print("\tx: "); println(x);
println("");

y = solve(A, verbose=false, genLatex = false);
print("\ty: "); println(y);
