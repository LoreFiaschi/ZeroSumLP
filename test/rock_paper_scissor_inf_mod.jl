include("../Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/BAN.jl")
using .BAN


A = zeros(Ban, 3, 3);

#=
A[1,2] =  one(Ban) + 2*(one(Ban)>>1); A[1,3] = -one(Ban) - 2*(one(Ban)>>1);
A[2,1] = -one(Ban) -   (one(Ban)>>1); A[2,3] =  one(Ban) + 2*(one(Ban)>>1);
A[3,1] =  one(Ban) + 2*(one(Ban)>>1); A[3,2] = -one(Ban) - 2*(one(Ban)>>1);
=#

A[1,2] =  one(Ban); A[1,3] = -one(Ban);
A[2,1] = -one(Ban) -   (one(Ban)>>1); A[2,3] =  one(Ban);
A[3,1] =  one(Ban); A[3,2] = -one(Ban);

x = solve(-A')#, genLatex = true);
print("\tx: "); println(x);
#println(x'*A);
println("");

y = solve(A);
print("\ty: "); println(y);
#println(A*y);
