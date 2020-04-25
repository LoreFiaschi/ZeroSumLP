include("../Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/BAN.jl")
using .BAN


A = zeros(Ban, 4, 3);

#=
A[1,2] =  one(Ban); A[1,3] = -one(Ban);
A[2,1] = -one(Ban); A[2,3] =  one(Ban);
A[3,1] =  one(Ban); A[3,2] = -one(Ban);
A[4,1] =  0.5*one(Ban); A[4,2] = -one(Ban); A[4,3] = 0.5one(Ban);
=#
#
A[1,2] =  one(Ban) + 2*(one(Ban)>>1); A[1,3] = -one(Ban) - 2*(one(Ban)>>1);
A[2,1] = -one(Ban) - 2*(one(Ban)>>1); A[2,3] =  one(Ban) + 2*(one(Ban)>>1);
A[3,1] =  one(Ban) + 2*(one(Ban)>>1); A[3,2] = -one(Ban) - 2*(one(Ban)>>1);
A[4,1] =  0.5*one(Ban) - one(Ban)>>1; A[4,2] = -one(Ban) - 2*one(Ban)>>1; A[4,3] = 0.5*one(Ban) + one(Ban)>>1;
#

x = solve(-A', genLatex = true);
print("\tx: "); println(x);
println("");

y = solve(A, genLatex = true);
print("\ty: "); println(y);
