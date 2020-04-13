include("../Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/BAN.jl")
using .BAN


A = zeros(Ban, 3, 4);

A[1,2] =  one(Ban) + 2*(one(Ban)>>1); A[1,3] = -one(Ban) - 2*(one(Ban)>>1);
A[2,1] = -one(Ban) - 2*(one(Ban)>>1); A[2,3] =  one(Ban) + 2*(one(Ban)>>1);
A[3,1] =  one(Ban) + 2*(one(Ban)>>1); A[3,2] = -one(Ban) - 2*(one(Ban)>>1);

A[1,4] = one(Ban)>>1; A[2,4] = -one(Ban)>>1;# A[3,4] = 2*one(Ban)>>1;


x = solve(-A');
print("\tx: "); println(x);
println("");

y = solve(A, verbose=false, genLatex = false);
print("\ty: "); println(y);
