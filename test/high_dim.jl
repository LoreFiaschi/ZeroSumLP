include("../Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/BAN.jl")
using .BAN

n = 100;

A = rand(Ban, (n, n)).*rand([-1, 1], (n, n));

compute_time(A) = @time solve(A);

compute_time(-A');
nothing