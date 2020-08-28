include("../src/Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
using .BAN

n = 4;
m = 25;

A = zeros(Ban, n, m);

tau = 3;

# row 1: upper left
# row 2: lower left
# row 3: right
# row 4: center

val_1 = 10;
val_2 = 50;
val_3 = 100;
val_4 = 1;

#=
val_1 = 10;
val_2 = 50;
val_3 = α;
val_4 = η;
=#

val = [val_1,val_2,val_3,val_4];

room_1 = [1,2,3,4,5,6];
room_2 = [16, 17, 18, 19, 20, 21, 22, 23, 24, 25];
room_3 = [11, 12, 13, 14, 15];
room_4 = [3,7,8,9,10];

room = [room_1, room_2, room_3, room_4];

n_room_1 = 6;
n_room_2 = 10;
n_room_3 = 5;
n_room_4 = 5

n_room = [n_room_1, n_room_2, n_room_3, n_room_4];

# notice: path 4 overwrites room 3 value
for i=1:n
	A[:,room[i]] .= val[i];
end

###########
tmp = A[:,20];
A[:,20] = A[:,12];
A[:,12] = tmp;

tmp = A[:,21];
A[:,21] = A[:,13];
A[:,13] = tmp;
###########

for i=1:n
	A[i,room[i]] .*= max(0, n_room[i]-tau)/n_room[i];
end

# To discount room 3 value on path 1 after overwriting
A[1,3] *= max(0, n_room[1]-tau)/n_room[1];

###########
#=
A[:,1] .= val_3;
A[1,1] *= max(0, n_room[1]-tau)/n_room[1];
=#
###########

#println(A);

tol = Ban(0, ones(SIZE).*1e-5)

x = solve(-A', verbose = true, eps=tol, genLatex = false);
print("\tx: "); println(x);
println("");


y = solve(A, verbose=false, eps=tol);
print("\ty: "); println(y);
println("");

print("Value: "); println(x'*A*y);
