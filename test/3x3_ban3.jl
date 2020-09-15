include("../src/Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/matrix2latex.jl")
using .BAN

n = 3;

#A = rand(Ban,n,n);

A = [0.528857+0.159226η+0.514453η*η 0.0342537 + 0.506636η*η + 0.505368η 0.835119 + 0.568376η*η + 0.106677η;
     0.642693 + 0.902055η*η + 0.67181η 0.218219 + 0.664161η*η + 0.824841η 0.400041 + 0.64016η*η + 0.96803η;
     0.317266 + 0.584096η*η + 0.442962η 0.930651 + 0.554241η*η + 0.718353η 0.0387195 + 0.362937η*η + 0.602939η];

#matrix2latex(A);

tol = Ban(0, ones(SIZE).*1e-5)

verbose = false;
genLatex = true;

#=
x = solve(-A', verbose=verbose, eps=tol, genLatex=genLatex);
print("\tx: "); println(x);
#println(x'*A);
println("");
=#


y = solve(A, verbose=verbose, eps=tol)#, genLatex=genLatex);
print("\ty: "); println(y);



#println(A*y);
#=
active_y = findall(z->z>0, y);
n_active_y = length(active_y);


active_x = findall(z->z>0, x);
n_active_x = length(active_x);

b = zeros(Ban,n_active_y,1);
b[n_active_y] = 1;

C = Matrix{Ban}(undef,n_active_y,n_active_x);
AA = copy(A[active_x,active_y]);

for i=1:n_active_x-1
	C[i,:] = AA[:,i]-AA[:,i+1];
end

C[end, :] = ones(Ban,n_active_x);

x_C = C\b;
print("x_C: "); println(x_C);
println("");
println("");


b = zeros(Ban,n_active_x,1);
b[n_active_x] = 1;

R = Matrix{Ban}(undef,n_active_x,n_active_y);

for i=1:n_active_y-1
	R[i,:] = AA[i,:]-AA[i+1,:];
end

R[end,:] = ones(Ban,n_active_y);

y_R = R\b;
print("y_R: "); println(y_R);
println("");
println("");


print("{");
for i=1:n_active_y
	print("{");
	for j=1:n_active_x
		for h=1:SIZE
			if h==1
				print(C[i,j][h]);
			else
				if C[i,j][h]>=0
					print("+ $(C[i,j][h])/(g^$(h-1))");
				else
					print("$(C[i,j][h])/(g^$(h-1))");
				end
			end
		end
		if j!= n
			print(", ");
		elseif i!=n
			println("},");
		else
			println("}");
		end
	end
end
=#
