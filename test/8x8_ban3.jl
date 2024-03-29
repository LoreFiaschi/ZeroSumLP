include("../src/Zero_Sum_LP.jl")
include("../../ArithmeticNonStandarNumbersLibrary/src/BAN.jl")
include("../../Utils/src/matrix2latex.jl")
using .BAN

n = 8;
#=
A = [0.650117+0.220683*η+0.0668177*η*η 0.0457915-0.0773716*η+0.00819438*η*η 0.123401+0.389841*η+0.024824*η*η 0.0382203+0.0850644*η+0.49097*η*η 0.192714-0.0349776*η+0.0893282*η*η 0.720847-0.245288*η+0.0284651*η*η -0.237809+0.580813*η+0.0438903*η*η -0.402455-0.276738*η+0.119441*η*η;
	-0.0457915+0.0773716*η-0.00819438*η*η 0.383273+0.612845*η*η+0.119763*η -0.0292925+0.213581*η*η+0.21913*η -0.00201577+0.127476*η*η-0.395549*η 0.202985-0.46358*η*η-0.152319*η 0.161758+0.0123527*η*η+0.532597*η 0.232+0.218402*η*η+0.221529*η -0.380975+0.0877249*η*η+0.758218*η;
	-0.123401-0.024824*η*η-0.389841*η 0.0292925-0.213581*η*η-0.21913*η 0.616721+0.836273*η*η+0.556601*η -0.176404-0.376775*η*η-0.00108258*η 0.0408527-0.674757*η*η-0.0129336*η -0.0456273+0.525442*η*η-0.167939*η -0.00251493+0.554466*η*η-0.0570219*η -0.0426986+0.474371*η*η+0.271489*η;
	-0.0382203-0.49097*η*η-0.0850644*η 0.00201577-0.127476*η*η+0.395549*η 0.176404+0.376775*η*η+0.00108258*η 0.0391397+0.327783*η*η+0.709921*η 0.330483+0.249965*η*η+0.165464*η 0.0603448-0.0210193*η*η-0.379307*η -0.0443614+0.0257967*η*η+0.214651*η -0.0832959-0.600408*η*η+0.149322*η;
	-0.192714-0.0893282*η*η+0.0349776*η -0.202985+0.46358*η*η+0.152319*η -0.0408527+0.674757*η*η+0.0129336*η -0.330483-0.249965*η*η-0.165464*η 0.800134+0.93119*η*η+0.0249308*η 0.0203673-0.270025*η*η-0.0517137*η -0.344858+0.0815255*η*η-0.235128*η 0.101667-0.0350513*η*η+0.67706*η;
	-0.720847-0.0284651*η*η+0.245288*η -0.161758-0.0123527*η*η-0.532597*η 0.0456273-0.525442*η*η+0.167939*η -0.0603448+0.0210193*η*η+0.379307*η -0.0203673+0.270025*η*η+0.0517137*η 0.545697+0.518061*η*η+0.167808*η 0.0199385-0.379621*η*η-0.537941*η -0.322996+0.0625491*η*η+0.728795*η;
	0.237809-0.0438903*η*η-0.580813*η -0.232-0.218402*η*η-0.221529*η 0.00251493-0.554466*η*η+0.0570219*η 0.0443614-0.0257967*η*η-0.214651*η 0.344858-0.0815255*η*η+0.235128*η -0.0199385+0.379621*η*η+0.537941*η 0.396595+0.982678*η*η+0.21497*η -0.225129+0.174438*η*η+0.157057*η;
	0.402455+0.276738*η-0.119441*η*η 0.380975-0.0877249*η*η-0.758218*η 0.0426986-0.474371*η*η-0.271489*η 0.0832959+0.600408*η*η-0.149322*η -0.101667+0.0350513*η*η-0.67706*η 0.322996-0.0625491*η*η-0.728795*η 0.225129-0.174438*η*η-0.157057*η -0.635429+0.320281*η*η-0.4090273*η];
=#

A = rand(Ban,n,n);

matrix2latex(A);

tol = Ban(0, ones(SIZE).*1e-5)

verbose = false;
genLatex = true;

x = solve(-A', verbose=verbose, eps=tol, genLatex=genLatex);
print("\tx: "); println(x);
#println(x'*A);
println("");

y = solve(A, verbose=verbose, eps=tol)#, genLatex=genLatex);
print("\ty: "); println(y);
#println(A*y);

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