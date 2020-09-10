include("../../I-Big-M/src/I_Big_M.jl")
using Combinatorics

function solve(A::Union{Matrix{T}, Transpose{T, M}, Adjoint{T, M}}; 
                eps::Number=convert(promote_type(T,Float64),1e-5),verbose::Bool=false,genLatex::Bool=false) where {T, M <: AbstractMatrix}

    #                                                #
    # Matrix A contains the losses of the row player #
    #                                                #
    
    _A = copy(A);
    
    # Make all the matrix entries non-negative
    any(x->x<0, _A) && (_A = _A.-minimum(_A));
    
    n_row, n_col = size(_A);
   
    _b =  ones(T, n_row, 1);
    _c = -ones(T, n_col, 1);
    _t =  ones(Int64, n_row);
    
    obj, x, base = I_Big_M(_A, _b, _c, _t, eps=eps, verbose=verbose, genLatex=genLatex);
    
    x /= sum(x);
    
    return x;
end


function solve_patrolling(P::Vector{T}, D::AbstractMatrix{Bool}, τ::Integer;
                eps::Number=convert(promote_type(T,Float64),1e-5),verbose::Bool=false,genLatex::Bool=false) where T<:Number


#### with disposition

	n = length(P);
	
	######################
	# Variables creation #
	######################
	
	states = multiset_permutations([i for i=1:n for j=1:τ],τ);
	to_keep = [];
	
	# delete (some) unconsistent states
	for (i,s) in enumerate(states)
	
		consistent = true;

		for j=1:τ-1
			!D[s[j+1],s[j]] && (consistent=false; break;)
		end
		
		#consistent && all(x->x==false, D[setdiff(1:n,s),s[end]]) && (consistent=false);
		consistent && push!(to_keep,i); 
	end
	
	states = collect(states)[to_keep];
	
	n_v = length(to_keep);
	
	# get the states without transitions
	states_no_transition = unique(map(x->x[2:end], states));
	n_s = length(states_no_transition);
	
	
	
	######################
	# Payoff constraints #
	######################
	
	
	A = zeros(T,n+n_s,n_v);
	
	for i=1:n_v
		A[setdiff(1:n,states[i][2:end]).+n_s,i] = copy(P[setdiff(1:n,states[i][2:end])]);
	end
    
    # Make all the matrix entries non-negative
    any(x->x<0, A[n_s+1:n_s+n,:]) && (A[n_s+1:n_s+n,:] = A[n_s+1:n_s+n,:].-minimum(A[n_s+1:n_s+n,:]));
	
	
	
	
	##################
	# Flux balancing #
	##################
	
	for i=1:n_s
		A[i,findall(x->x[2:end]==states_no_transition[i], states)] .= 1;
		A[i,findall(x->x[1:end-1]==states_no_transition[i], states)] .= -1;
	end
	
	
   
    b =  [zeros(T,n_s,1); ones(T,n, 1)];
	
	t = [zeros(Int64,n_s); ones(Int64,n)];
	
	
	
	
	#################
	# Cost function #
	#################
	
	
	c = -ones(T,n_v,1);

	A = convert(Matrix{Ban},A);
	b = convert(Matrix{Ban},b);
	c = convert(Matrix{Ban},c);
	
	
	obj, x, base = I_Big_M(A, b, c, t, eps=eps, verbose=verbose, genLatex=genLatex);
	
	#=
	@show states
	println("");
	=#
	
	return x/sum(x), states

end






























#=
function solve_patrolling(P::Vector{T}, D::AbstractMatrix{Bool}, τ::Integer;
                eps::Number=convert(promote_type(T,Float64),1e-5),verbose::Bool=false,genLatex::Bool=false) where T<:Number

    #                                                #
    # Vector P contains the losses of the row player #
    #                                                #
    # Matrix D is the adjacency matrix 				 #
	#												 #

	n_nodes = length(P);
	
	n_variables = n_nodes^2+n_nodes*(τ+1);
	
	# Variables -- u -- they are n_nodes and represent the probability the patroller enters in node i after more than τ istants from the last time
	
	# Variables -- x -- they are n_nodes^2 and represent the probability the patroller moves from i to j
	
	# Variables -- z -- they are τ*n_nodes and represent the probability the patroller enters in node i after less 1,2,...,τ istants (they are auxiliary)

	
	######################
	# Payoff constraints #
	# Variables -- u --  #
	######################
	
    A_payoff = Matrix{T}(undef,n_nodes,n_nodes);
	
	for i=1:n_nodes
		A_payoff[i,:] = copy(P);
		A_payoff[i,i]=0;
	end
    
    # Make all the matrix entries non-negative
    any(x->x<0, A_payoff) && (A_payoff = A_payoff.-minimum(A_payoff));
   
    b_payoff =  ones(T,n_nodes, 1);
	
	t_payoff = ones(n_nodes);
	
	#=
	@show A_payoff;
	println("");
    =#
	
	
	
	##########################
	# Traversing constraints #
	# Variables -- x --      #
	##########################
	
	holes = findall(.!D);
	n_holes = length(holes);
	A_traversing = zeros(T,n_holes,n_nodes^2);
	for i=1:n_holes
		A_traversing[i,(holes[i][1]-1)*n_nodes+holes[i][2]] = 1;
	end
	
	b_traversing = zeros(T,n_holes);
	
	t_traversing = zeros(n_holes);
	
	#=
	@show A_traversing;
	println("");
	=#
	
	
	
	#############################
	# Normalization constraints #
	# Variables -- x --			#
	#############################
	
	A_normalization = ones(T,1,n_nodes^2);
	
	b_normalization = 1;
	
	t_normalization = 0;
	
	
	
	
	#####################
	# Flux constraints  #
	# Variables -- x -- #
	#####################
	
	A_flux = zeros(T,n_nodes,n_nodes^2);
	for i=1:n_nodes
		A_flux[i, [(j-1)*n_nodes+i for j=1:n_nodes]] .= -1
		A_flux[i,(i-1)*n_nodes+1:i*n_nodes] .= 1;
		A_flux[i,(i-1)*n_nodes+1] = 0;
	end
	
	b_flux = zeros(T,n_nodes,1);
	
	t_flux = zeros(n_nodes);
	
	#=
	@show A_flux;
	println("");
	=#
	
	
	
	###########################
	# Reciprocity constraints #
	# Variables -- x,z --     #
	###########################
	
	A_reciprocity = zeros(T,n_nodes,n_nodes^2+τ*n_nodes);
	for i=1:n_nodes
		A_reciprocity[i,n_nodes^2+i] = -1;
		A_reciprocity[i,(i-1)*n_nodes+i] = 1;
	end
	
	b_reciprocity = zeros(T,n_nodes,1);
	
	t_reciprocity = zeros(n_nodes);
	
	#=
	@show A_reciprocity;
	println("");
	=#
	
	
	
	#########################
	# Longrun constraints   #
	# Variables -- x,z -- #
	#########################
	
	A_longrun = zeros(T,n_nodes,n_nodes^2+τ*n_nodes);
	
	for i=1:n_nodes
		A_longrun[i,[n_nodes^2+(j-1)*n_nodes+i for j=1:τ]] .= -1;
		A_longrun[i,[(j-1)*n_nodes+i for j=1:n_nodes]] .= 1;
	end
	
	b_longrun = zeros(T,n_nodes,1);
	
	t_longrun = zeros(n_nodes);
	
	
	
	
	#######################
	# Density constraints #
	# Variables -- z --   #
	#######################
	
	A_density = zeros(T,n_nodes,τ*n_nodes);
	coef = [i for i=1:τ];
	for i=1:n_nodes
		A_density[i,[(j-1)*n_nodes+i for j=1:τ]] = coef;
	end
	
	b_density = ones(T,n_nodes,1);
	
	t_density = -ones(n_nodes);
	
	#=
	@show A_density;
	println("");
	=#
	
	
	###############################
	# Complementarity constraints #
	# Variables -- u,z --         #
	###############################
	
	A_complementarity_z = zeros(T,n_nodes,τ*n_nodes);
	for i=1:n_nodes
		A_complementarity_z[i,[(j-1)*n_nodes+i for j=1:τ-1]] .= -1;
	end
	
	A_complementarity_u = I;
	
	b_complementarity = zeros(T,n_nodes,1);
	
	t_complementarity = zeros(n_nodes);
	
	#=
	@show A_complementarity_z;
	println("");
	=#
	
	
	
	#################
	# Cost function #
	#################
	
	
	c = zeros(T,n_variables,1);
	c[1:n_nodes] .= -1;
	c = convert(Matrix{Ban},c);
	
	
	
	
	#####################
	# Constraints types #
	#####################

    t =  [t_payoff; t_traversing; t_normalization; t_flux; t_reciprocity; t_longrun; t_density; t_complementarity];
	t = convert(Array{Int64},t)
	
	
	
	
	#####################
	# Constraint matrix #
	#####################	
	
	#        u                                x                          z
	A = [A_payoff                 zeros(T,n_nodes,n_nodes^2) zeros(T,n_nodes,τ*n_nodes);
	     zeros(T,n_holes,n_nodes) A_traversing               zeros(T,n_holes,τ*n_nodes);
		 zeros(T,1,n_nodes) 	  A_normalization            zeros(T,1,τ*n_nodes);
		 zeros(T,n_nodes,n_nodes) A_flux                     zeros(T,n_nodes,τ*n_nodes);
		 zeros(T,n_nodes,n_nodes)                   A_reciprocity;
		 zeros(T,n_nodes,n_nodes)                     A_longrun
		 zeros(T,n_nodes,n_nodes) zeros(T,n_nodes,n_nodes^2) A_density;
		 A_complementarity_u      zeros(T,n_nodes,n_nodes^2) A_complementarity_z];
	
	A = convert(Matrix{Ban},A);
	
	
	########################
	# Costant terms vector #
	########################
	
	b = [b_payoff; b_traversing; b_normalization; b_flux; b_reciprocity; b_longrun; b_density; b_complementarity];
	b = convert(Matrix{Ban},b);
	
	#=
	println(size(A));
	println(size(b));
	println(size(c));
	println(size(t));
	=#
	
    
    obj, x, base = I_Big_M(A, b, c, t, eps=eps, verbose=verbose, genLatex=genLatex);
    
    u = x[1:n_nodes]/sum(x[1:n_nodes]);
	f = x[n_nodes+1:n_nodes^2];
    
    return u, f;
	
end
=#