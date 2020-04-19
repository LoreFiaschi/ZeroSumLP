include("../I-Big-M/I_Big_M.jl")

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
    
    obj, x, base = I_Big_M(_A, _b, _c, _t, eps, verbose, genLatex);
    
    x /= sum(x);
    
    return x;
end