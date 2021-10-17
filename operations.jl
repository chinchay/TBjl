#*#############################################################################
"""
    smart_prod(A, B, C)

Compute the product of three matrices efficiently by determining the best 
factor order either `(A * B) * C` or `A * (B * C)`
"""
function smart_prod(A::Any, B::Any, C::Any)
    a1, a2 = size(A)
    b1, b2 = size(B)
    c1, c2 = size(C)
    @assert a2 == b1
    @assert b2 == c1
    
    # (A * B) * C
    complexity1 = (a1 * a2 * b2) + (a1 * b2 * c2)

    if a1 == a2 == b2 == c2
        return (A * B) * C, complexity1, -1, true
    end

    # A * (B * C)
    complexity2 = (b1 * b2 * c2) + (a1 * a2 * c2)

    if complexity1 < complexity2
        return (A * B) * C, complexity1, complexity2, true
    else
        return A * (B * C), complexity1, complexity2, false
    end
    
    # return (complexity1 <= complexity2) ? (A * B) * C : A * (B * C)
    
end

#*#############################################################################
"""
    ×(A, B, C)

Ternary operator definition for smart_prod(A, B, C)
# Example
```julia-repl
julia> A × B × C

```
"""
×(A::Any, B::Any, C::Any) = smart_prod(A, B, C)

