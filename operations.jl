#*#############################################################################
"""
    multiplySmart(A, B)

Compute the product of two matrices efficiently by determining the best 
factor order either `A * B` or `B * A`
"""
function multiplySmart(A::Any, B::Any)
    return A * B
end

#*#############################################################################
"""
    ×(A, B)

Binary operator definition for multiplySmart(A, B)
# Example
```julia-repl
julia> A × B

```
"""
×(A::Any, B::Any) = multiplySmart(A, B)

