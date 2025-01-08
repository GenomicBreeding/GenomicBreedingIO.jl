"""
    levenshteindistance(a::String, b::String)::Int64

Calculate the Levenshtein distance between 2 strings.
TODO: optimise as we only need to compute 2 rows at a time and do not need the full matrix

# Examples
```jldoctest; setup = :(using GBIO)
julia> levenshteindistance("populations", "populations")
0

julia> levenshteindistance("populations", "poplation")
2

julia> levenshteindistance("populations", "entry")
3
```
"""
function levenshteindistance(a::String, b::String)::Int64
    # a = "populations"; b = "entry";
    n::Int64 = length(a)
    m::Int64 = length(b)
    d::Matrix{Int64} = fill(0, n, m)
    for j = 2:m
        d[1, j] = j - 1
    end
    cost::Int64 = 0
    deletion::Int64 = 0
    insertion::Int64 = 0
    substitution::Int64 = 0
    for j = 2:m
        for i = 2:n
            if a[i] == b[j]
                cost = 0
            else
                cost = 1
            end
            deletion = d[i-1, j] + 1
            insertion = d[i, j-1] + 1
            substitution = d[i-1, j-1] + cost
            d[i, j] = minimum([deletion, insertion, substitution])
        end
    end
    d[end, end]
end

"""
    isfuzzymatch(a::String, b::String; threshold::Float64=0.3)::Bool

Fuzzy string matching using the Levenshtein distance, and a threshold as a fraction of the smaller string.

# Examples
```jldoctest; setup = :(using GBIO)
julia> isfuzzymatch("populations", "populations")
true

julia> isfuzzymatch("populations", "poplation")
true

julia> isfuzzymatch("populations", "entry")
false
```
"""
function isfuzzymatch(a::String, b::String; threshold::Float64 = 0.3)::Bool
    # a = "populations"; b = "populatins"; threshold = 0.3;
    n::Int64 = length(a)
    m::Int64 = length(b)
    dist::Int64 = levenshteindistance(a, b)
    dist < Int64(round(minimum([n, m]) * threshold))
end
