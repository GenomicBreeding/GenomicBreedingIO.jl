"""
    levenshteindistance(a::String, b::String)::Int64

Calculate the Levenshtein distance (edit distance) between two strings.

The Levenshtein distance is a measure of the minimum number of single-character edits 
(insertions, deletions, or substitutions) required to change one string into another.

# Arguments
- `a::String`: First input string
- `b::String`: Second input string

# Returns
- `Int64`: The minimum number of edits needed to transform string `a` into string `b`

# Examples
```jldoctest; setup = :(using GenomicBreedingIO)
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

Determines if two strings approximately match each other using Levenshtein distance.

The function compares two strings and returns `true` if they are considered similar enough
based on the Levenshtein edit distance and a threshold value. The threshold is applied as
a fraction of the length of the shorter string.

# Arguments
- `a::String`: First string to compare
- `b::String`: Second string to compare
- `threshold::Float64=0.3`: Maximum allowed edit distance as a fraction of the shorter string length

# Returns
- `Bool`: `true` if the strings match within the threshold, `false` otherwise

# Examples
```jldoctest; setup = :(using GenomicBreedingIO)
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
