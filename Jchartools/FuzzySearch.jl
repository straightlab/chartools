__precompile__()
module FuzzySearch
    include("FuzzySearch_.jl")
    export brute_force_search, brute_force_search_bestfuzzy, bestof2searches, brute_force_search_wildcard
end
