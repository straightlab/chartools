"""
    brute_force_search(w::Vector{UInt8}, s::Vector{UInt8}, ix_start, ix_stop)
Searches for word `w` in template `s`. `w` and `s` are strings represented as UInt8 vectors for speed. `ix_start` and `ix_stop` limit the search to within these indices. Returns the number of exact matches `N_matches`. Also returns `pos`, the offset position of the word with respect to `ix_start` to find the match if the match is unique, or 0 otherswise (multiple matches or no match)
# Examples
```jldoctest
julia> w=b"julia"; s=b"I love julia language";

julia> brute_force_search(w,s,1,length(s))
(8, 1)

julia> brute_force_search(w,s,7,length(s))
(2, 1)

julia> brute_force_search(w,s,9,length(s))
(0, 0)

julia> brute_force_search(w,s,1,10)
(0, 0)

julia> brute_force_search(b"l",s,1,length(s))
(0, 3)
```
"""
function brute_force_search(w::Vector{UInt8}, s::Vector{UInt8}, ix_start::Int64, ix_stop::Int64) #pattern, template
    nw=length(w)
    s_cursor=ix_start
    s_cursor_max=ix_stop-nw+2
    j=s_cursor
    N_match=0
    pos=0
    while s_cursor<s_cursor_max
        for i=1:nw
            w[i]==s[j] ? j+=1 : break
        end
        if j==(s_cursor+nw)
            N_match+=1
            N_match==1 ? pos=s_cursor : pos=0
        end
        s_cursor+=1
        j=s_cursor
    end
    pos>0 ? pos-=(ix_start-1) : true
    return pos, N_match
end

function brute_force_search_wildcard(w::Vector{UInt8}, s::Vector{UInt8}, ix_start::Int64, ix_stop::Int64) #pattern, template, wildcard in word is 0x4e (N)
    nw=length(w)
    s_cursor=ix_start
    s_cursor_max=ix_stop-nw+2
    j=s_cursor
    N_match=0
    pos=0
    while s_cursor<s_cursor_max
        for i=1:nw
                ((w[i]==s[j]) || (w[i]==0x4e)) ? j+=1 : break
        end
        if j==(s_cursor+nw)
            N_match+=1
            N_match==1 ? pos=s_cursor : pos=0
        end
        s_cursor+=1
        j=s_cursor
    end
    pos>0 ? pos-=(ix_start-1) : true
    return pos, N_match
end


function brute_force_search_bestfuzzy(w::Vector{UInt8}, s::Vector{UInt8}, ix_start::Int64, ix_stop::Int64, maxerr::Int64=0)::NTuple{4,Int64} #pattern, template
    nw=length(w)
    s_cursor=ix_start
    s_cursor_max=ix_stop-nw+2
    j=s_cursor
    N_match=0
    N_fuzzy=0
    pos=0
    N_err=0
    while s_cursor<s_cursor_max
        n_err=0
        for i=1:nw
            if w[i]==s[j]
                j+=1
            else
                n_err+=1
                n_err>maxerr ? break : j+=1
            end
        end
        if j==(s_cursor+nw) && n_err<=maxerr
            if n_err==0
                N_err=0
                N_match+=1
                N_match==1 ? pos=s_cursor : pos=0
            else #n_err>0
                N_fuzzy+=1
                if N_match==0
                    if N_fuzzy==1 || N_err>n_err #if improve err then on
                        pos=s_cursor
                        N_err=n_err
                    else #N_fuzzy>1 and no improvement
                        pos=0
                    end
                end
            end
        end
        s_cursor+=1
        j=s_cursor
    end
    pos>0 ? pos-=(ix_start-1) : true
    return pos, N_err, N_match, N_fuzzy #pos is 0 unless N_match=1 | (N_match=0 & N_fuzzy=1) | (there is a single fuzzy match at the lowest err value). Return unique exact match pos (then n_err=0) or if none pos of unique best fuzzy match (then n_err>0). N_err is 0 if exact match, or the best error level found (even if there are multiple match at this error level)
end

function bestof2searches(s1::NTuple{4,Int64},s2::NTuple{4,Int64})::Tuple{UInt8,Int64,Int64,Int64,Int64}
    if s1[1]>0 && s2[1]>0
        if s1[2]<s2[2] #s1 wins
            bof=(0x02,s1[1],s1[2],s1[3]+s2[3],s1[4]+s2[4])
        elseif s1[2]>s2[2] #s2 wins
            bof=(0x03,s2[1],s2[2],s1[3]+s2[3],s1[4]+s2[4])
        else #tie, turns inot confusing out
            bof=(0x04,0,s1[2],s1[3]+s2[3],s1[4]+s2[4])
        end
    elseif s1[1]>0 #s1 has single best, s2 has either no match or confusing (s2[1]=0)
        if s2[2]<=s1[2]
            s2[3]+s2[4]>0 ? bof=(0x04,s2[1],s2[2],s1[3]+s2[3],s1[4]+s2[4]) : bof=(0x02,s1[1],s1[2],s1[3]+s2[3],s1[4]+s2[4]) # s2 has confusion and wins : s2 has no match so s1 wins
        else #s2[2]>s1[2] #s1 has no confusion and is better than s2 anyway, ezpz
            bof=(0x02,s1[1],s1[2],s1[3]+s2[3],s1[4]+s2[4])
        end
    elseif s2[1]>0 #s2 has single best, s1 has either no match or confusing (s1[1]=0)
        if s1[2]<=s2[2]
            s1[3]+s1[4]>0 ? bof=(0x04,s1[1],s1[2],s1[3]+s2[3],s1[4]+s2[4]) : bof=(0x03,s2[1],s2[2],s1[3]+s2[3],s1[4]+s2[4]) # s2 has confusion and wins : s2 has no match so s1 wins
        else #s2[2]>s1[2] #s1 has no confusion and is better than s2 anyway, ezpz
            bof=(0x03,s2[1],s2[2],s1[3]+s2[3],s1[4]+s2[4])
        end
    else #both s1 and s2 are either confusing or have no match
        if (s1[3]+s1[4]==0) && (s2[3]+s2[4]==0) #no match
            bof=(0x01,0,0,0,0)
        elseif (s1[3]+s1[4]>0) && (s2[3]+s2[4]>0) # both are confusing
            bof=(0x04,0,min(s1[2],s2[2]),s1[3]+s2[3],s1[4]+s2[4])
        else #only one of them is confusin
            (s1[3]+s1[4]==0) ? bof=(0x04,0,s2[2],s2[3],s2[4]) : bof=(0x04,0,s1[2],s1[3],s1[4])
        end
    end
    return bof
end
