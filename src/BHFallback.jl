"""
    arbitraryEjectaBH(m2_i, frac)

#TODO Description

# Arguments:
#TODO
- m2_i:
- frac:

# Output:
#TODO
- 
"""
function arbitraryEjectaBH(m2_i, frac)
    return m2_i*(1-frac)
end

"""
    restrictedEjectaBH(m2_i, frac; Ma=10, Mb=15, max_frac=1.0, min_frac=0.1)

#TODO Description

# Arguments:
#TODO
- m2_i:
- frac:
- Ma=10:
- Mb=15:
- max_frac=1.0:
- min_frac=0.1:

# Output:
#TODO
- 
"""
function restrictedEjectaBH(m2_i, frac; Ma=10, Mb=15, max_frac=1.0, min_frac=0.1)
    frac_limit = max_frac
    if m2_i > Ma && m2_i < Mb
        frac_limit = max_frac + (min_frac-max_frac)*(m2_i-Ma)/(Mb-Ma)
    elseif m2_i > Mb
        frac_limit = min_frac
    end
    return m2_i*(1-frac*frac_limit)
end
