import Random123

"""
Return a tuple of 2 random UInt64 
Uses the philox random CBRNG 4x32_10
"""
function randbits_2x64(key::UInt64, ctr1::UInt64, ctr2::UInt64)::NTuple{2, UInt64}
    r= Random123.Philox4x{UInt32, 10}(0,0,0,0,
                                    key%UInt32, UInt32(key>>>32),
                                    ctr1%UInt32, UInt32(ctr1>>>32),
                                    ctr2%UInt32, UInt32(ctr2>>>32), 0)
    x= Random123.random123_r(r)
    (UInt64(x[1])|UInt64(x[2])<<32, UInt64(x[3])|UInt64(x[4])<<32)
end


"""
Return a tuple of 2 normal random Float64 
Uses the philox random CBRNG 4x32_10 and box-muller transform
"""
function randn_2x64(key::UInt64, ctr1::UInt64, ctr2::UInt64)::NTuple{2, Float64}
    x= randbits_2x64(key, ctr1, ctr2)
    a= signed(x[1]) * 2.0^-63 # uniform [-1,1]
    b= x[2] * 2.0^-64 + 2.0^-65 # uniform (0,1]
    r= sqrt(-2*log(b))
    r .* sincos(a*pi)
end
