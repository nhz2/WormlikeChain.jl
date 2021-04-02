import Random123
using KernelAbstractions.Extras


"""
Return a tuple of 4 random UInt32 
Uses the philox random CBRNG 4x32_10
"""
function randbits_4x32(key::UInt64, ctr1::UInt64, ctr2::UInt64)::NTuple{4, UInt32}
    r= Random123.Philox4x{UInt32, 10}(0,0,0,0,
                                    key%UInt32, UInt32(key>>>32),
                                    ctr1%UInt32, UInt32(ctr1>>>32),
                                    ctr2%UInt32, UInt32(ctr2>>>32), 0)
    Random123.random123_r(r)
end

"""
Return a tuple of 2 random UInt64 
Uses the philox random CBRNG 4x32_10
"""
function randbits_2x64(key::UInt64, ctr1::UInt64, ctr2::UInt64)::NTuple{2, UInt64}
    x= randbits_4x32(key, ctr1, ctr2)
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


"""
Return a tuple of 4 normal random Float32 
Uses the philox random CBRNG 4x32_10 and box-muller transform
"""
function randn_4x32(key::UInt64, ctr1::UInt64, ctr2::UInt64)::NTuple{4, Float32}
    x= randbits_4x32(key, ctr1, ctr2)
    a= signed(x[1]) * 2.0f0^-31 # uniform [-1,1]
    b= x[2] * 2.0f0^-32 + 2.0f0^-33 # uniform (0,1]
    r= sqrt(-2*log(b))
    n1= r*sin(a*pi)
    n2= r*cos(a*pi)
    a= signed(x[3]) * 2.0f0^-31 # uniform [-1,1]
    b= x[4] * 2.0f0^-32 + 2.0f0^-33 # uniform (0,1]
    r= sqrt(-2*log(b))
    n3= r*sin(a*pi)
    n4= r*cos(a*pi)
    (n1,n2,n3,n4)
end

"""
Return a tuple of n normal random Float32, incrementing ctr2 for each 4.
Uses the philox random CBRNG 4x32_10 and box-muller transform
"""
function randn_32(key::UInt64, ctr1::UInt64, ctr2::UInt64, n)
    out= ()
    nrngs= (n-1)÷4 + 1
    @unroll for i in 0:(nrngs-1)
        out= (out..., randn_4x32(key, ctr1, ctr2+i)...)
    end
    out[1:n]
end