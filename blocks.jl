## System model: Functions and Types ##

using Interpolations, Unitful, QuadGK, CSV, Roots

# Importing the file that contains parameters
include("parameters.jl")

# Creating heat capacity functions
"""
** Heat capacity of liquid mercury **\n
Sourced from H.G. Lee materials thermodynamics
"""
function CpHg_l(T)
    (30.39u"J*mol^-1*K^-1" - 11.47e-3u"J*mol^-1*K^-2" * T)/M_Hg
end

function CpHg(T)
    if T<Tsat
        cp = CpHg_l(T)
    else
        cp = CpHg_g
    end
    return cp
end

"""
** Heat capacity of yttrium oxide powder **\n
Uses linear interpolation (with extrapolation) over data
sourced from [scientificgroupthermodataeuropesgteThermodynamicPropertiesCompounds2001]
"""
function CpPowder(T)
    # Heat capacity of Y₂O₃ from 
    # [scientificgroupthermodataeuropesgteThermodynamicPropertiesCompounds2001]
    temps =[300, 390, 480, 570, 660, 750, 840, 930, 1020, 1110, 1200]u"K"
    caps = [101.85,	110, 116.9, 119.2, 121.5, 123.8,
            124.95, 126.1, 127.25, 128.4, 128.4]u"J * mol^-1 * K^-1"/M_p
    interp = linear_interpolation(temps, caps, extrapolation_bc=Line())
    return interp(T)
end

function CpWater(T)
    temps = (convert(Array, 0:10:100) .+273.15)u"K"
    caps = [4.2176, 4.1921, 4.1818, 4.1784, 4.1785,
            4.1806, 4.1843, 4.1895, 4.1963, 4.2050, 4.2159]u"J*g^-1*K^-1"
    interp = linear_interpolation(temps, caps, extrapolation_bc=Line())
    return interp(T)
end

# Heat capacity of Al from H.G. Lee book 
function CpAl(T)
    (20.67u"J*mol^-1*K^-1" + 12.39e-3u"J*mol^-1*K^-2"*T)/M_Al
end

function CpGlass(T)
    (46.95u"J*mol^-1*K^-1" + 34.31e-3u"J*mol^-1*K^-2"*T - 11.3e-5u"J*mol^-1*K"*T^(-2))/M_glass
end

# Vapour pressure from [kroschwitzKirkOthmerEncyclopediaChemical2004] (only works from 0-150C, in kPa)
function v_P_Hg(T)
    exp((-3212.5/ustrip(T))+7.150)u"kPa"
end

# A super rough estimate of the condenser recovery based on the ratio of vapour pressure to partial pressure of mercury
# This approximates that the partial pressure of mercury is equivalent to the condenser pressure (low volume of nitrogen) 
"""
**recovery(Tout, Pcond) = r_condenser**\n
This is based on the concept of vapour pressure at the temperature of the
condenser as a fraction of total pressure
"""
function recovery(Tout, Pcond)
    1-(v_P_Hg(Tout)/(Pcond-v_P_Hg(Tout)))
end

# Creating the component type
"""
# Component

**comp(m, h, sp, pz)**\n
Structure with fields mass, specific enthalpy, species (string), phase (string)
"""
mutable struct comp
    m                       # Mass
    h                       # Specific enthalpy
    sp::AbstractString      # Species name
    pz::AbstractString      # Phase
    cp
    H
    # This picks the heat capacity based on sp and pz if I don't do it manually
    function comp(m, h, sp, pz)
        if sp == "Hg"
            cp = CpHg
        elseif sp == "N2"
            cp = CpN2
        elseif sp == "powder"
            cp = CpPowder
        elseif sp == "glass"
            cp = CpGlass
        elseif sp == "water"
            cp = CpWater
        elseif sp == "Al"
            cp = CpAl
        else println("Species or phase is named incorrectly; please define Cp")
        end
        # Also calculating the total energy for convenience
        H = m*h |> u"J"
        return new(m, h, sp, pz, cp, H) 
    end
end

# Creating the stream type
"""
# Stream

**strm([a::comp, b::comp], T)**\n
takes a vector of components and a temperature to make a stream object\n
(x::strm).cmp returns the component vector; T returns the temperature
"""
mutable struct strm
    cmp::Vector{comp}       # Vector of components
    T                       # Temperature of stream (assumed all streams are well mixed)
    ID
    function strm(cmp, T)
        new(cmp, T, 0)
    end
end

function streamenergy(stream::strm)
    total = 0u"J"
    for component in stream.cmp
        total += component.H
    end
    return total
end

function streamenergy(streams::Vector{strm})
    total = sum(streamenergy.(streams))
    return total
end

function streammass(stream::strm)
    total = 0u"g"
    for component in stream.cmp
        total += component.m
    end
    return total
end

function streammass(streams::Vector{strm})
    total = sum(streammass.(streams))
    return total
end

function energy_check(input::Vector{strm}, output::Vector{strm})
    Hin = streamenergy(input)
    Hout = streamenergy(output)
    ΔH = streamenergy(input) - streamenergy(output)
    if isapprox(ustrip(ΔH), 0, atol = 0.01)
    else
        throw(error("Energy balance does not close; Δ = $ΔH"))
    end
end

function energy_check(input::Vector{strm}, output::Vector{strm}, Q)
    ΔH = streamenergy(input) - streamenergy(output) + Q |> u"J"

    if isapprox(ustrip(ΔH), 0, atol = 0.01)
    else
        throw(error("Energy balance does not close; Δ = $ΔH"))
    end
end

function mass_check(input::Vector{strm}, output::Vector{strm})
    Δm = streammass(input) - streammass(output)
    if isapprox(streammass(input), streammass(output), rtol = 0.01)
    else
        throw(error("Mass balance does not close; Δ = $Δm"))
    end
end

# Each piece of equipment will be represented by a constructor that takes stream(s) and creates stream(s)
# If the output temperature is given and heat removal is calculated
"""
# Condenser
**condenser(input::strm, Tout, efficiency, Pcond) = [strm(vapour), strm(liquid), Qout]**
"""
function condenser(input::strm, Tout, efc, Pcond)
    # Tout is the output temperature of all streams (assuming well mixed)
    # efc is the combined efficiency of the condenser and chiller
    # r_C is the percentage of the mercury recovered by the condenser
    # This assumes that there are only 2 components in the input stream (Hg and N2 in this case)
    # Vapour output is given as 2, and liquid output as 3

    # Recovery fraction for mercury
    r_C = recovery(Tout, Pcond)

    # Parameters of component 1 in stream 2:
    comp21 = comp(input.cmp[1].m*(1-r_C) |> u"g", 
                  quadgk(input.cmp[1].cp, T0, Tout)[1],
                  input.cmp[1].sp,
                  "g")
    # Parameters of component 2 in stream 2:
    comp22 = comp(input.cmp[2].m |> u"g",
                  CpN2 * (Tout - T0), 
                  input.cmp[2].sp,
                  "g")
    # Parameters of component 1 in stream 3 (this assumes that only mercury condenses out):
    comp31 = comp(input.cmp[1].m*r_C |> u"g",
                  quadgk(input.cmp[1].cp, T0, Tout)[1] + dh_v_Hg,
                  input.cmp[1].sp,
                  "l")

    # Creating streams from components
    stream2 = strm([comp21, comp22], Tout)
    stream3 = strm([comp31], Tout)
    
    Qout = (comp21.m*(comp21.h - input.cmp[1].h) + comp22.m*(comp22.h - input.cmp[2].h)
            + comp31.m*(comp31.h-input.cmp[1].h))/efc

    # Checking energy balance and throwing error if it doesn't close
    Δ = (streamenergy([stream3, stream2]) - (streamenergy(input) + Qout*efc))
    if (streamenergy([stream3, stream2]) ≈ (streamenergy(input) + Qout*efc))
        else throw(error("energy does not balance, Δ = $Δ"))
    end
    # Return a vector of streams and the output temperature
    return [stream2, stream3, Qout]
end

"""
# Furnace
**furnace(feed::strm, flowgas::strm, Tf, efficiency, recovery) = [strm(vapour), strm(solid), Qf]**
"""
function furnace(in1::strm, in2::strm, Tf, eff, r_F)
    # Tf is the furnace temperature
    # r_F is the fraction of mercury removed from the solid waste
    # eff is the efficiency of the furnace
    # Input of powder or glass is in1, input of nitrogen is in2, output vapor is 3, output powder or glass is 4
    # Component 1 is mercury, 2 is nitrogen, 3 is the powder
    # Parameters of component 1 in stream 3
    comp31 = comp(in1.cmp[1].m*r_F,
                  dh_v_Hg + quadgk(in1.cmp[1].cp, T0, Tf)[1],
                  in1.cmp[1].sp,
                  "g")

    # Parameters of component 2 in stream 3
    comp32 = comp(in2.cmp[1].m,
                  in2.cmp[1].h + CpN2*(Tf - in2.T),
                  in2.cmp[1].sp,
                  "g")

    # Parameters of component 1 in stream 4
    comp41 = comp(in1.cmp[1].m*(1-r_F), 
                  quadgk(in1.cmp[1].cp, T0, Tf)[1],
                  in1.cmp[1].sp,
                  "l") # it is assumed that this is liquid, but it's really solid solution; this heat cap is technically wrong

    # Parameters of component 3 in stream 4
    comp43 = comp(in1.cmp[2].m,
                  # I don't know why the interpolation needs the (1) but if I remove it, it breaks
                  # It works with any number????
                  quadgk(in1.cmp[2].cp, T0, Tf)[1], 
                  in1.cmp[2].sp,
                  "s")
    if length(in1.cmp) == 3
        comp44 = comp(in1.cmp[3].m,
                      quadgk(in1.cmp[3].cp, T0, Tf)[1],
                      in1.cmp[3].sp,
                     "s")
        stream3 = strm([comp31, comp32], Tf)
        stream4 = strm([comp41, comp43, comp44], Tf)
    else
        stream3 = strm([comp31, comp32], Tf)
        stream4 = strm([comp41, comp43], Tf)
    end

    Qf = (streamenergy([stream3, stream4]) - streamenergy([in2]))/eff |> u"kJ"
    energy_check([in1, in2], [stream3, stream4], Qf*eff)

    # Note: Neglecting the nitrogen coming out with the powder when the system opens

    # Calculating the total heat input of the furnace
    # Don't need to subtract the input enthalpy because the input is at the ref temp T0
    return [stream3, stream4, Qf]
end

# Function for the carbon filter
"""
# Activated carbon filter
**cfilter(input) = [strm(output), strim(collection)]** \\
Calculates the output parameters for the 
[1] gas pass through (nitrogen) stream and
[2] collection stream (mercury)
of an activated carbon filter (assumed 100% capture rate) and returns them as a vector of streams
"""
function cfilter(input::strm)
    stream2 = strm([input.cmp[2]], input.T)
    stream3 = strm([input.cmp[1]], input.T)
    energy_check([input], [stream2, stream3])
    return [stream2, stream3]
end

"""
# Join two streams of mercury and nitrogen ideal gas

**mergestream(stream1::strm, stream2::strm) = stream3::strm**\n
- There can only be the two components in each stream
- Assume well mixed
"""
function mergestream(stream1::strm, stream2::strm)
    # see join_algebra.jl notbook for the equation used here
    # Pulling variables out of the stream structs 
    Cp_a = stream1.cmp[1].cp
    Cp_b = stream1.cmp[2].cp
    T1 = stream1.T; T2 = stream2.T
    m1a = stream1.cmp[1].m; m1b = stream1.cmp[2].m
    m2a = stream2.cmp[1].m; m2b = stream2.cmp[2].m
    h1a = stream1.cmp[1].h; h2a = stream2.cmp[1].h
    h1b = stream1.cmp[2].h; h2b = stream2.cmp[2].h

    # Calculating the output parameters
    m3a = m1a + m2a
    m3b = m1b + m2b

    # Calculating the output temperature
    # Objective function; enthalpy of mixing is neglected; assumed adiabatic mixing
    f(T) = m3a*(quadgk(Cp_a, T0, T)[1] + dh_v_Hg) + m3b*Cp_b*(T - T0) - 
    (m1a*(quadgk(Cp_a, T0, T1)[1] + dh_v_Hg) +
     m2a*(quadgk(Cp_a, T0, T2)[1] + dh_v_Hg) +
     m1b*Cp_b*(T1-T0) + m2b*Cp_b*(T2 - T0))
    # Finding the root
    T3 = find_zero(f, 500u"K")

    h3a = (quadgk(Cp_a, T0, T3)[1] + dh_v_Hg)
    h3b = (Cp_b*(T3 - T0))

    # Creating the stream output
    stream3 = strm([comp(m3a, h3a, stream1.cmp[1].sp, stream1.cmp[1].pz),
                    comp(m3b, h3b, stream1.cmp[2].sp, stream1.cmp[2].pz)], T3)
    # Throw an error if the energy balance doesn't close
    energy_check([stream1, stream2], [stream3])
    # Mass balance close check 
    mass_check([stream1, stream2], [stream3])
    return stream3
end

# Chiller takes coolant and required cooling as inputs and returns mass flow rate and streams
function chiller(coolant::comp, Qc)
end

"""
# Chiller

**chiller(Q, Cp_coolant, T_c, T_h, effChiller) = [stream_c::strm, strream_h::strm, Qchiller]**\n
- It may not be reasonable to assume that T_h is equivalent to the Tout of the condenser
- T_c must be less than T_h, which must be less than or equal to Tout
"""
function chiller(Q, Cpc, Tc, Th, effChiller)
    Δh = quadgk(Cpc, Th, Tc)[1]
    m = Q/Δh
    stream_c = strm([comp(m, absolute_h(Cpc, Tc) , "water", "l")], Tc)
    stream_h = strm([comp(m, absolute_h(Cpc, Th) , "water", "l")], Th)
    Qchiller = m*Δh/effChiller

    # Testing for energy close:
    if isapprox(m*(stream_c.cmp[1].h - stream_h.cmp[1].h), Qchiller*effChiller, rtol = 0.0005) 
    else
        Δ = m*(stream_c.cmp[1].h - stream_h.cmp[1].h) - Qchiller*effChiller
        throw(error("Chiller enthalpy does not close; Δ = $Δ, Qchiller = $Qchiller,
                    Δh = $Δh, Δh2 = $(stream_c.cmp[1].h - stream_h.cmp[1].h)"))
    end
    return [stream_c, stream_h, Qchiller]
end

function absolute_h(Cp::Function, T)
    h = quadgk(Cp, T0, T)[1]
end

function absolute_h(Cp, T)
    h = Cp*(T - T0)
end
