## Basic system parameters ##
# These are more or less constant
R = 8.314u"J*mol^-1 * K^-1"         # Gas constant
dh_v_Hg = 294.68u"J*g^-1"           # Heat of vaporization for mercury
T0 = 300u"K"                        # This is the reference temperature (arbitrary)
h0 = 0u"J*g^-1"
M_N2 = 28.014u"g*mol^-1"
M_Hg = 200.59u"g*mol^-1"
M_p = 225.809u"g*mol^-1"            # Molar mass of the powder (just the mass of yttrium oxide)
M_glass = 60.08u"g/mol"
M_Al = 26.98u"g/mol"
CpN2 = 30u"J*mol^-1*K^-1"/M_N2
CpHg_g = (20786u"J*kmol^-1*K^-1")/M_Hg |> u"J*g^-1*K^-1"

# Specifying the parameters that are likely to change
mfeed = 2550.69u"kg"/3                          # Total mass of the feed; 8502 bulbs; each about 300g 
                                                # According to [kelleherFluorescentLightingOntario2007]
x_Al_feed = 0.045                               # Mass fraction Al in the feed 
x_Hg_Al = upreferred(3.3u"mg/kg")               # Mass fraction of mercury in the Al
x_Hg_powder = upreferred((7500)u"mg"/1000u"g")  # Mass fraction of mercury in the powder
                                                # from [nobleFluorescentLightTube]
x_feed_powder = 0.015                           # Mass fraction of powder in the feed
x_Hg_glass = upreferred(2.605u"mg*kg^-1")       # Mass fraction of mercury in the glass
                                                # from [nobleFluorescentLightTube]
x_glass_feed = 1 - x_Al_feed - x_feed_powder    # Mass fraction of glass in the feed
mglass = x_glass_feed*mfeed                     # Mass of glass in the feed
mpowder = mfeed*x_feed_powder                   # Mass of yttrium oxide in the feed powder
mAl = mfeed*x_Al_feed                           # Mass of aluminium in the feed
mHg_glass = mglass*x_Hg_glass
mHg_Al = mAl * x_Hg_Al
mHg_powder = x_Hg_powder * mpowder
m_Hg_total = mHg_glass + mHg_Al + mHg_powder

# Equipment efficiency
# Equipment does not have listed efficiency; leaving these as they are for now
η_condenser = 1       # Condenser efficiency
η_fp = 0.5             # Powder furnace efficiency
η_fg = 0.5              # Glass furnace efficiency
η_chiller = 0.5         # Chiller efficiency based on the 0.6 COP of the thermoelectric cooler; little bit extra for pump

Tf = (600+273.15)u"K"               # Powder furnace temperature
Tf_glass = (400+273.15)u"K"         # Glass furnace temperature
Pf = 10u"kPa"                       # Furnace pressure (system pressure)
Tout = (10+273.15)u"K"              # Condenser output temperature

pmax_fg = 60u"kW"                   # Maximum power of the chosen glass furnace
pmax_fp = 7u"kW"                    # Maximum power of the chosen powder furnace

# Tsat is a function of pressure, so this will need to change if the pressure does
# This is from FactSage
Tsat = 521u"K"

# Total cycle time
# Assuming these are equivalent for the time being
tcyc_powder = 6.5u"hr" 
tcyc_glass = 6.5u"hr"

# Calculating the mass flow rate of nitrogen based on the volumetric flow
# PV = nRT
# MPV/RT = m
# Volumetric flow rate of nitrogen gas according to Lee 2020
V_N2 = 2u"L*minute^-1"
m_N2 = M_N2*Pf*V_N2/(R*Tf) * tcyc_powder |> u"kg"
V_N2_satp = (m_N2/M_N2)*R*298.15u"K" / 1u"bar" |> u"L"

# Chiller stuff
Tc = (273.15+8.)u"K"         # Temperature of the cold water
Th = (273.15+9)u"K"          # Temperature of the hot water:w

carbon_eff = 0.95            # Removal rate of Hg in the carbon filter 
ccap = 0.2                   # Mercury capacity of activated carbon (mass frac)
