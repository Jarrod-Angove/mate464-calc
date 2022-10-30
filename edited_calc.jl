using Markdown
using InteractiveUtils

using Unitful 		# This lets us use units in calculations

# Creating an interpolation function for the heat capacity of water
begin
using Interpolations
# Temps from 0 to 100C converted to K
temps = (convert(Array, 0:10:100) .+273.15)u"K"

# Cp in J/(g K) from CRC handbook
caps = [4.2176, 4.1921, 4.1818, 4.1784, 4.1785, 4.1806, 4.1843, 4.1895, 4.1963, 4.2050, 4.2159]u"J * g^-1 * K^-1"

# Create a linear interpolation function from the data above
Cpc_l = linear_interpolation(temps, caps)

# Note that you must input the temp as Kelvin; Here is an example using this function
Cpc_l((55+273.15)u"K")
end;

using DataFrames, LaTeXStrings, Latexify, CSV

md"""
## Establishing some variables
"""

begin
# All uncertainties here are rough estimates based on variations in the literature
m_feed = (500)u"kg"; 						# Feed mass (arbitrary)
X_Hg = upreferred((3575)u"mg"/1000u"g"); 	# Mass fraction of Hg in feed powder 													(Average from Park 2016)
m_Hg = X_Hg * m_feed; 							# Total mass of mercury
M_Hg = 200.59u"g*mol^-1"; 						# Molar mass of mercury
P_f = 10.0u"kPa"; 								# Furnace pressure (Lee 2020)
R = 8.314u"J*mol^-1 * K^-1"; 	 				# Gas constant
T_f = (600+273.15)u"K"; 						# Furnace temperature (Lee 2020)
V_sys = 5u"m^3"; 								# Assumed system volume (arbitrary)
M_N2 = 28.014u"g*mol^-1"; 						# Molar mass of N₂
t_cyc = 390u"minute" |> u"s" 					# Total cycle time from Lee 2020
P_atm = 1.0u"atm"|>u"kPa" 						# Atmospheric pressure convert to kPa
ṁ_Hg = m_Hg / t_cyc 							# Assume constant evaporation
	
# Taking molar mass of powder as molar mass of Y₂O₃ (rough estimate based on the fact that Y₂O₃ makes up ~70% of powder mass)
M_p = 225.809u"g*mol^-1" 
m_p = m_feed*(1-X_Hg)	# Mass of the powder without the mercury

T0 = 300.0u"K" 			# Temp at which everything enters the system
Tout = (10+273.15)u"K"; # Assumed condenser output temperature is 10°C (Lee 2020)

# Volumetric flow rate will be dictated by our pump choice for the assumed furnace pressure P_f; temporarily using value from Lee 2020 as an assumption
V̇_N2 = (2)u"L*minute^-1" |> u"m^3 * s^-1";

end;

Tout

V_tot = m_Hg/(M_Hg * P_f) * R * T_f |> u"m^3" # Total volume of Hg gas

# Total mass of nitrogen required to fill the system 
# (rough estimate since we don't know the total system volume yet)
m_n2_sys = (M_N2 * P_f * V_sys)/(R * Tout) |> u"g" 

m_Hg |>u"g" # Total mass of mercury in grams

md"""
# Mass Balances
"""

md"""
## Heat capacity stuff
"""

# Heat capacity of liquid mercury as a function of temperature
CpHg_l(T) = 30.39u"J*mol^-1*K^-1" - 11.47e-3u"J*mol^-1*K^-2" * T ;
# From Hae-Geon Lee Materials thermodynamics pg. 431

# Ideal gas heat capacity of mercury
CpHg_g = 20786u"J*kmol^-1*K^-1"

# Ideal gas heat capacity of nitrogen from CRC (rough approximation)
CpN2 = 30u"J*mol^-1*K^-1"

# Creating an interpolation function based on data from scientificgroupthermodataeuropesgteThermodynamicPropertiesCompounds2001
# This is for Y2O3 but we are using it as an assumption
begin
capp = [101.85,	110, 116.9, 119.2, 121.5, 123.8, 124.95, 126.1, 127.25, 128.4, 128.4]u"J * mol^-1 * K^-1";
tempp = [300, 390, 480, 570, 660, 750, 840, 930, 1020, 1110, 1200]u"K";
	
Cpp = linear_interpolation(tempp, capp);
end;

md"""
## Furnace
"""

# I am assuming the pump flow rate is calculated at T=Tout as this is what the pump is pulling at inlet.
# This is assuming ideal gas (literally just the ideal gas law here)
# Should not be using P_f for this, as the pressure at the pump inlet will be lower than in the furnace, so we will change this when we do pressure loss calculations for equipment selection.
ṁ_N2 = (M_N2 * P_f * V̇_N2)/(R * Tout) |> u"g * s^-1"

# Convert volumetric to mass flow rate for nitrogen, assuming 2L/min
begin
# Saturation temperature at 10kPa as found by Jessica from factsage
# This will need to be changed if we use a different pressure!!!
# I would like to create a function for this in the future...
Tsat = 521u"K"

# Assume constant flow rate
m_N2 = ṁ_N2 * t_cyc

# Mass fraction of nitrogen and mercury in the stream leaving the furnace
X_N2 = (ṁ_N2)/(ṁ_Hg+ṁ_N2)
X_Hg_F = 1-X_N2
end;

begin
using QuadGK # This is a numerical integration package
# Using numerical integration over heat capacity of liquid Hg
ΔH_Hg_l = quadgk(CpHg_l, T0, Tsat, rtol=1e-8)[1] * m_Hg / M_Hg |> u"kJ"
end

# Volumetric flow rate in the furnace is different becuase the temperature changes while the mass and pressure remain constant
V̇_N2_F = (ṁ_N2 * R * T_f)/(M_N2 * P_f) |> u"m^3 * s^-1"

md"""
### Enthalpy change of mercury
"""

ΔH_Hg_g = CpHg_g * (T_f-Tsat) * m_Hg / M_Hg |> u"kJ"

# Heat of vaporization of mercury from CRC handbook
# This is at 1atm pressure, so it may not be accurate
ΔH_v_Hg = 59.11u"kJ*mol^-1" * m_Hg / M_Hg |> u"kJ"

ΔH_Hg_F = ΔH_Hg_l + ΔH_Hg_g + ΔH_v_Hg

md"""
### Enthalpy change of powder
"""

# numerically integrating over the interpolation function created above
# Assuming there are no phase changes
ΔH_powder = quadgk(Cpp, T0, T_f)[1] * m_p / M_p |> u"kJ"
# Note that the poweder leaves at the same temperature, so this enthalpy is not relavent to the steady state system

md"""
### Enthalpy change of nitrogen
"""

ΔH_N2_F = m_N2 * CpN2 * (T_f - Tout) / M_N2 |> u"kJ"

md"""
### Enthalpy of the furnace
"""

# Total energy input of the furnace for a single batch
ΔH_F = ΔH_powder + ΔH_N2_F + ΔH_Hg_F

# Power required during continuous portion of the process
ΔḢ_F_steady = (ΔH_N2_F + ΔH_Hg_F)/t_cyc |> u"W"

md"""
This is an extremely small furnace, even if we account for efficiency losses!
"""

md"""
## Condenser Balance
"""

md"""
Since everything else is known or assumed, this will be a calculation of our coolant (water) mass flow rate for a given temperature
"""

# Assuming the coolant enters the system at 5C
Tc_in = (5+273.15)u"K"

ΔH_N2_c = CpN2 * m_N2 * (Tout - T_f) / M_N2 |> u"kJ" 
# Input temp is the temp of the furnace T_f

ΔH_Hg_c = (-ΔH_Hg_g - ΔH_v_Hg + m_Hg * quadgk(CpHg_l, Tsat, Tout)[1]/M_Hg)  |> u"kJ" 

Δh_c = quadgk(Cpc_l, Tout, Tc_in)[1]

# Mass flow rate of coolant required for given parameters
m_c = (ΔH_N2_c + ΔH_Hg_c)/Δh_c 

ΔH_c = m_c*quadgk(Cpc_l, Tout, Tc_in)[1]

ṁ_c = m_c/t_cyc |> u"g*s^-1"

# Checking if the condenser energy balance closes
ΔH_N2_c + ΔH_Hg_c ≈ ΔH_c

md"""
Therefore, the mass flow rate of the coolant is $0.34 \frac{g}{s}$
"""

md"""
Alternatively, we can find the condenser power:
"""

ΔH_chiller = -(ΔH_N2_c + ΔH_Hg_c) |> u"kJ" 
# This is the power  required by the chiller unit before efficiency losses

md"""
 $\therefore$ power requirements of the condensor are trivial
"""

md"""
## Pump stuff
"""

# Rough estimate of the pressure difference generated by the pump
# This will change if we consider pressure losses during equipment selection
ΔP_p = P_atm - P_f |> u"kPa"

# Work required by pump before accounting for efficiency
W_p = V̇_N2 * ΔP_p |> u"W"

# Enthalpy change from pump (i.e) ΔH(s₄ --> s₅)
ΔH_p = W_p * t_cyc |> u"kJ"
# This isnt being considered in the enthalpy balance as we assume the valve K-112 reverses this input perfectly

md"""
## Table creation
"""

md"""
For all streams:
- Temperature
- Species; Phase
- Mass flow rate
"""

enths = DataFrame(
"Name" =>
[L"\Delta H_{\mathrm{Hg,C}}",
L"\Delta H_{\mathrm{N_2, C}}",
L"\Delta H_\mathrm{Hg, F}",
L"\Delta H_{\mathrm{N_2, F}}",
L"\Delta H_{\mathrm{p}}",
L"\Delta H_\mathrm{C}", 
L"\Delta H_\mathrm{F}" ], 
"Value (kJ)" => 
Measurements.value.(ustrip.(u"kJ", [
	ΔH_Hg_c,
	ΔH_N2_c,
	ΔH_Hg_F,
	ΔH_N2_F,
	ΔH_powder,
	ΔH_c,
	ΔH_F
])),
"Description" =>
[
	"Mercury in condenser",
	"Nitrogen in condenser",
	"Mercury in furnace",
	"Nitrogen in furnace",
	"Powder in furnace (batch)",
	"Total condenser change",
	"Total furnace change"
]
)

md"""
Mass flow rates for each stream:
"""

begin
	mn2 = ṁ_N2
mfrs::Vector{Any} = round.(ustrip.(u"mg*s^-1", upreferred.([
		m_feed/t_cyc,
		ṁ_N2,
		mn2,
		mn2,
		ṁ_c, 
		ṁ_c,
		m_n2_sys/t_cyc,
		mn2,
		mn2,
		mn2,
		m_n2_sys/t_cyc,
		ṁ_Hg,
		m_p/t_cyc])), sigdigits = 4)
insert!(mfrs, 2, round.(ustrip.(u"mg*s^-1",[ṁ_N2, ṁ_Hg]), sigdigits = 4))
end

# Exporting the enthalpy table
CSV.write("enthalpy.csv", enths);

md"""
Putting the calculated variables into a stream summary table:
"""

begin
	n2 = "N₂"
	h20 = "H₂O"
	hg = "Hg"
	powder = "Y₂O₃ (powder)"
	nd = ", "
	s = "solid"
	l = "liquid"
	g = "gas"
	mn3 = ṁ_N2
stream_table = DataFrame(
	"ID" => 1:14,
	"Species" => [
		powder*nd*hg, n2*nd*hg, n2, n2, n2, 
		h20, h20, n2, n2, n2, n2, n2, hg, powder],
	"Phase" => [
		s, g, g, g, g, l, l, g, g, g, g, g, l, s],
	"Temperature (K)" => ustrip.([
		T0, T_f, Tout, Tout,
		Tout, Tout, Tc_in, Tout, 
		Tout, Tout,
		Tout, 					#Assuming the pressure valve does not change the temp
		T0, Tout, T_f]),
	"Mass flow rate (mg/s)" => mfrs,
	"Pressure (kPa)" => ustrip.([
		P_atm, P_f, P_f, P_f,
		P_atm, P_atm, P_atm, P_atm,
		P_atm, P_atm, P_f, P_atm,
		P_f, P_atm])
)
end

# Writing it to a CSV
CSV.write("stream_table.csv", stream_table);

