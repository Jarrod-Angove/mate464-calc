### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 58e6fc8e-e575-4b09-acad-c2ef0f6207cc
begin
using Unitful 		# This lets us use units in calculations
using Measurements# This lets us use uncertainties in calculations
end

# ╔═╡ 6ab18a9a-4e0c-4436-a7d2-f21d33fb1975
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
end

# ╔═╡ 86462929-5492-4dff-b220-f3029ad7a2bf
using DataFrames, LaTeXStrings, Latexify, CSV

# ╔═╡ ff57b090-8d8e-43c2-b3f0-f89ad56f57a4
using PrettyTables

# ╔═╡ 748e9a55-5034-47b3-bb8d-cc0bb948b570
using LatexPrint

# ╔═╡ 3d1ba4f5-1cc0-48fd-ae5c-b92662438220
md"""
## Establishing some variables
"""

# ╔═╡ c49232ac-0858-4961-ad50-39c106564f34
begin
# All uncertainties here are rough estimates based on variations in the literature
m_feed = (250)u"kg"; 						# Feed mass (arbitrary)
X_Hg = upreferred((3575)u"mg"/1000u"g"); 	# Mass fraction of Hg in feed powder 													(Average from Park 2016)
m_Hg = X_Hg * m_feed; 							# Total mass of mercury
M_Hg = 200.59u"g*mol^-1"; 						# Molar mass of mercury
P_f = 10.0u"kPa"; 								# Furnace pressure (Lee 2020)
R = 8.314u"J*mol^-1 * K^-1"; 	 				# Gas constant
T_f = (600+273.15)u"K"; 						# Furnace temperature (Lee 2020)
V_sys = 5u"m^3"; 								# Assumed system volume (arbitrary)
M_N2 = 28.014u"g*mol^-1"; 						# Molar mass of N₂
t_cyc = 390u"minute" |> u"s" 					# Total cycle time from Lee 2020
P_atm = 1.0u"atm"|>u"kPa" 						# Atmospheric pressure
ṁ_Hg = m_Hg / t_cyc 							# Assume constant evaporation
	
# Taking molar mass of powder as molar mass of Y₂O₃ (rough estimate based on the fact that Y₂O₃ makes up ~70% of powder mass)
M_p = 225.809u"g*mol^-1" 
m_p = m_feed*(1-X_Hg)	# Mass of the powder without the mercury

T0 = 300.0u"K" 			# Temp at which everything enters the system
Tout = (10+273.15)u"K"; # Assumed condenser output temperature is 10°C (Lee 2020)

# Volumetric flow rate will be dictated by our pump choice for the assumed furnace pressure P_f; temporarily using value from Lee 2020 as an assumption
V̇_N2 = (2)u"L*minute^-1" |> u"m^3 * s^-1";

end;

# ╔═╡ a11d24ce-8da8-485c-afbe-1adaf98464b0
Tout

# ╔═╡ ec08e747-e51c-49fb-88f2-447fa6a350b9
V_tot = m_Hg/(M_Hg * P_f) * R * T_f |> u"m^3" # Total volume of Hg gas

# ╔═╡ 8c40a519-fe74-49a6-8b62-a199d4f07940
m_Hg |>u"g" # Total mass of mercury in grams

# ╔═╡ 13283540-3e4f-4728-a8f4-0526694e1854
md"""
# Mass Balances
"""

# ╔═╡ 63692945-afaf-49bf-97a1-6d01e07e06e1
md"""
## Heat capacity stuff
"""

# ╔═╡ ad91a107-e6fc-424b-9f71-ec6d30b16ad4
# Heat capacity of liquid mercury as a function of temperature
CpHg_l(T) = 30.39u"J*mol^-1*K^-1" - 11.47e-3u"J*mol^-1*K^-2" * T ;
# From Hae-Geon Lee Materials thermodynamics pg. 431

# ╔═╡ db471787-0295-4281-b556-ede54b5d8c72
# Ideal gas heat capacity of mercury
CpHg_g = 20786u"J*kmol^-1*K^-1"

# ╔═╡ 37c5088c-d37a-412c-89a9-4419984a90f8
# Ideal gas heat capacity of nitrogen from CRC (rough approximation)
CpN2 = 30u"J*mol^-1*K^-1"

# ╔═╡ c056aea2-4ae1-4073-8ae3-d118b510d1fb
# Creating an interpolation function based on data from scientificgroupthermodataeuropesgteThermodynamicPropertiesCompounds2001
# This is for Y2O3 but we are using it as an assumption
begin

capp = [101.85,	110, 116.9, 119.2, 121.5, 123.8, 124.95, 126.1, 127.25, 128.4, 128.4]u"J * mol^-1 * K^-1";
tempp = [300, 390, 480, 570, 660, 750, 840, 930, 1020, 1110, 1200]u"K";
	
Cpp = linear_interpolation(tempp, capp);
end;

# ╔═╡ 0497e058-0f3d-4967-a5f8-7896d40d9035
md"""
## Furnace
"""

# ╔═╡ b75437d1-877c-4054-b715-ed89d8f59262
# I am assuming the pump flow rate is calculated at T=Tout as this is what the pump is pulling at inlet.
# This is assuming ideal gas (literally just the ideal gas law here)
# Should not be using P_f for this, as the pressure at the pump inlet will be lower than in the furnace, so we will change this when we do pressure loss calculations for equipment selection.
ṁ_N2 = (M_N2 * P_f * V̇_N2)/(R * Tout) |> u"g * s^-1"

# ╔═╡ 74f84d29-239d-4b08-bdd5-d860d768edba
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

# ╔═╡ 5e3179e3-247d-4d40-bf71-80823c46bcb3
begin
using QuadGK # This is a numerical integration package
# Using numerical integration over heat capacity of liquid Hg
ΔH_Hg_l = quadgk(CpHg_l, T0, Tsat, rtol=1e-8)[1] * m_Hg / M_Hg |> u"kJ"
end

# ╔═╡ 31925dc5-7381-45ba-ad8b-082c0abd123b
# Volumetric flow rate in the furnace is different becuase the temperature changes while the mass and pressure remain constant
V̇_N2_F = (ṁ_N2 * R * T_f)/(M_N2 * P_f) |> u"m^3 * s^-1"

# ╔═╡ 05acda5a-c4fc-4098-99a4-6d1cef7bfff0
md"""
### Enthalpy change of mercury
"""

# ╔═╡ 4770ad68-766a-441e-a002-df5c7cd5989a
ΔH_Hg_g = CpHg_g * (T_f-Tsat) * m_Hg / M_Hg |> u"kJ"

# ╔═╡ 73a6edec-abf8-4067-afda-598f5de25f86
# Heat of vaporization of mercury from CRC handbook
# This is at 1atm pressure, so it may not be accurate
ΔH_v_Hg = 59.11u"kJ*mol^-1" * m_Hg / M_Hg |> u"kJ"

# ╔═╡ bbe04a52-bdc9-4559-bc1d-743e48f92cb3
h_v_Hg = 59.11u"kJ*mol^-1" / M_Hg |> u"J*g^-1"

# ╔═╡ 21cb1220-029d-49df-9fcc-ed93cd4cd269
ΔH_Hg_F = ΔH_Hg_l + ΔH_Hg_g + ΔH_v_Hg

# ╔═╡ 669fc0e7-1e36-4f15-bcce-fc77c2d00405
md"""
### Enthalpy change of powder
"""

# ╔═╡ 729c8de5-dacb-473c-bd88-e26a9a047d56
Cpp

# ╔═╡ 410abdf4-06bf-4394-9498-75d25071e262
# numerically integrating over the interpolation function created above
# Assuming there are no phase changes
ΔH_powder = quadgk(Cpp, T0, T_f)[1] * m_p / M_p |> u"kJ"
# Note that the poweder leaves at the same temperature, so this enthalpy is not relavent to the steady state system

# ╔═╡ bab13792-36b3-48cb-83c9-0ffd1e956214
md"""
### Enthalpy change of nitrogen
"""

# ╔═╡ 3c7bc9fc-686e-4437-97ef-93b71953c58a
ΔH_N2_F = m_N2 * CpN2 * (T_f - Tout) / M_N2 |> u"kJ"

# ╔═╡ 410016f1-6697-4f99-a56c-09f7acd66f9e
md"""
### Enthalpy of the furnace
"""

# ╔═╡ 536f4265-c9f1-4247-b53b-68aaaf99eafe
# Total energy input of the furnace for a single batch
ΔH_F = ΔH_powder + ΔH_N2_F + ΔH_Hg_F

# ╔═╡ 0655f266-7107-4bb3-9860-4a5a4a9b173c
# Power required during continuous portion of the process
ΔḢ_F_steady = (ΔH_N2_F + ΔH_Hg_F)/t_cyc |> u"W"

# ╔═╡ 5b8249b7-dd85-44bc-ab9b-3704034e6065
md"""
This is an extremely small furnace, even if we account for efficiency losses!
"""

# ╔═╡ 6f7d518a-8d6a-4752-943e-255c4959dae3
md"""
## Condenser Balance
"""

# ╔═╡ d6321fbe-bb95-42d6-b2ff-2277dd3b01bc
md"""
Since everything else is known or assumed, this will be a calculation of our coolant (water) mass flow rate for a given temperature
"""

# ╔═╡ 6e8c20c8-f7c6-40b4-948b-ff700fbf29f6
# Assuming the coolant enters the system at 5C
Tc_in = (5+273.15)u"K"

# ╔═╡ 3fbb9964-9195-4884-8de4-4be527471326
ΔH_N2_c = CpN2 * m_N2 * (Tout - T_f) / M_N2 |> u"kJ" 
# Input temp is the temp of the furnace T_f

# ╔═╡ 6658372e-2e59-4da1-97fd-9b865e3ef890
ΔH_Hg_c = (-ΔH_Hg_g - ΔH_v_Hg + m_Hg * quadgk(CpHg_l, Tsat, Tout)[1]/M_Hg)  |> u"kJ" 

# ╔═╡ a7d7c699-99c2-49fa-bb83-f5163c58fe13
Δh_c = quadgk(Cpc_l, Tout, Tc_in)[1]

# ╔═╡ 9fdb5e93-d631-42e6-b700-f00a6fdd08d7
# Mass flow rate of coolant required for given parameters
ṁ_c = (ΔH_N2_c + ΔH_Hg_c)/Δh_c 

# ╔═╡ bb8057a4-7376-4c65-b5a9-e4b61244a95e
ΔH_c = ṁ_c*quadgk(Cpc_l, Tout, Tc_in)[1]

# ╔═╡ cedbabde-5a6a-4438-bcea-a0a77ed90f7d
# Checking if the condenser energy balance closes
ΔH_N2_c + ΔH_Hg_c ≈ ΔH_c

# ╔═╡ 68970f40-d1e6-4b70-8923-f81ae4ffe21d
md"""
Therefore, the mass flow rate of the coolant is $0.34 \frac{g}{s}$
"""

# ╔═╡ dabbd1ec-0da7-4670-a7c5-c3c8d978c930
md"""
Alternatively, we can find the condenser power:
"""

# ╔═╡ 3a76a00a-050a-457c-99d4-46b770214a39
ΔH_chiller = -(ΔH_N2_c + ΔH_Hg_c) |> u"kJ" 
# This is the power  required by the chiller unit before efficiency losses

# ╔═╡ 3392f0f9-c2db-4df0-8c06-fcfd1cbaad8a
md"""
 $\therefore$ power requirements of the condensor are trivial
"""

# ╔═╡ e8e5019d-9512-47be-8c5c-16f782dd4a40
md"""
## Pump stuff
"""

# ╔═╡ 59708a48-5501-4c9f-9be0-f0825adda7b7
# Rough estimate of the pressure difference generated by the pump
# This will change if we consider pressure losses during equipment selection
ΔP_p = P_atm - P_f |> u"kPa"

# ╔═╡ c287482d-6abc-4854-973a-e00c882bff67
# Work required by pump before accounting for efficiency
W_p = V̇_N2 * ΔP_p |> u"W"

# ╔═╡ f2f6c0af-5e46-498c-b8e9-5badf70f9a29
# Enthalpy change from pump (i.e) ΔH(s₄ --> s₅)
ΔH_p = W_p * t_cyc |> u"kJ"

# ╔═╡ 4af8d03d-b89a-40d1-898c-9432606c84f5
md"""
## System evacuation and batch features
"""

# ╔═╡ 658dc976-b526-4bed-b9b0-60cdb29472b6
md"""
At the end of each cycle, all gas must be evacuated from the system before the powder can be removed. This is based on the total system volume, which is not yet known. However, it is known that the mass of flue gas is equivalent to the mass of the nitrogen that must be injected into the system before the next cycle. For this, we will assume that half of the system is at the temperature of the furnace, and the other half is at the temperature of the output.
"""

# ╔═╡ e84d48a8-0601-4456-a378-39e0a3c2d9cc
md"""
## Table creation
"""

# ╔═╡ 7dea7e8f-13c6-4acc-a22f-b1cb0244585a
md"""
For all streams:
- Temperature
- Species; Phase
- Mass flow rate
"""

# ╔═╡ dde1e743-4c51-45ab-ac16-27f9ce09d57a
enths = DataFrame(
"Name" =>
["\\Delta H_{\\mathrm{Hg,C}}",
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

# ╔═╡ c594262a-76f9-4730-b9c5-ac3716daaefd
CSV.write("enthalpy.csv", enths)

# ╔═╡ ade0ca4a-3d6e-459a-8ec3-a53565adc979
typeof(enths)

# ╔═╡ 02c405b8-3ed3-4ae5-8f10-aafbc79450bd
latexify(L"$ \Delta H_{\mathrm{Hg,C}} $")

# ╔═╡ 7f311fd7-e550-4a63-9b41-1f56228ccb6a
pretty_table(enths)

# ╔═╡ d8b01ee6-5de3-437e-a119-64b2e132c27f
enths.Name[1] |> typeof

# ╔═╡ 160ede70-9ec8-4fa1-a8a6-2e4409a1eb12
print(latexify.([L"\Delta H_{\mathrm{Hg,C}}",
L"\Delta H_{\mathrm{N_2, C}}",
L"\Delta H_\mathrm{Hg, F}",
L"\Delta H_{\mathrm{N_2, F}}",
L"\Delta H_{\mathrm{p}}",
L"\Delta H_\mathrm{C}", 
L"\Delta H_\mathrm{F}" ]))

# ╔═╡ b222d890-e01e-42c2-951c-3b2405640daa
# ╠═╡ show_logs = false
lap(enths)

# ╔═╡ 12d70b56-df27-4db4-a222-13dbcee3b249
latexify(L"\left[
\begin{array}{c}
\text{$\Delta H_{\mathrm{Hg,C}}$} \\
\text{$\Delta H_{\mathrm{N_2, C}}$} \\
\text{$\Delta H_\mathrm{Hg, F}$} \\
\text{$\Delta H_{\mathrm{N_2, F}}$} \\
\text{$\Delta H_{\mathrm{p}}$} \\
\text{$\Delta H_\mathrm{C}$} \\
\text{$\Delta H_\mathrm{F}$} \\
\end{array}
\right]")

# ╔═╡ 252d92f4-02c6-45e2-8584-8c2c7bc06a4a
lap("\Delta H_{\mathrm{Hg,C}")

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LatexPrint = "d2208f48-c256-5759-9089-c25ed2a93924"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
CSV = "~0.10.4"
DataFrames = "~1.4.1"
Interpolations = "~0.14.6"
LaTeXStrings = "~1.3.0"
LatexPrint = "~1.1.0"
Latexify = "~0.15.17"
Measurements = "~2.8.0"
PrettyTables = "~2.1.2"
QuadGK = "~2.5.0"
Unitful = "~1.12.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "da4799f2e9b1d4d34fac2c0c6355c46e51a03418"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "3ca828fe1b75fa84b021a7860bd039eaea84d2f2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.3.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "46d2680e618f8abd007bce0c3026cb0c4a8f2032"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.12.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "558078b0b78278683a7445c626ee78c86b9bb000"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "db619c421554e1e7e07491b85a8f4b96b3f04ca0"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "842dd89a6cb75e02e85fdd75c760cdc43f5d6863"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.6"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LatexPrint]]
deps = ["Requires"]
git-tree-sha1 = "67019e47920747446c28c9ff4eadef0306c14031"
uuid = "d2208f48-c256-5759-9089-c25ed2a93924"
version = "1.1.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "12950d646ce04fb2e89ba5bd890205882c3592d7"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.8.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "f71d8950b724e9ff6110fc948dff5a329f901d64"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.8"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "6c01a9b494f6d2a9fc180a08b182fcb06f0958a0"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "460d9e154365e058c4d886f6f7d6df5ffa1ea80e"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.1.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "3c009334f45dfd546a16a57960a821a1a023d241"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.5.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "d12e612bba40d189cead6ff857ddb67bd2e6a387"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "efd23b378ea5f2db53a55ae53d3133de4e080aa9"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.16"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "f86b3a049e5d05227b10e15dbb315c5b90f14988"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.9"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d57a4ed70b6f9ff1da6719f5f2713706d57e0d66"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.12.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"
"""

# ╔═╡ Cell order:
# ╟─3d1ba4f5-1cc0-48fd-ae5c-b92662438220
# ╠═58e6fc8e-e575-4b09-acad-c2ef0f6207cc
# ╠═c49232ac-0858-4961-ad50-39c106564f34
# ╠═a11d24ce-8da8-485c-afbe-1adaf98464b0
# ╠═ec08e747-e51c-49fb-88f2-447fa6a350b9
# ╠═8c40a519-fe74-49a6-8b62-a199d4f07940
# ╟─13283540-3e4f-4728-a8f4-0526694e1854
# ╟─63692945-afaf-49bf-97a1-6d01e07e06e1
# ╠═ad91a107-e6fc-424b-9f71-ec6d30b16ad4
# ╠═6ab18a9a-4e0c-4436-a7d2-f21d33fb1975
# ╠═db471787-0295-4281-b556-ede54b5d8c72
# ╠═37c5088c-d37a-412c-89a9-4419984a90f8
# ╠═c056aea2-4ae1-4073-8ae3-d118b510d1fb
# ╟─0497e058-0f3d-4967-a5f8-7896d40d9035
# ╠═b75437d1-877c-4054-b715-ed89d8f59262
# ╠═74f84d29-239d-4b08-bdd5-d860d768edba
# ╠═31925dc5-7381-45ba-ad8b-082c0abd123b
# ╟─05acda5a-c4fc-4098-99a4-6d1cef7bfff0
# ╠═5e3179e3-247d-4d40-bf71-80823c46bcb3
# ╠═4770ad68-766a-441e-a002-df5c7cd5989a
# ╠═73a6edec-abf8-4067-afda-598f5de25f86
# ╠═bbe04a52-bdc9-4559-bc1d-743e48f92cb3
# ╠═21cb1220-029d-49df-9fcc-ed93cd4cd269
# ╟─669fc0e7-1e36-4f15-bcce-fc77c2d00405
# ╠═729c8de5-dacb-473c-bd88-e26a9a047d56
# ╠═410abdf4-06bf-4394-9498-75d25071e262
# ╟─bab13792-36b3-48cb-83c9-0ffd1e956214
# ╠═3c7bc9fc-686e-4437-97ef-93b71953c58a
# ╟─410016f1-6697-4f99-a56c-09f7acd66f9e
# ╠═536f4265-c9f1-4247-b53b-68aaaf99eafe
# ╠═0655f266-7107-4bb3-9860-4a5a4a9b173c
# ╟─5b8249b7-dd85-44bc-ab9b-3704034e6065
# ╟─6f7d518a-8d6a-4752-943e-255c4959dae3
# ╟─d6321fbe-bb95-42d6-b2ff-2277dd3b01bc
# ╠═6e8c20c8-f7c6-40b4-948b-ff700fbf29f6
# ╠═3fbb9964-9195-4884-8de4-4be527471326
# ╠═6658372e-2e59-4da1-97fd-9b865e3ef890
# ╠═a7d7c699-99c2-49fa-bb83-f5163c58fe13
# ╠═9fdb5e93-d631-42e6-b700-f00a6fdd08d7
# ╠═bb8057a4-7376-4c65-b5a9-e4b61244a95e
# ╠═cedbabde-5a6a-4438-bcea-a0a77ed90f7d
# ╟─68970f40-d1e6-4b70-8923-f81ae4ffe21d
# ╟─dabbd1ec-0da7-4670-a7c5-c3c8d978c930
# ╠═3a76a00a-050a-457c-99d4-46b770214a39
# ╟─3392f0f9-c2db-4df0-8c06-fcfd1cbaad8a
# ╟─e8e5019d-9512-47be-8c5c-16f782dd4a40
# ╠═59708a48-5501-4c9f-9be0-f0825adda7b7
# ╠═c287482d-6abc-4854-973a-e00c882bff67
# ╠═f2f6c0af-5e46-498c-b8e9-5badf70f9a29
# ╟─4af8d03d-b89a-40d1-898c-9432606c84f5
# ╟─658dc976-b526-4bed-b9b0-60cdb29472b6
# ╟─e84d48a8-0601-4456-a378-39e0a3c2d9cc
# ╟─7dea7e8f-13c6-4acc-a22f-b1cb0244585a
# ╠═86462929-5492-4dff-b220-f3029ad7a2bf
# ╠═dde1e743-4c51-45ab-ac16-27f9ce09d57a
# ╠═c594262a-76f9-4730-b9c5-ac3716daaefd
# ╠═ade0ca4a-3d6e-459a-8ec3-a53565adc979
# ╠═02c405b8-3ed3-4ae5-8f10-aafbc79450bd
# ╠═ff57b090-8d8e-43c2-b3f0-f89ad56f57a4
# ╠═7f311fd7-e550-4a63-9b41-1f56228ccb6a
# ╠═d8b01ee6-5de3-437e-a119-64b2e132c27f
# ╠═160ede70-9ec8-4fa1-a8a6-2e4409a1eb12
# ╠═748e9a55-5034-47b3-bb8d-cc0bb948b570
# ╠═b222d890-e01e-42c2-951c-3b2405640daa
# ╠═12d70b56-df27-4db4-a222-13dbcee3b249
# ╠═252d92f4-02c6-45e2-8584-8c2c7bc06a4a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
