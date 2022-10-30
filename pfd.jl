## Building the process ##
using DataFrames

include("blocks.jl")

stream1 = strm([comp(mHg_glass + mHg_powder + mHg_Al, h0, "Hg", "s"),
                comp(mpowder, h0, "powder", "s"), comp(mglass, h0, "glass", "s"),
                comp(mAl, h0, "Al", "s")], T0)

stream2 = deepcopy(stream1)

# After sieving; glass stream
stream3 = strm([comp(mHg_glass, h0, "Hg", "s"),
                comp(mglass, h0, "glass", "s"),
                comp(mAl, h0, "Al", "s")], T0)

# After sieving; powder stream
stream4 = strm([comp(mHg_powder, h0, "Hg", "s"), comp(mpowder, h0, "powder", "s")], T0)

# Assuming the flow gas streams are the same for both furnaces
stream16 = strm([comp(m_N2, CpN2*(Tout-T0), "N2", "s")], Tout)
stream17 = deepcopy(stream16)

# Furnace for the powder
(stream5, stream19, Qf_powder) = furnace(stream4, stream17, Tf, η_fp, 0.83)

# Furnace for the glass (just using same temp as the powder furnace for now)
(stream6, stream18, Qf_glass) = furnace(stream3, stream16, Tf_glass, η_fg, 0.83)

# Stream merge balance
stream7 = mergestream(stream5, stream6)

# Condenser balance
(stream8, stream21, Qc) = condenser(stream7, Tout, η_condenser, Pf)

# Carbon filter balance
(stream9, streamCCollect) = cfilter(stream8)

stream10 = deepcopy(stream9)

# Chiller balance
(stream12, stream11, Qchiller) = chiller(Qc, CpWater, Tc, Th, η_chiller)

# These depend on the system volume
# stream14 = stream20
stream14 = deepcopy(stream10)
stream20 = deepcopy(stream10)

stream13 = deepcopy(stream10)

stream15 = deepcopy(stream13)

streams = [stream1, stream2, stream3, stream4, stream5, stream6, stream7, stream8,
          stream9, stream10, stream11, stream12, stream13, stream14, stream15, stream16,
          stream17, stream18, stream19, stream20, stream21, streamCCollect]

# Stream ID labels, in order of the streams listed above
IDs = ["1", "2", "3", "4", "5", "6", "7", "8",
          "9", "10", "11", "12", "13", "14", "15", "16", "17",
          "18", "19", "20", "21", "CCollect"]

# Updating the stream ID based on the labels above
for i in 1:length(streams)
    streams[i].ID = IDs[i]
end

# Converting stream objects to simple vectors of parameters for each component
function stream2vec(stream::strm)
    vecs = []
    for component in stream.cmp
        vec = []
        for i in 1:nfields(component)
            push!(vec, getfield(component, i))
        end
        push!(vec, stream.T)
        push!(vec, stream.ID)
        push!(vecs, vec)
    end
    return vecs
end

# Flattening the arrays
vecofvecs = collect(Iterators.flatten(stream2vec.(streams)))

# Including temperature and stream ID in the list of parameters in the table
namelist = String.(collect(fieldnames(comp)))
push!(namelist, "T")
push!(namelist, "ID")

# Pushing the parameter vectors into a dataframe
df = DataFrame()
for (i, name) in enumerate(namelist)
    df[!, name] = [vecofvecs[j][i] for j in 1:length(vecofvecs)]
end

# Creating some summary tables
equipment_power = DataFrame(variables = ["Qf_powder", "Qf_glass", "Qcondense", "Qchiller"],
                      values = [Qf_powder/tcyc_powder, Qf_glass/tcyc_glass, Qc/tcyc_powder,
                                Qchiller/tcyc_powder].|>u"W")
equipment_energy = DataFrame(variables = ["Qf_powder", "Qf_glass", "Qcondense", "Qchiller"],
                             energy = [Qf_powder, Qf_glass, Qc, Qchiller].|>u"kJ",
                             power = [Qf_powder/tcyc_powder, Qf_glass/tcyc_glass, Qc/tcyc_powder,
                                      Qchiller/tcyc_powder] .|>u"W",
                             efficiency = [η_fp, η_fg, η_condenser, η_chiller], 
                             cycle_time = [tcyc_powder, tcyc_glass, tcyc_powder, tcyc_powder])

system_parameters = DataFrame(parameter = ["Feed/batch", "T₀", "X Al in feed", "X powder in feed",
                                           "X Hg in Al", "X Hg in powder", "X Hg in glass",
                                           "η condenser", "η powder furnace", "η glass furnace",
                                           "η chiller", "Tout", "Tf glass", "Tf powder", "Tc chiller",
                                           "Th chiller", "V̇ nitrogen", "Cycle time"],
                              value = [mfeed, T0, x_Al_feed, x_feed_powder, x_Hg_Al,
                                       x_Hg_powder, x_Hg_glass, η_condenser, η_fp, η_fg,
                                       η_chiller, Tout, Tf_glass, Tf, Tc, Th, V_N2, tcyc_powder]
                             )

# Selecting, renaming, and unifying the units for each parameter
cleandf = select(df, :ID=>"Stream", :sp => "Species",
                :m=> (x-> ustrip.(x.|>u"g")) =>"Mass (g)",
                :T => ByRow(ustrip) => "Temp (K)",
                :h => (x->ustrip.(u"J*g^-1", x)) => "Enthalpy (J/g)",
                :H => (x->ustrip.(u"J", x)) => "Enthalpy (J)",
                :pz => "Phase")


# Checking mass balances
# Just for the sieve
mass_check([stream3, stream4], [stream2])
# For the whole input/output of solids (nitrogen not accounted for)
mass_check([stream2], [stream18, stream19, stream21])
# Checking the energy
# I am not including the chiller power as this is redundant with the condenser power
input_enths = [Qf_glass*η_fg, Qf_powder*η_fp, Qc*η_condenser]
energy_check([stream2, stream14],
             [stream18, stream19, stream21, streamCCollect, stream20],
             sum(input_enths))

# Calculating some output parameters
g_heat_time = streamenergy(stream18)/pmax_fg |> u"minute" # heating time for the glass furnace
p_heat_time = streamenergy(stream19)/pmax_fp |> u"minute" # heating time for the powder furnace

cuse = streamCCollect.cmp[1].m/ccap |> u"mg"              # Activated carbon use per batch

system_results = DataFrame(variable = ["Mass Hg", "Mass glass", "Mass powder", "Mass Al", "Mass N₂", "Hg Collected in filter", "Glass heating time", "Powder heating time"],
                           values = [m_Hg_total, mglass, mpowder, mAl, m_N2, cuse, g_heat_time, p_heat_time])

# Writing the dataframes to CSV files
CSV.write("cleandf.csv", cleandf)
CSV.write("system_parameters.csv", system_parameters)
CSV.write("system_results.csv", system_results)
CSV.write("equipment_energy.csv", equipment_energy)
