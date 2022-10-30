### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ 43d8e566-4a9a-11ed-0776-3db8daafb728
using SymPy, Unitful, QuadGK, Latexify

# ╔═╡ 24f118b6-7112-4aea-9d4c-eecddbb4da26
begin
@vars m_A m_B m_C H_A H_B H_C h__Hg_C h__Hg_B h__N2_B h__N2_A h__Hg_A x__N2_B x__N2_A x__Hg_A  r__Hg_B Q_out T_sat T_A T_B T_C T Cp__Hg_g Cp_N2 Δh_v__Hg T_out nonzero=true x__Hg_B


Cp_Hg_l = sympy.Function("Cp__Hg_l")(T)

end

# ╔═╡ 1fd5b132-ea10-4814-821d-8ff0818ffe14
eqs = [
	m_A ~ m_B + m_C,
	m_A*x__Hg_A ~ m_C + m_B*x__Hg_B,
	m_A*x__N2_A ~ m_B*x__N2_B,
	x__Hg_B + x__N2_B ~ 1,
	r__Hg_B ~ (m_B*x__Hg_B)/(m_A*x__Hg_A),
	H_A ~ H_B + H_C + Q_out,
	H_A ~ m_A*(x__Hg_A*h__Hg_A + x__N2_A*h__N2_A),
	H_B ~ m_B*(x__Hg_B*h__Hg_B + x__N2_B*h__N2_A),
	H_C ~ m_C*h__Hg_C,
	h__Hg_C ~ sympy.Integral(Cp_Hg_l, (T, T_sat, T_out)) + Cp__Hg_g*(T_out - T_A) - Δh_v__Hg,
	h__N2_B ~ Cp_N2*(T_out - T_A),
	h__Hg_B ~ Cp__Hg_g*(T_out - T_A),
]

# ╔═╡ f8f0e5a5-dcbf-4d5f-8ce3-f830d1b3807a
latexify.(eqs)

# ╔═╡ 4c61ec1e-88c6-491b-9b20-3d83dd703094
length(eqs)

# ╔═╡ 85fd0915-4e7f-499c-8dcb-297ba2def69d
sols = nonlinsolve(eqs, (m_B, m_C, H_B, H_C, Q_out, x_Hg_B, h_Hg_B));

# ╔═╡ 19a6296c-100a-4dec-873d-9f7d67e888b2
free_symbols(sols)

# ╔═╡ 3e197862-4a08-42ff-89d8-b22f51281bac
solutions = collect.(collect.(simplify(collect.(simplify.(collect(collect(Set(sols...))[1].args)), r_Hg_B)), Cp_Hg_g),m_A)

# ╔═╡ 2f8e0b1e-7543-44a5-884a-40c05bcc27f6
(m_B, m_C, H_B, H_C, h_Hg_C, h_Hg_B, h_N2_B, x_N2_B, x_Hg_B, Q_out)

# ╔═╡ d2fbdd96-22c0-4832-ad69-a524c00b82f0
inputs = free_symbols(collect(Set(sols...))[1])

# ╔═╡ 4b20039f-bf47-4b70-9f98-49ae08e395ee
outputs = [m_B, m_C, H_B, H_C, h_Hg_C, h_Hg_B, h_N2_B, x_N2_B, x_Hg_B, Q_out]

# ╔═╡ 6b3dff0e-7fb1-45f6-b3d3-bd8037fe66f2
explain = [outputs[i]~ solutions[i] for i=1:length(solutions)]

# ╔═╡ 64531c25-582a-46a3-af29-7bc2aa7722c7
actual_Cp(T) = (30.39 - 11.47e-3 * T)/200.59 # J/(g K)

# ╔═╡ f940c11a-1ece-4c50-a7ca-cc37cb473a72
fn = lambdify((expand.(solutions.subs(Cp_Hg_l, actual_Cp(T)))))

# ╔═╡ d2746359-19af-4206-9906-c56d395c02bd
expand.(solutions.subs(Cp_Hg_l, actual_Cp(T)))

# ╔═╡ c2a13e8e-42b2-4f7e-ac80-be0dc98fd7b0
free_symbols(collect(Set(sols...))[1].subs(sympy.Integral(Cp_Hg_l, (T, T_sat, T_out)), actual_Cp(T)))

# ╔═╡ eb3d475f-53a7-4069-9b0e-7739873b6a21
fn(ustrip.([15u"g*s^-1", 5u"J * g^-1 * K^-1", 1u"J * g^-1 * K^-1", 10u"J*s^-1", (600+273)u"k",  (10+273)u"K", 573u"K", 2u"J*g^-1", 0.1, 48u"J*g^-1"])...)

# ╔═╡ 3f5c9c7e-5df2-4f1a-a6d9-2574017a3aa3
simplify(test_int.subs(Cp_Hg_l, actual_Cp(T)))

# ╔═╡ 8c41c714-a329-4f20-9eee-1ed056237945
collect(1:10)

# ╔═╡ 6fb75877-e399-4416-888d-9ca6b5717167
convert(Expr, solutions[10])

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
Latexify = "~0.15.17"
QuadGK = "~2.5.0"
SymPy = "~1.1.7"
Unitful = "~1.12.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CommonEq]]
git-tree-sha1 = "d1beba82ceee6dc0fce8cb6b80bf600bbde66381"
uuid = "3709ef60-1bee-4518-9f2f-acd86f176c50"
version = "0.2.0"

[[deps.CommonSolve]]
git-tree-sha1 = "332a332c97c7071600984b3c31d9067e1a4e6e25"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.1"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "3ca828fe1b75fa84b021a7860bd039eaea84d2f2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.3.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6e47d11ea2776bc5627421d59cdcc1296c058071"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.7.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

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

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "6c01a9b494f6d2a9fc180a08b182fcb06f0958a0"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "53b8b07b721b77144a0fbbbc2675222ebf40a02d"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.94.1"

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

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "612a4d76ad98e9722c8ba387614539155a59e30c"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.SymPy]]
deps = ["CommonEq", "CommonSolve", "Latexify", "LinearAlgebra", "Markdown", "PyCall", "RecipesBase", "SpecialFunctions"]
git-tree-sha1 = "de83b8c89b2744fee5279326fe8e3f4a9b94d1e1"
uuid = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
version = "1.1.7"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═43d8e566-4a9a-11ed-0776-3db8daafb728
# ╠═24f118b6-7112-4aea-9d4c-eecddbb4da26
# ╠═1fd5b132-ea10-4814-821d-8ff0818ffe14
# ╠═f8f0e5a5-dcbf-4d5f-8ce3-f830d1b3807a
# ╠═4c61ec1e-88c6-491b-9b20-3d83dd703094
# ╠═85fd0915-4e7f-499c-8dcb-297ba2def69d
# ╠═19a6296c-100a-4dec-873d-9f7d67e888b2
# ╠═3e197862-4a08-42ff-89d8-b22f51281bac
# ╠═2f8e0b1e-7543-44a5-884a-40c05bcc27f6
# ╠═d2fbdd96-22c0-4832-ad69-a524c00b82f0
# ╠═4b20039f-bf47-4b70-9f98-49ae08e395ee
# ╠═6b3dff0e-7fb1-45f6-b3d3-bd8037fe66f2
# ╠═64531c25-582a-46a3-af29-7bc2aa7722c7
# ╠═f940c11a-1ece-4c50-a7ca-cc37cb473a72
# ╠═d2746359-19af-4206-9906-c56d395c02bd
# ╠═c2a13e8e-42b2-4f7e-ac80-be0dc98fd7b0
# ╠═eb3d475f-53a7-4069-9b0e-7739873b6a21
# ╠═3f5c9c7e-5df2-4f1a-a6d9-2574017a3aa3
# ╠═8c41c714-a329-4f20-9eee-1ed056237945
# ╠═6fb75877-e399-4416-888d-9ca6b5717167
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
