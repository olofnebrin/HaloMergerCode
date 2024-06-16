include("./MergerTree.jl")
using .MergerTree

z0         = 0.0
M0         = 1e13
z_max      = 30.0
M_res      = 1e7

#HostMassRedshiftVvir   = "HostMassRedshiftVvir.txt"
#SubMassSubRedshift     = "SubMassSubRedshift_M0_MW.txt"
#StarFormingHalos       = "StarFormingHalos_MW1.1e12_Mres7e5_29Dec2021.txt"

HostMassRedshiftVvir   = "HostMassRedshiftVvir_13.txt"
SubMassSubRedshift     = "SubMassSubRedshift_M013.txt"

@time M_main( M0, z0, z_max, M_res,
              "HostMassRedshiftVvir_M01e13.txt", SubMassSubRedshift, true )

#@time FindStarFormingHalos( M0, z0, z_max, M_res, StarFormingHalos, true )
