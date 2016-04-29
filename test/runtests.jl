
include("Tests.jl")
include("PTests.jl")

using Tests: tests
using PTests: ptests

tests()
ptests()
