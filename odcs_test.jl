#
# Test routines for the exact separation of (lifted) odd degree cut set inequalities
#
# by Artur Pessoa (2019)
#

include("OddDegreeCutSep.jl")

const max_cuts = 50

# instance data
ps = Vector{Int}()
pt = Vector{Int}()
pv = Vector{Float64}()
ri = Vector{Int}()
rj = Vector{Int}()

# check the number of program arguments
if (length(ARGS) != 1)
	println("Use: OddDegreeCutSep <filename>")
	exit(-1)
end

# read all the input file to the vector "numbers"
f = open(ARGS[1], "r")
numbers = map(y->parse(Float64,y),split(read(f,String)))
close(f)

# get the data from "numbers"
n = Int(numbers[1])
m = Int(numbers[2])
r = Int(numbers[3])
pos = 4
@show numbers
@show n, m, r
for a = 1:m
	push!(ps, Int(numbers[pos]) + 1)
	push!(pt, Int(numbers[pos + 1]) + 1)
	push!(pv, numbers[pos + 2])
	pos += 3
end
for e = 1:r
	push!(ri, Int(numbers[pos]) + 1)
	push!(rj, Int(numbers[pos + 1]) + 1)
	pos += 2
end

# call the separation function
ncuts, cutSets = odcs_separate(n, m, r, ps, pt, pv, ri, rj, max_cuts, 0.001)

# print the cuts found
if ncuts == 0
	println("No cuts found")
else
	println("Found $ncuts cuts!")
	for c = 1:ncuts
		print("$c", ":")
		ll = length(cutSets[c])
		for l = 1:ll
			print(" ", cutSets[c][l])
		end
		println("")
	end
end

