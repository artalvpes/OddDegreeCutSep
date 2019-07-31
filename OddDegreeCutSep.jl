#
# Julia implementation of an exact separation routine for the (lifted) odd degree cut set inequalities
#
# by Artur Pessoa (2019)
#

using DataStructures

const PRINT_LEVEL = 1

level(l::Int) = (l <= PRINT_LEVEL) 

const ValueScale = 1e6
const Infinity = 1000000000

mutable struct Graph
	_n::Int
	adj::Vector{Vector{Int}}
	val::Vector{Vector{Int}}
	mark::Vector{Bool}

	# position of the opposite arc in the adjacency list of its head node
	opp::Vector{Vector{Int}}
end

# constructor
function Graph(n::Int, m::Int, r::Int, ps::Vector{Int}, pt::Vector{Int}, pv::Vector{Float64},
				ri::Vector{Int}, rj::Vector{Int})
	g = Graph(n, [Vector{Int}() for i=1:n], [Vector{Int}() for i=1:n], zeros(Bool, n), [Vector{Int}() for i=1:n])

	# add the arcs to the graph
	for a = 1:m
		push!(g.adj[ps[a]], pt[a])
		push!(g.val[ps[a]], Int(trunc(ValueScale * pv[a])))
		push!(g.opp[ps[a]], length(g.opp[pt[a]])+1)
		push!(g.adj[pt[a]], ps[a])
		push!(g.val[pt[a]], Int(trunc(ValueScale * pv[a])))
		push!(g.opp[pt[a]], length(g.opp[ps[a]]))
    end

	# set the vertex marks for those with odd number of adjacent required edges
	deg = zeros(Int, n)
	for e = 1:r
		deg[ri[e]] += 1
		deg[rj[e]] += 1
	end
	resize!(g.mark, n)
	for v = 1:n
		g.mark[v] = ((deg[v] % 2) == 1)
	end
	
	return g
end

# find a path p from s to t using only non-zero arcs and return its minimum arc value, and p
# * the path is represented by a sequence of vertex indices in adjacency lists starting at
#   t until reach s.
# * if no such a path exist, store in p a list of vertices that are reachable from s and
#   return 0, and p.
function findPath(g::Graph, s::Int, t::Int)::Tuple{Int, Vector{Int}}
	# do a bfs in the graph
	prev = fill(-2, g._n)
	q = CircularDeque{Int}(g._n)
	push!(q, s)
	prev[s] = -1
	while !isempty(q)
		i = shift!(q)

		# check if the sink was reached
		if i == t
			# fill the path and return the minimum value
			minVal = Infinity
			p = Vector{Int}()
			while prev[i] != -1
				push!(p, prev[i])
				k = g.opp[i][prev[i]]
				i = g.adj[i][prev[i]]
				if (minVal > g.val[i][k])
					minVal = g.val[i][k]
				end
			end
			return minVal, p
		end

		# check the adjacent vertices
		l = length(g.adj[i])
		for k = 1:l
			if g.val[i][k] == 0
				continue
			end
			j = g.adj[i][k]
			if prev[j] != -2
				continue
			end
			prev[j] = g.opp[i][k]
			push!(q, j)
		end
	end

	# store in p a list of vertices that are reachable from s and return 0
	p = Vector{Int}()
	for i = 1:g._n
		if prev[i] != -2
			push!(p, i)
		end
	end
	return 0, p
end

# find a minimum s-t cut in the graph instance and return the cut value
function findMinCut(instance::Graph, s::Int, t::Int)::Tuple{Int, Vector{Bool}}
	# make a copy of the instance graph as the residual graph
	residual = deepcopy(instance)

	# iterate finding paths from s to t
	p = Vector{Int}()
	flow, p = findPath(residual, s, t)
	maxFlow = 0
	while (flow > 0)
		# pass flow through the path
		i = t
		ll = length(p)
		for l = 1:ll
			ki = p[l]
			j = residual.adj[i][ki]
			kj = residual.opp[i][ki]
			residual.val[j][kj] -= flow
			residual.val[i][ki] += flow
			i = j
		end
		@assert i == s
		maxFlow += flow
		flow, p = findPath(residual, s, t)
	end

	# here, p contains a list of vertices in the cut... convert it the boolean vector cut
	cut = zeros(Bool, residual._n)
	ll = length(p)
	for l = 1:ll
		cut[p[l]] = true
	end
	return maxFlow, cut
end

# calculate the intersection between sets A and B of integers where B is represented as a
# boolean vector and save the resulting set in A
# @param complB indicate that B is complemented when true
function intersect!(A::Vector{Int}, B::Vector{Bool}, complB::Bool)
	l = 0
	kk = length(A)
	for k = 1:kk
		if B[A[k]] != complB
			l += 1
			A[l] = A[k]
		end
	end
	resize!(A, l)
	level(3) && @show A
end

# check if a belongs to set A represented as a sorted list
function inSet(a::Int, A::Vector{Int})::Bool
	if isempty(A)
		return false
	end
	if A[1] == a
		return true
	end
	if A[1] > a
		return false
	end

	# do a binary search
	l = 1
	r = length(A) + 1
	while r > (l + 1)
		m = div(l + r, 2)
		if A[m] == a
			return true
		end
		if A[m] > a
			r = m
		else
			l = m
		end
	end
	return false;
end

# main separation routine
# @param instance separation problem instance
# @param max_cuts maximum number of cuts to be separated
# @param minViol minimum violation required to add a new cut
# @return number of separated cuts, list of vertices indices for each separated cut.
function odcs_separate(instance::Graph, maxCuts::Int, minViol::Float64)::Tuple{Int, Vector{Vector{Int}}}
	# build a list of marked vertices
	mv = Vector{Int}()
	for i = 1:instance._n
		if instance.mark[i]
			push!(mv, i)
		end
	end

	# return if there are no marked vertices
	cutSets = Vector{Vector{Int}}()
	@assert (length(mv) % 2) == 0	# need an even number of odd degrees in a graph
	if isempty(mv)
		return 0
	end
	ncuts = 0
	
	level(2) && @show instance

	# add the first minimum cut between the first two marked vertices
	treeEdges = Vector{Tuple{Int, Int, Int}}()
	minCuts = Vector{Vector{Bool}}()
	cutVal, currCut = findMinCut(instance, mv[1], mv[2])
	push!(minCuts, currCut)
	push!(treeEdges, (mv[1], mv[2], cutVal))

	level(2) && @show mv
	level(2) && @show cutVal
	level(2) && @show currCut

	# iterate adding one remaining vertex at a time
	ll = length(mv)
	for l = 3:ll
		# build a list of candidates to connect to mv[l]
		cand = Vector{Int}()
		for k = 1:(l - 1)
			push!(cand, mv[k])
		end

		level(2) && @show l
		level(2) && @show cand
		level(2) && @show treeEdges

		# find the tree vertex to connect to mv[l]
		while length(cand) > 1
			# find the smallest edge value connecting the candidates
			minVal = Infinity
			best_q = -1
			qq = length(treeEdges)
			for q = 1:qq
				if !inSet(treeEdges[q][1], cand) || !inSet(treeEdges[q][2], cand)
					continue
				end
				if minVal > treeEdges[q][3]
					minVal = treeEdges[q][3]
					best_q = q
				end
			end

			level(2) && @show minVal
			level(2) && @show best_q

			# filter the candidate by the vertices in the same side of the cut as mv(l)
			intersect!(cand, minCuts[best_q], !minCuts[best_q][mv[l]])
			
			level(2) && @show cand
		end
		@assert !isempty(cand)

		# add the first minimum cut between the mv[l] and cand[1]
		cutVal, currCut = findMinCut(instance, cand[1], mv[l])
		push!(minCuts, currCut)
		push!(treeEdges, (cand[1], mv[l], cutVal))

		level(2) && @show cutVal
		level(2) && @show currCut
		
		# check if the current cut is valid and violated
		count = 0
		kk = length(mv)
		for k = 1:kk
			if currCut[mv[k]]
				count += 1
			end
		end
		if (count % 2) == 1
			viol = 1.0 - (Float64(cutVal) / ValueScale)
			if viol >= minViol
				level(1) && println("Found cut with violation $viol", "!")
				
				# add the cut
				push!(cutSets, Vector{Int}());
				ncuts += 1
				for i = 1:instance._n
					if currCut[i]
						push!(cutSets[ncuts], i)
					end
				end

				# stop if reached the maximum number of separated cuts
				if ncuts == maxCuts
					return maxCuts
				end
			end
		end
	end
	return ncuts, cutSets
end
  
# function that separates odd degree cutset cuts via construction of separator trees (according to the
# book Network Flows: theory, algorithms, and applications by Ahuja et al)
# @param n number of vertices
# @param m number of (shortest path) arcs with non-zero values
# @param r number of required edges
# @param ps array of start vertices by (shortest path) arc index
# @param pt array of end vertices by (shortest path) arc index
# @param pv array of arc values by (shortest path) arc index
# @param ri array of first vertices by required edge index
# @param ri array of second vertices by required edge index
# @param max_cuts maximum number of cuts to be separated
# @param min_viol minimum violation required to add a new cut
# @return number of separated cuts, list of vertices indices for each separated cut.
function odcs_separate(n::Int, m::Int, r::Int, ps::Vector{Int}, pt::Vector{Int}, pv::Vector{Float64},
			ri::Vector{Int}, rj::Vector{Int}, max_cuts::Int, min_viol::Float64)::Tuple{Int, Vector{Vector{Int}}}

	# construct a graph and call the separation routine
	g = Graph(n, m, r, ps, pt, pv, ri, rj);
	return odcs_separate(g, max_cuts, min_viol);
end
