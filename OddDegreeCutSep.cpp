//
// C++ implementation of an exact separation routine for the (lifted) odd degree cut set inequalities
//
// by Artur Pessoa (2019)
//
#include "OddDegreeCutSep.h"

#include <queue>
#include <iostream>

#define PRINT_VIOLATIONS

const double ValueScale = 1e6;
const int Infinity = 1000000000;

Graph::Graph(int n, int m, int r, int* ps, int* pt, double* pv,
	  int* ri, int* rj)
{
	_n = n;

	// add the arcs to the graph
	adj.resize(n);
	val.resize(n);
	opp.resize(n);
	for (int a = 0; a < m; a++)
	{
		adj[ps[a]].push_back(pt[a]);
		val[ps[a]].push_back(int(ValueScale * pv[a]));
		opp[ps[a]].push_back(opp[pt[a]].size());
		adj[pt[a]].push_back(ps[a]);
		val[pt[a]].push_back(int(ValueScale * pv[a]));
		opp[pt[a]].push_back(opp[ps[a]].size() - 1);
	}

	// set the vertex marks for those with odd number of adjacent required edges
	std::vector<int> deg(n, 0);
	for (int e = 0; e < r; e++)
	{
		deg[ri[e]]++;
		deg[rj[e]]++;
	}
	mark.resize(n);
	for (int v = 0; v < n; v++)
		mark[v] = ((deg[v] % 2) == 1);
}

int Graph::findPath(int s, int t, std::vector<int>& p)
{
	// do a bfs in the graph
	std::vector<int> prev(_n, -2);
	std::queue<int> q;
	q.push(s);
	prev[s] = -1;
	while (!q.empty())
	{
		int i = q.front();
		q.pop();

		// check if the sink was reached
		if (i == t)
		{
			// fill the path and return the minimum value
			int minVal = Infinity;
			p.clear();
			while (prev[i] != -1)
			{
				p.push_back(prev[i]);
				int k = opp[i][prev[i]];
				i = adj[i][prev[i]];
				if (minVal > val[i][k])
					minVal = val[i][k];
			}
			return minVal;
		}

		// check the adjacent vertices
		for (int k = 0; k < (int)adj[i].size(); k++)
		{
			if (val[i][k] == 0) continue;
			int j = adj[i][k];
			if (prev[j] != -2) continue;
			prev[j] = opp[i][k];
			q.push(j);
		}
	}

	// store in p a list of vertices that are reachable from s and return 0
	p.clear();
	for (int i = 0; i < _n; i++)
		if (prev[i] != -2)
			p.push_back(i);
	return 0;
}

//===================================================================================

OddDegreeCutSep::OddDegreeCutSep(int n, int m, int r, int* ps, int* pt, double* pv,
		int* ri, int* rj) :
		instance(n, m, r, ps, pt, pv, ri, rj)
{
}

int OddDegreeCutSep::separate(int max_cuts, std::vector< std::vector<int> >& cutSets,
		double minViol)
{
	// build a list of marked vertices
	std::vector<int> mv;
	for (int i = 0; i < instance._n; i++)
		if (instance.mark[i]) mv.push_back(i);

	// return if there are no marked vertices
	cutSets.clear();
	assert((mv.size() % 2) == 0);	// need an even number of odd degrees in a graph
	if (mv.empty()) return 0;
	int ncuts = 0;

	// add the first minimum cut between the first two marked vertices
	std::vector< std::pair<std::pair<int, int>, int> > treeEdges;
	std::vector< std::vector<bool> > minCuts;
	std::vector<bool> currCut;
	int cutVal = findMinCut(mv[0], mv[1], currCut);
	minCuts.push_back(currCut);
	treeEdges.push_back(std::make_pair(std::make_pair(mv[0], mv[1]), cutVal));

	// iterate adding one remaining vertex at a time
	for (int l = 2; l < (int)mv.size(); l++)
	{
		// build a list of candidates to connect to mv[l]
		std::vector<int> cand;
		for (int k = 0; k < l; k++) cand.push_back(mv[k]);

		// find the tree vertex to connect to mv[l]
		while (cand.size() > 1)
		{
			// find the smallest edge value connecting the candidates
			int minVal = Infinity;
			int best_q = -1;
			for (int q = 0; q < (int)treeEdges.size(); q++)
			{
				if (!inSet(treeEdges[q].first.first, cand)) continue;
				if (!inSet(treeEdges[q].first.second, cand)) continue;
				if (minVal > treeEdges[q].second)
				{
					minVal = treeEdges[q].second;
					best_q = q;
				}
			}

			// filter the candidate by the vertices in the same side of the cut as mv(l)
			intersect(cand, minCuts[best_q], !minCuts[best_q][mv[l]]);
		}
		assert(!cand.empty());

		// add the first minimum cut between the mv[l] and cand[0]
		cutVal = findMinCut(cand[0], mv[l], currCut);
		minCuts.push_back(currCut);
		treeEdges.push_back(std::make_pair(std::make_pair(cand[0], mv[l]), cutVal));

		// check if the current cut is valid and violated
		int count = 0;
		for (int k = 0; k < (int)mv.size(); k++)
			if (currCut[mv[k]]) count++;
		if ((count % 2) == 1)
		{
			double viol = 1.0 - (double(cutVal) / ValueScale);
			if (viol >= minViol)
			{
#ifdef PRINT_VIOLATIONS
				std::cout << "Found cut with violation " << viol << "!" << std::endl;
#endif
				// add the cut
				cutSets.push_back(std::vector<int>());
				for (int i = 0; i < instance._n; i++)
					if (currCut[i]) cutSets[ncuts].push_back(i);
				ncuts++;

				// stop if reached the maximum number of separated cuts
				if (ncuts == max_cuts) return max_cuts;
			}
		}
	}
	return ncuts;
}

int OddDegreeCutSep::findMinCut(int s, int t, std::vector<bool>& cut)
{
	// make a copy of the instance graph as the residual graph
	Graph residual = instance;

	// iterate finding paths from s to t
	std::vector<int> p;
	int flow;
	int maxFlow = 0;
	while ((flow = residual.findPath(s, t, p)) > 0)
	{
		// pass flow through the path
		int i = t;
		int l = 0;
		while (l < (int)p.size())
		{
			int ki = p[l];
			int j = residual.adj[i][ki];
			int kj = residual.opp[i][ki];
			residual.val[j][kj] -= flow;
			residual.val[i][ki] += flow;
			i = j;
			l++;
		}
		assert(i == s);
		maxFlow += flow;
	}

	// here, p contains a list of vertices in the cut... convert it the boolean vector cut
	cut.clear();
	cut.resize(residual._n, false);
	for (int l = 0; l < (int)p.size(); l++)
		cut[p[l]] = true;
	return maxFlow;
}

void OddDegreeCutSep::intersect(std::vector<int>& A, std::vector<bool>& B, bool complB)
{
	int l = 0;
	for (int k = 0; k < (int)A.size(); k++)
		if (B[A[k]] != complB)
			A[l++] = A[k];
	A.resize(l);
}

bool OddDegreeCutSep::inSet(int a, std::vector<int> A)
{
	if (A.empty()) return false;
	if (A[0] == a) return true;
	if (A[0] > a) return false;

	// do a binary search
	int l, r, m;
	l = 0; r = (int)A.size();
	while (r > (l + 1))
	{
		m = (l + r) / 2;
		if (A[m] == a) return true;
		if (A[m] > a)
			r = m;
		else
			l = m;
	}
	return false;
}

//===================================================================================

int odcs_separate(int n, int m, int r, int* ps, int* pt, double* pv,
		int* ri, int* rj, int max_cuts, int* cut_sets, double min_viol)
{
	// construct a separator and call the separation routine
	OddDegreeCutSep sep(n, m, r, ps, pt, pv, ri, rj);
	std::vector< std::vector<int> > cutSets;
	int ncuts = sep.separate(max_cuts, cutSets, min_viol);

	// convert the cuts and return
	int l = 0;
	for (int c = 0; c < ncuts; c++)
	{
		for (int k = 0; k < (int) cutSets[c].size(); k++)
			cut_sets[l++] = cutSets[c][k];
		cut_sets[l++] = -1;
	}
	return ncuts;
}
