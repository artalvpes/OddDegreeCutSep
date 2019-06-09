//
// C++ implementation of an exact separation routine for the (lifted) odd degree cut set inequalities
//
// by Artur Pessoa (2019)
//
#include <vector>

struct Graph
{
	int _n;
	std::vector< std::vector<int> > adj;
	std::vector< std::vector<int> > val;
	std::vector<bool> mark;

	// position of the opposite arc in the adjacency list of its head node
	std::vector< std::vector<int> > opp;

	// constructor
	Graph(int n, int m, int r, int* ps, int* pt, double* pv,
		  int* ri, int* rj);

	// find a path p from s to t using only non-zero arcs and return its minimum arc value
	// * the path is represented by a sequence of vertex indices in adjacency lists starting at
	//   t until reach s.
	// * if no such a path exist, store in p a list of vertices that are reachable from s and
	//   return 0.
	int findPath(int s, int t, std::vector<int>& p);
};

class OddDegreeCutSep
{
public:
	// constructor
	OddDegreeCutSep(int n, int m, int r, int* ps, int* pt, double* pv,
					int* ri, int* rj);

	// separation routine
	// @param max_cuts maximum number of cuts to be separated
	// @param cutSets array to receive a list of vertices indices for each separated cut
	// @param minViol minimum violation required to add a new cut
	// @return number of separated cuts.
	int separate(int max_cuts, std::vector< std::vector<int> >& cutSets, double minViol);

private:
	Graph instance;

	// find a minimum s-t cut in the graph G and return the cut value
	int findMinCut(int s, int t, std::vector<bool>& cut);

	// calculate the intersection between sets A and B of integers where B is represented as a
	// boolean vector and save the resulting set in A
	// @param complB indicate that B is complemented when true
	void intersect(std::vector<int>& A, std::vector<bool>& B, bool complB);

	// check if a belongs to set A represented as a sorted list
	bool inSet(int a, std::vector<int> A);
};

// function that separates odd degree cutset cuts via construction of separator trees (according to the
// book Network Flows: theory, algorithms, and applications by Ahuja et al)
// @param n number of vertices
// @param m number of (shortest path) arcs with non-zero values
// @param r number of required edges
// @param ps array of start vertices by (shortest path) arc index
// @param pt array of end vertices by (shortest path) arc index
// @param pv array of arc values by (shortest path) arc index
// @param ri array of first vertices by required edge index
// @param ri array of second vertices by required edge index
// @param max_cuts maximum number of cuts to be separated
// @param cut_sets array to receive a sequence of -1 terminated lists of vertices indices corresponding
//                to the separated cuts (the allocated buffer size should be max_cuts * n)
// @param min_viol minimum violation required to add a new cut
// @return number of separated cuts.
int odcs_separate(int n, int m, int r, int* ps, int* pt, double* pv,
		int* ri, int* rj, int max_cuts, int* cut_sets, double min_viol);

