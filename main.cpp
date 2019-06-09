//
// Test routines for the exact separation of (lifted) odd degree cut set inequalities
//
// by Artur Pessoa (2019)
//

#include "OddDegreeCutSep.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

const int max_cuts = 50;

int main(int argc, char* argv[])
{
	// instance data
	int n;
	int m;
	int r;
	std::vector<int> ps;
	std::vector<int> pt;
	std::vector<double> pv;
	std::vector<int> ri;
	std::vector<int> rj;

	// data structure to receive the cuts
	std::vector<int> cut_sets;

	// check the number of program arguments
	if (argc != 2)
	{
		cout << "Use: OddDegreeCutSep <filename>" << endl;
		return -1;
	}

	// try to open the input file
	ifstream f;
	f.open(argv[1], ios_base::in);
	if (!f.is_open())
	{
		cout << "Fail to open instance file " << argv[1] << endl;
		return -2;
	}

	// read the file
	f >> n >> m >> r;
	for (int a = 0; a < m; a++)
	{
		int _ps, _pt;
		double _pv;
		f >> _ps >> _pt >> _pv;
		ps.push_back(_ps);
		pt.push_back(_pt);
		pv.push_back(_pv);
	}
	for (int e = 0; e < r; e++)
	{
		int _ri, _rj;
		f >> _ri >> _rj;
		ri.push_back(_ri);
		rj.push_back(_rj);
	}

	// fclose the input file
	f.close();

	// allocate space for the cuts
	cut_sets.resize(max_cuts * n);

	// call the separation function
	int ncuts = odcs_separate(n, m, r, &ps[0], &pt[0], &pv[0], &ri[0], &rj[0], max_cuts,
			&cut_sets[0], 0.001);

	// print the cuts found
	if (ncuts == 0)
		cout << "No cuts found" << endl;
	else
	{
		cout << "Found " << ncuts << " cuts!" << endl;
		int l = 0;
		for (int c = 0; c < ncuts; c++)
		{
			cout << (c + 1) << ":";
			for (; cut_sets[l] != -1; l++)
				cout << " " << cut_sets[l];
			cout << endl;
			l++;
		}
	}
	return 0;
}

