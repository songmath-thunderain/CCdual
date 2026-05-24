/*- mode: C++;
 * 
 * Date: July 3, 2014
 */

#include <ilcplex/ilocplex.h>
#include <vector>
#include <set>

using namespace std;

typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloNumVarArray> NumVarMatrix;

class IndexVal {
public:
	int ind;
	double val;
	IndexVal() {}
	IndexVal(const int& i, const double& v) 
		{ ind = i; val = v; }
};

bool operator<(const IndexVal& a, const IndexVal& b) {
	if (a.val > b.val)
		return true;
	else if (a.val < b.val)
		return false;
	else if (a.ind > b.ind)
		return true;
	else
		return false;
}

bool operator>(const IndexVal& a, const IndexVal& b) {
	if (a.val < b.val)
		return true;
	else if (a.val > b.val)
		return false;
	else if (a.ind < b.ind)
		return true;
	else
		return false;
}

struct Knapsack
{
    IloInt nbItems;
    IloInt nbScens;
	IloInt nbDims;
    NumMatrix weight;
    IloNumArray capacity;
	IloNumArray cost;
	double totalcost;

};

struct Subprob
{
	// subproblem
	IloCplex cplex;
	IloModel model;
	IloEnv env;
	IloNumVarArray alpha;
	IloNumVarArray lambda;
	IloNumVarArray gamma1;
	IloNumVarArray gamma2;
	IloNumVarArray beta1;
	IloNumVarArray beta2;
	IloNumVarArray xi;
	IloObjective obj;
	IloRangeArray rangelist;
};

struct BendersCoef
{
	// Coefficients for Benders
	vector<double> xcoef;
	int zind;
	double zUcoef;
	double zrest;
	double ycoef;
	double rhsUcoef;
	double rhsrest;
	double normviol;
	double viol;
};

struct Model
{
    IloNumVarArray x;
    IloNumVarArray z;
	IloNumVarArray u;
	IloNumVarArray w;
	IloNumVar y;
	int p;
	int option;
	int solveoption;
	int subproboption;
	int lastnode;
	vector<bool> fixscen;
	Subprob subprob;
	double U;
	IloEnv env;
	double subLPtime;
	double Utime;
	double rootlb;
	int ncuts;
	int nUsercuts;
	int nLazycuts;
	double heurtime;
	double zthreshold;
	double maxpar;
	double cutviol;
	bool dopriority;
	vector< vector<double> > bigM;
	double fivetime;
	double onetime;
	bool fiveflag;
	bool oneflag;
};

struct entiresubprob
{
	IloModel model;
	IloCplex cplex;
	IloEnv env;
	NumVarMatrix x;
};

struct entireprob
{
	IloModel model;
	IloCplex cplex;
	IloEnv env;
	IloNumVarArray x;
	IloNumVarArray pi;
};


double PPlus(IloEnv& env, const Knapsack& knap, const Model& mod, const vector< vector<bool> >& fixscen, int i, int k);
double single(const Knapsack& knap, const Model& mod, int i, int ii, int k, int kk);
double quantile(IloEnv& env, const Knapsack& knap, Model& mod);
void buildSub(const Knapsack& knap, Model& mod, int i);
void genCuts(IloEnv& env, const Knapsack& knap, Model& mod, const IloNumArray& xvals, const IloNumArray& zvals, double yval, double U, vector<BendersCoef>& bcoef, int option, bool bigMoption, double cutviol);
double Heuristic(const Knapsack& knap, const Model& mod, vector<double>& xsols);
bool checkFeas(int i, const vector<double>& xsol, const Knapsack& knap);
double PPlus(const Knapsack& knap, const Model& mod, int i, int k, vector< vector< vector<double> > >& total_pool);
void single_sort(const Knapsack& knap, const Model& mod, int i, vector< vector<double> >& pool);
void lp(IloEnv& env, const Knapsack& knap, const Model& mod, IloTimer& timer, double& lpval, double& lptime);

void rootLP(IloEnv& env, const Knapsack& knap, Model& mod, double U, int option, bool bigMoption, IloTimer& timer, double& val, double& time);
double V2LP(IloEnv& env, const Knapsack& knap, const Model& mod, double UB);
double V2LPiter(IloEnv& env, const Knapsack& knap, const Model& mod, double UB);

double V1LP(IloEnv& env, const Knapsack& knap, const Model& mod);

bool check_parallel(const Model& mod, vector<BendersCoef>& added, BendersCoef bcoef, vector<double>& norm_added);

double par(const BendersCoef& bcoef, const BendersCoef& added, double norm_added, double normA);

double rootLP2(IloEnv& env, const Knapsack& knap, Model& mod, double U, int option, bool bigMoption);
