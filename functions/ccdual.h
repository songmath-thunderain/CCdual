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

struct Model
{
    IloNumVarArray x;
    IloNumVarArray z;
	int p;
	double rootlb;
    double roottime;
    double rootlptime;
	int option;
	int exactNnodes;
	int V1IPNnodes; 
	int V2IPNnodes;
	int V2IPlocalNnodes;
	vector<double> objvals;
	vector<int> fix_z;
	IloEnv env;
	int lastnode;
};

struct Model2
{
    IloNumVarArray x;
    IloIntVarArray z;
	IloNumVarArray u;
	IloNumVarArray w;
	IloNumVar y;
};

 
struct subprob
{
	IloModel model;
	IloCplex cplex;
	IloEnv env;
	IloNumVarArray x;
	IloObjective obj;
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
bool checkFeas(const Knapsack& knap, const Model& mod, const vector<int>& xsols);
bool compare(const vector< vector<int> >& xsols, int k, int kk);
void checkUnique(const Knapsack& knap, const vector< vector<int> >& xsols, vector<bool>& label);
void buildMIP(IloEnv& env, IloModel& model, const Knapsack& knap, const Model& mod);
void buildLP(IloEnv& env, IloModel& model, const Knapsack& knap, const Model& mod);
void ccdual(IloEnv& env, vector<subprob>& submodels, double& LB, double& UB, const Knapsack& knap, const Model& mod, IloTimer& clock, double& worsttime, double& ldtime, int& niter, int option);
double lr1(const Knapsack& knap, const Model& mod, const vector< vector<double> >& lam, vector< vector<double> >& xvals, vector<subprob>& submodels);
double V1LD(IloEnv& env, const Knapsack& knap, const Model& mod, vector<subprob>& submodels);
double quantile(IloEnv& env, const Knapsack& knap, Model& mod);
double V2LPiter(IloEnv& env, const Knapsack& knap, const Model& mod, double UB);
double V2LP(IloEnv& env, const Knapsack& knap, const Model& mod, double U, int& iterV2);
double V3LD(IloEnv& env, const Knapsack& knap, const Model& mod, double U, double L, int& iterV3);
double V3LDiter(IloEnv& env, const Knapsack& knap, const Model& mod, double l);
double V1LP(IloEnv& env, const Knapsack& knap, Model& mod, int option);
double exact(IloEnv& env, const Knapsack& knap, Model& mod, int option);
double quantile_nogood(IloEnv& env, const Knapsack& knap, const Model& mod, const vector< vector<int> >& exploredxsols);
double V2LPiter_nogood(IloEnv& env, const Knapsack& knap, const Model& mod, double UB, const vector< vector<int> >& exploredxsols);
double V2LP_nogood(IloEnv& env, const Knapsack& knap, const Model& mod, double UB, const vector< vector<int> >& exploredxsols);
double Heuristic(const Knapsack& knap, const Model& mod);
double V2IP(IloEnv& env, const Knapsack& knap, Model& mod, double UB, int option);
double quantile_fixed(const Knapsack& knap, const Model& mod, const IloNumArray& ubs, const IloNumArray& lbs);
