/*- mode: C++;
 * Date: July 18, 2014
 */

#include <iostream>
#include <ilcplex/ilocplex.h>
#include "ccdual.h"
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

class DepthData : public IloCplex::MIPCallbackI::NodeData {
   unsigned int const depth; /**< Depth of node at which this data
                                 *   is stored. */
   public:
      DepthData(unsigned int idepth) : depth(idepth) {}
      unsigned int getDepth() const { return depth; }
};
 
DepthData const rootDepth(0);

struct BranchCallback : public IloCplex::BranchCallbackI {
   BranchCallback(IloEnv env) : IloCplex::BranchCallbackI(env) {}
 
   IloCplex::CallbackI *duplicateCallback() const {
         return new (getEnv()) BranchCallback(getEnv());
   }
 
   void main() {
         // How many branches would CPLEX create?
         IloInt const nbranch = getNbranches();
         if (nbranch > 0) {
            // CPLEX would branch. Get the branches CPLEX would create
            // and create exactly those branches. With each branch store
            // its depth in the tree.
            // Note that getNodeData() returns NULL for the root node.
            DepthData const *data = dynamic_cast<DepthData *>(getNodeData());
            if ( !data )
               data = &rootDepth;
 
            IloNumVarArray vars(getEnv());
            IloNumArray bounds(getEnv());
            IloCplex::BranchDirectionArray dirs(getEnv());
 
            for (IloInt i = 0; i < nbranch; ++i) {
               IloNum const est = getBranch(vars, bounds, dirs, i);
               makeBranch(vars, bounds, dirs, est,
                          new DepthData(data->getDepth() + 1U));
            }
            dirs.end();
            bounds.end();
            vars.end();
         }
         else {
            // CPLEX would not create any branch here. Prune the node.
            prune();
         }
   }
};

struct UserCutCallback : public IloCplex::UserCutCallbackI {
	Model2& modlocal;
	const Knapsack& knap;
	Model& mod;
	public:
	UserCutCallback(IloEnv env, Model2& modlocalee, const Knapsack& knapee, Model& modee) : IloCplex::UserCutCallbackI(env), modlocal(modlocalee), knap(knapee), mod(modee) {	
	}
	IloCplex::CallbackI *duplicateCallback() const {
	   return new (getEnv()) UserCutCallback(*this);
	}

	void main() {
		// start cutcallback
		DepthData const *data = dynamic_cast<DepthData *>(getNodeData());
        if ( !data )
           data = &rootDepth;
		int depth = data->getDepth();

		IloEnv env = getEnv();
		int nNodes = getNnodes();
		if (nNodes > 0 && depth <= 4 && nNodes != mod.lastnode)
		{
			double incumbent = 0;
			if (hasIncumbent())
				incumbent = getIncumbentObjValue();
			IloNumArray xvals(env);
			IloNumArray zvals(env);
			IloNumArray uvals(env);
			IloNumArray wvals(env);
			getValues(xvals, modlocal.x);
			getValues(zvals, modlocal.z);
			getValues(uvals, modlocal.u);
			getValues(wvals, modlocal.w);
			double yval = getValue(modlocal.y);
			IloNumArray ubs(env);
			getUBs(ubs, modlocal.z);
			IloNumArray lbs(env);
			getLBs(lbs, modlocal.z);
			mod.fix_z = vector<int>(knap.nbScens, -1);
			for (int i = 0; i < knap.nbScens; ++i)
			{
				if (ubs[i] > 1-1e-5 && lbs[i] > 1-1e-5)
					mod.fix_z[i] = 1;
				if (ubs[i] < 1e-5 && lbs[i] < 1e-5)
					mod.fix_z[i] = 0;
			}
			double upper = quantile_fixed(knap, mod, ubs, lbs);
			ubs.end();
			lbs.end();
			int iterV3 = 0;
			cout << "upper = " << upper << ", incumbent = " << incumbent << endl;
			double newbound = V3LD(mod.env, knap, mod, upper, incumbent, iterV3);
			cout << "newbound = " << newbound << endl;
			mod.fix_z.clear();
			for (int i = 0; i < knap.nbScens; ++i)
			{
				// local cut for u
				double lhsval = 0;
				for (int j = 0; j < knap.nbItems; ++j)
					lhsval += uvals[i*knap.nbItems+j]*knap.cost[j];
				lhsval -= zvals[i]*newbound;
				lhsval -= yval;
				if (lhsval + newbound < -1e-3)
				{
					// add the local cut
					IloExpr lhs(env);
					for (int j = 0; j < knap.nbItems; ++j)
						lhs += modlocal.u[i*knap.nbItems+j]*knap.cost[j];
					lhs -= modlocal.z[i]*newbound;
					lhs -= modlocal.y;
					IloRange range(env, -newbound, lhs, IloInfinity);
					addLocal(range);
					lhs.end();
				}
				// local cut for w
				lhsval = 0;
				for (int j = 0; j < knap.nbItems; ++j)
					lhsval += wvals[i*knap.nbItems+j]*knap.cost[j];
				lhsval += zvals[i]*newbound;
				lhsval -= yval;
				if (lhsval < -1e-3)
				{
					// add the local cut
					IloExpr lhs(env);
					for (int j = 0; j < knap.nbItems; ++j)
						lhs += modlocal.w[i*knap.nbItems+j]*knap.cost[j];
					lhs += modlocal.z[i]*newbound;
					lhs -= modlocal.y;
					IloRange range(env, 0, lhs, IloInfinity);
					addLocal(range);
					lhs.end();
				}
			}
			xvals.end();
			zvals.end();
			uvals.end();
			wvals.end();
		}
		mod.lastnode = getNnodes();
	}
};


int main(int argc, char **argv)
{
    IloEnv env;
    Knapsack knap;
	Model mod;
	mod.option = atoi(argv[3]);
	// option = 1: solve ccdual by quantile bound
	// option = 2: solve ccdual by v2LP using the LP relaxation
	// option = 3: Get the LP bound and LD bound of the multi-dim [0,1] knapsack
    knap.weight = NumMatrix(env);
    knap.cost = IloNumArray(env);
    knap.capacity = IloNumArray(env);
	double eps = atof(argv[2]);
    const char* filename;
    if (argc > 1)
		filename = argv[1];
    ifstream file(filename);
    if (!file)
    {
	cerr << "ERROR: could not open file '" << filename
	     << "' for reading" << endl;
	cerr << "usage:  " << argv[0] << " <file>" << endl;
        throw(-1);
    }
    file >> knap.cost >> knap.capacity >> knap.weight;
    knap.nbItems = knap.cost.getSize();
	knap.nbDims = knap.capacity.getSize();
    knap.nbScens = knap.weight.getSize();
	mod.p = int(eps*knap.nbScens+1e-5);
    cout << "Reliability level:" << eps << endl;
    cout << "Number of items:" << knap.nbItems << endl;
    cout << "Number of scenarios:" << knap.nbScens << endl;	
	cout << "Number of dims:" << knap.nbDims << endl;
    cout << "##################################" << endl;

	knap.totalcost = 0;
	for (int j = 0; j < knap.nbItems; ++j)
		knap.totalcost += knap.cost[j];
    
    IloTimer clock(env);
    clock.start();
    IloNum buildtime, solvetime, lasttime;

	/*
		Build MIP model
	*/   
	if (mod.option == 1 || mod.option == 2)
	{
		mod.x = IloNumVarArray(env, knap.nbItems, 0, 1, ILOINT);
    	mod.z = IloNumVarArray(env, knap.nbScens, 0, 1, ILOINT);
    	
		IloModel model(env);
 		IloCplex cplex(model);
    	// This is equivalent to IloCplex cplex(env) + cplex.extract(model)
		buildMIP(env, model, knap, mod);
    
  	  	cplex.setParam(IloCplex::TiLim,7200);
		cplex.setParam(IloCplex::Threads, 1);

		double startcplex = clock.getTime();
		cplex.solve();
		double cplextime = clock.getTime()-startcplex;
		double cplexobj = cplex.getObjValue();
		vector<subprob> submodels;
		for (int k = 0; k < knap.nbScens; ++k)
		{
			subprob submodel;
			IloEnv subenv;
			submodel.env = subenv;
			submodel.x = IloNumVarArray(submodel.env, knap.nbItems, 0, 1, ILOINT);
			submodel.model = IloModel(submodel.env);
			for (int d = 0; d < knap.nbDims; ++d)
			{
				IloExpr lhs(submodel.env);			
				for (int j = 0; j < knap.nbItems; ++j)
					lhs += submodel.x[j]*knap.weight[k][d*knap.nbItems+j];
				submodel.model.add(lhs <= knap.capacity[d]);
				lhs.end();
			}
			submodel.obj = IloObjective(submodel.env);
			IloExpr obj(submodel.env);
			for (int j = 0; j < knap.nbItems; ++j)
				obj += submodel.x[j]*knap.cost[j];
			submodel.obj = IloMaximize(submodel.env, obj);
			submodel.model.add(submodel.obj);
			obj.end();

			submodel.cplex = IloCplex(submodel.model);
			submodel.cplex.setParam(IloCplex::MIPDisplay, 0);
			submodels.push_back(submodel);
		}
		/*
		   Decomposed Lagrangian
		*/
		cout << "cplexobj = " << cplexobj << endl;
		cout << "cplextime = " << cplextime << endl;
		exit(0);
		double starttime = clock.getTime();
		cout << "start the main loop..." << endl;
		// Main loop
		double LB = 0; // Later on, LB could be replace by Weijun's heuristic
		double UB = 1e8; 
		double worsttime = 0; // record the worst case run time for subproblem to "mimic" a parallel
		double ldtime = 0; // record the time spent on solving the V2LP
		int niter = 0; // record the number of iterations
		if (mod.option == 1)
			ccdual(env, submodels, LB, UB, knap, mod, clock, worsttime, ldtime, niter, 1);
		if (mod.option == 2)
			ccdual(env, submodels, LB, UB, knap, mod, clock, worsttime, ldtime, niter, 2);
		cout << "UB = " << UB << ", LB = " << LB << endl;
		double soltime = clock.getTime()-starttime;
		cout << "soltime = " << soltime << endl;
		cout << "worsttime = " << worsttime << endl;
		cout << "ldtime = " << ldtime << endl;
		cout << "niter = " << niter << endl;

		ofstream out(argv[4], ios::app);
		if (out)  {
				out << setw(24) << argv[1];
				out << ", ";
				out << setw(8) << argv[2];
				out << ", ";
				out << setw(12) << cplexobj;
				out << ", ";
				out << setw(12) << UB;
				out << ", ";
				out << setw(12) << LB;
				out << ", ";
				out << setw(12) << cplextime;
				out << ", ";
				out << setw(12) << soltime;
				out << ", ";
				out << setw(12) << worsttime;
				out << ", ";
				out << setw(12) << ldtime;
				out << ", ";
				out << setw(8) << niter;
				out << endl;
				out.close();
		}	
		for (int k = 0; k < knap.nbScens; ++k)
		{
			submodels[k].cplex.end();
			submodels[k].model.end();
			submodels[k].x.end();
			submodels[k].env.end();
		}
	}

	if (mod.option == 3)
	{
		double exactobj = exact(env, knap, mod,1);
		cout << "exact = " << exactobj << endl;
		double lpobj = exact(env, knap, mod,0);
		cout << "lpobj = " << lpobj << endl;
		double quantileval = quantile(env, knap, mod);
		cout << "quantile = " << quantileval << endl;
		double starttime = clock.getTime();
		// ccdual: standard NAC dual
		//double ld1 = V1LD(env, knap, mod, submodels);
		//double soltime = clock.getTime()-starttime;
		//cout << "V1 dual = " << ld1 << endl;
		//cout << "soltime = " << soltime << endl;
		double lp1 = V1LP(env, knap, mod, 0);
		cout << "V1 primal = " << lp1 << endl;
		int iterV2 = 0;
		double lp2 = V2LP(env,knap, mod, quantileval, iterV2);
		cout << "V2 primal = " << lp2 << ", iterV2 = " << iterV2 << endl;
		double lb = 0;
		int iterV3 = 0;
		double ld3 = V3LD(env,knap, mod, quantileval, lb, iterV3);
		cout << "V3 dual = " << ld3 << ", iterV3 = " << iterV3 << endl;
		double heuristic = Heuristic(knap, mod);
		cout << "heuristic = " << heuristic << endl;
		ofstream out(argv[4], ios::app);
		if (out)  {
				out << setw(24) << argv[1];
				out << ", ";
				out << setw(8) << argv[2];
				out << ", ";
				//out << setw(12) << exactobj;
				//out << ", ";
				out << setw(12) << lpobj;
				out << ", ";
				out << setw(12) << quantileval;
				out << ", ";
				out << setw(12) << lp1;
				out << ", ";
				out << setw(12) << lp2;
				out << ", ";
				out << setw(8) << iterV2;
				out << ", ";
				out << setw(12) << ld3;
				out << ", ";
				out << setw(8) << iterV3;
				out << ", ";
				out << setw(12) << heuristic;
				out << endl;
				out.close();
		}	
	}
	if (mod.option == 4)
	{
		double starttime = clock.getTime();
		double exactobj = exact(env, knap, mod, 1);
		double exacttime = clock.getTime()-starttime;
		cout << "exactobj = " << exactobj << ", exactNnodes = " << mod.exactNnodes << endl;
		cout << "exacttime = " << exacttime << endl;

		double quantileval = quantile(env, knap, mod);

		starttime = clock.getTime();
		double V1IPobj = V1LP(env, knap, mod, 1);
		cout << "V1IPobj = " << V1IPobj << ", V1IPNodes = " << mod.V1IPNnodes << endl;
		double V1IPtime = clock.getTime()-starttime;
		cout << "V1IPtime = " << V1IPtime << endl;
		starttime = clock.getTime();
		int iterV3 = 0;
		double rootlp = V3LD(env,knap, mod, quantileval, 0, iterV3);
		// UB comes from V2LP at root, option 0: no local cuts, just use the primal IP formulation with a root quantile-based LD bound 
		double V2IPobj = V2IP(env, knap, mod, rootlp, 0);
		cout << "V2IPobj = " << V2IPobj << ", V2IPNodes = " << mod.V2IPNnodes << endl;
		double V2IPtime = clock.getTime()-starttime;
		cout << "V2IPtime = " << V2IPtime << endl;
		IloEnv env2;
		mod.env = env2;
		starttime = clock.getTime();
		double V2IPlocalobj = V2IP(env, knap, mod, rootlp, 1);
		double V2IPlocaltime = clock.getTime()-starttime;
		cout << "V2IPlocalobj = " << V2IPlocalobj << ", V2IPlocalNodes = " << mod.V2IPlocalNnodes << endl;
		cout << "V2IPlocaltime = " << V2IPlocaltime << endl;

		ofstream out(argv[4], ios::app);
		if (out)  {
				out << setw(24) << argv[1];
				out << ", ";
				out << setw(8) << argv[2];
				out << ", ";
				out << setw(12) << exactobj;
				out << ", ";
				out << setw(12) << V1IPobj;
				out << ", ";
				out << setw(12) << V2IPobj;
				out << ", ";
				out << setw(12) << exacttime;
				out << ", ";
				out << setw(12) << mod.exactNnodes;
				out << ", ";
				out << setw(12) << V1IPtime;
				out << ", ";
				out << setw(12) << mod.V1IPNnodes;
				out << ", ";
				out << setw(12) << V2IPtime;
				out << ", ";
				out << setw(12) << mod.V2IPNnodes;
				out << ", ";
				out << setw(12) << V2IPlocaltime;
				out << ", ";
				out << setw(12) << mod.V2IPlocalNnodes;
				out << endl;
				out.close();
		}	
	}
 	env.end();
    return 0;
}   

double single(const Knapsack& knap, const Model& mod, int i, int ii, int k, int kk)
{
	double threshold = 1e-5;
	vector<IndexVal> cost;
	double returnval = 0.0;
	for (int j = 0; j < knap.nbItems; ++j)
	{
		if (knap.weight[ii][kk*knap.nbItems+j] < threshold)
			returnval += knap.weight[i][k*knap.nbItems+j];
		else
		{
			IndexVal indexval;
			indexval.ind = j;
			indexval.val = knap.weight[i][k*knap.nbItems+j]*1.0/knap.weight[ii][kk*knap.nbItems+j];
			cost.push_back(indexval);
		}
	}
	sort(cost.begin(), cost.end(), less<IndexVal>());
	double currentweight = 0.0;
	int iter = 0;
	while (currentweight < knap.capacity[kk] && iter < cost.size())
	{
		if (currentweight + knap.weight[ii][kk*knap.nbItems+cost[iter].ind] <= knap.capacity[kk])
		{
			currentweight += knap.weight[ii][kk*knap.nbItems+cost[iter].ind];
			returnval += knap.weight[i][k*knap.nbItems+cost[iter].ind];
		}
		else
		{
			returnval += knap.weight[i][k*knap.nbItems+cost[iter].ind]*(knap.capacity[kk]-currentweight)*1.0/knap.weight[ii][kk*knap.nbItems+cost[iter].ind];
			break;
		}
		iter++;
	}
	returnval -= knap.capacity[k];
	return returnval;
}

double PPlus(IloEnv& env, const Knapsack& knap, const Model& mod, const vector< vector<bool> >& fixscen, int i, int k)
{
	double returnval = 1e8;
	for (int kk = 0; kk < knap.nbDims; ++kk)
	{
		vector<IndexVal> eta;
		for (int ii = 0; ii < knap.nbScens; ++ii)
		{
			if (ii == i && kk == k)
			{
				IndexVal indexval;
				indexval.ind = i;
				double temp = 0.0;
				for (int j = 0; j < knap.nbItems; ++j)
					temp += knap.weight[i][k*knap.nbItems+j];
				if (temp-knap.capacity[k] > 0)
					indexval.val = 0;
				else
					indexval.val = temp-knap.capacity[k];
				eta.push_back(indexval);
			}
			else
			{
				if (fixscen[kk][ii] == 0)
				{
					IndexVal indexval;
					indexval.ind = ii;
					indexval.val = single(knap, mod, i, ii, k, kk);
					eta.push_back(indexval);
				}
			}
		}		
		sort(eta.begin(), eta.end(), greater<IndexVal>());
		double subval = eta[mod.p].val;
		if (subval < returnval)
			returnval = subval;
	}
	return returnval;
}

void checkUnique(const Knapsack& knap, const vector< vector<int> >& xsols, vector<bool>& label)
{
	// Check the uniqueness of the subproblem solutions, label them as "1"
	label[0] = 1;
	for (int k = 1; k < knap.nbScens; ++k)
	{
		bool flag = 1;
		for (int kk = 0; kk < k; ++kk)
		{
			if (compare(xsols, k, kk) == 1)
			{
				flag = 0;
				break;
			}
		}
		if (flag == 1)
			label[k] = 1;
	}
}


bool compare(const vector< vector<int> >& xsols, int k, int kk)
{
	bool flag = 1;
	for (int j = 0; j < xsols[k].size(); ++j)
	{
		if (xsols[k][j] != xsols[kk][j])
		{
			flag = 0;
			break;
		}
	}
	return flag;
}

bool checkFeas(const Knapsack& knap, const Model& mod, const vector<int>& xsols)
{
	int succ = 0;
	for (int k = 0; k < knap.nbScens; ++k)
	{
		bool flag = 1;
		for (int d = 0; d < knap.nbDims; ++d)
		{
			double lhsval = 0;
			for (int j = 0; j < knap.nbItems; ++j)
			{
				if (xsols[j] == 1)
					lhsval += knap.weight[k][d*knap.nbItems+j];
			}
			if (lhsval > knap.capacity[d]+1e-5)
			{
				flag = 0;
				break;
			}
		}
		if (flag == 1)
			succ++;
	}
	if (succ >= knap.nbScens-mod.p)
		return 1;
	else
		return 0;
}


void buildMIP(IloEnv& env, IloModel& model, const Knapsack& knap, const Model& mod)
{
	model.add(IloSum(mod.z) >= (knap.nbScens-mod.p));
    // reliability constraint

	vector< vector<double> > bigM;
	
	vector<bool> fixscen(knap.nbScens, 0);

	for (int i = 0; i < knap.nbScens; ++i)
	{
		vector<double> tempBigM(knap.nbDims);
		bigM.push_back(tempBigM);
		bool fflag = 1;
    	for (int k = 0; k < knap.nbDims; k++)
    	{
			double totalweight = 0.0;
        	for (int j = 0; j < knap.nbItems; j++)
				totalweight += knap.weight[i][k*knap.nbItems+j];
			bigM[i][k] = totalweight - knap.capacity[k];
			if (bigM[i][k] > 1e-5)
				fflag = 0;
    	}
		if (fflag == 1)
			fixscen[i] = 1;
	}

	// Add big-M valid inequalities for (x,z) space
	for (int i = 0; i < knap.nbScens; ++i)
	{
		if (fixscen[i] == 1)
			model.add(mod.z[i] >= 1);
		else
		{
			for (int k = 0; k < knap.nbDims; ++k)
			{
				IloExpr lhs(env);
				for (int j = 0; j < knap.nbItems; ++j)
					lhs += mod.x[j]*knap.weight[i][k*knap.nbItems+j];
				lhs += mod.z[i]*bigM[i][k];
				model.add(lhs <= (knap.capacity[k]+bigM[i][k]));
				lhs.end();
			}
		}
	}

    // big-M constraint
    
    IloExpr obj(env);

    for(int j = 0; j < knap.nbItems; j++)
        obj += mod.x[j]*knap.cost[j];
    model.add(IloMaximize(env, obj));
    obj.end();
}

void buildLP(IloEnv& env, IloModel& model, const Knapsack& knap, const Model& mod)
{
	model.add(IloSum(mod.z) >= (knap.nbScens-mod.p));
    // reliability constraint

	vector< vector<double> > bigM;
	
	vector<bool> fixscen(knap.nbScens, 0);

	for (int i = 0; i < knap.nbScens; ++i)
	{
		vector<double> tempBigM(knap.nbDims);
		bigM.push_back(tempBigM);
		bool fflag = 1;
    	for (int k = 0; k < knap.nbDims; k++)
    	{
			double totalweight = 0.0;
        	for (int j = 0; j < knap.nbItems; j++)
				totalweight += knap.weight[i][k*knap.nbItems+j];
			bigM[i][k] = totalweight - knap.capacity[k];
			if (bigM[i][k] > 1e-5)
				fflag = 0;
    	}
		if (fflag == 1)
			fixscen[i] = 1;
	}

	// Add big-M valid inequalities for (x,z) space
	for (int i = 0; i < knap.nbScens; ++i)
	{
		if (fixscen[i] == 1)
			model.add(mod.z[i] >= 1);
		else
		{
			for (int k = 0; k < knap.nbDims; ++k)
			{
				IloExpr lhs(env);
				for (int j = 0; j < knap.nbItems; ++j)
					lhs += mod.x[j]*knap.weight[i][k*knap.nbItems+j];
				lhs += mod.z[i]*bigM[i][k];
				model.add(lhs <= (knap.capacity[k]+bigM[i][k]));
				lhs.end();
			}
		}
	}

    // big-M constraint
    
    IloExpr obj(env);

    for(int j = 0; j < knap.nbItems; j++)
        obj += mod.x[j]*knap.cost[j];
    model.add(IloMaximize(env, obj));
    obj.end();
}


void ccdual(IloEnv& env, vector<subprob>& submodels, double& LB, double& UB, const Knapsack& knap, const Model& mod, IloTimer& clock, double& worsttime, double& ldtime, int& niter, int option)
{
	// Just use the quantile bound as the UB
	double starttime = clock.getTime();
	vector< vector<int> > exploredxsols;
	while ((UB-LB)*1.0/UB > 1e-3)
	{
		niter++;
		cout << "UB = " << UB << ", LB = " << LB << endl;
		// Solve each scenario subproblem
		vector< vector<int> > xsols;
		vector<double> objs;
		double subtime = 0;		
		for (int k = 0; k < knap.nbScens; ++k)
		{
			double start = clock.getTime();
			submodels[k].cplex.solve();
			double temp = clock.getTime()-start;
			if (temp > subtime)
				subtime = temp;
			if (submodels[k].cplex.getStatus() == IloAlgorithm::Feasible || submodels[k].cplex.getStatus() == IloAlgorithm::Optimal)
			{
				double obj = submodels[k].cplex.getObjValue();
				objs.push_back(obj);
				IloNumArray xvals(submodels[k].env);
				submodels[k].cplex.getValues(submodels[k].x, xvals);
				vector<int> xsol(knap.nbItems, 0);
				for (int j = 0; j < knap.nbItems; ++j)
				{
					if (xvals[j] > 1e-5)
						xsol[j] = 1;
				}
				xvals.end();
				xsols.push_back(xsol);
			}
			else
			{
				cout << "Subproblem " << k << " is infeasible!" << endl;
				exit(0);
			}
		}
		worsttime += subtime;
		vector<IndexVal> subcost;
		for (int k = 0; k < knap.nbScens; ++k)
		{
			IndexVal indexval;
			indexval.ind = k;
			indexval.val = objs[k];
			subcost.push_back(indexval);
		}
		sort(subcost.begin(), subcost.end(), less<IndexVal>());

		if (option == 1)
		{
			// use the quantile bound as the UB
			double currentUB2 = subcost[knap.nbScens-mod.p-1].val;
			if (currentUB2 < UB)
				UB = currentUB2;
		}
	
		if (option == 2)
		{
			// use the iterative LP solution as the bound
			// first calculte quantile
			double ldstarttime = clock.getTime();
			double upper = quantile_nogood(env, knap, mod, exploredxsols);
			double currentUB = V2LP_nogood(env, knap, mod, upper, exploredxsols);
			ldtime += (clock.getTime()-ldstarttime);
			if (currentUB < UB)
				UB = currentUB;
		
			// TEST BEGIN - Just for a comparison
			// use the quantile bound as the UB
			double currentUB2 = subcost[knap.nbScens-mod.p-1].val;
			cout << "quantile = " << currentUB2 << ", V2LP = " << currentUB << endl;
			// TEST END
		}
		vector<bool> label(knap.nbScens, 0);
		checkUnique(knap, xsols,label);
		// Update submodels by excluding these unique solutions
		int cutcount = 0;
		for (int kk = 0; kk < knap.nbScens; ++kk)
		{
			if (label[kk] == 1)
			{
				exploredxsols.push_back(xsols[kk]);
				// Add no-good cuts to each scenario
				cutcount++;
				for (int k = 0; k < knap.nbScens; ++k)
				{
					// unique solution
					IloExpr lhs(submodels[k].env);
					int cardpos = 0;
					for (int j = 0; j < knap.nbItems; ++j)
					{
						if (xsols[kk][j] == 0)
							lhs += submodels[k].x[j];
						else
						{
							lhs -= submodels[k].x[j];
							cardpos++;
						}
					}
					submodels[k].model.add(lhs >= 1-cardpos);
					lhs.end();
				}
			}
		}
		// Check the feasibility of these unique solutions, and update the LB accordingly
		// Check from the sorted list, start from the largest, if feasible, then forget about the rest
		for (int k = 0; k < knap.nbScens; ++k)
		{
			int ind = subcost[k].ind;
			if (subcost[k].val <= LB)
				break;
			if (label[ind] == 1)
			{
				if (checkFeas(knap, mod, xsols[ind]) == 1)
				{
					LB = subcost[k].val;
					break;
				}
			}
		}
		if (clock.getTime()-starttime > 7200)
			break;
	}
}

double lr1(const Knapsack& knap, const Model& mod, const vector< vector<double> >& lam, vector< vector<double> >& xvals, vector<subprob>& submodels)
{
	// LR served for V1
	double returnval = 0;
	vector<double> u(knap.nbScens, 0);
	vector<double> v(knap.nbScens, 0);
	for (int i = 0; i < knap.nbScens; ++i)
	{
		for (int j = 0; j < knap.nbItems; ++j)
		{
			double coef = knap.cost[j]*1.0/knap.nbScens+lam[i][j];
			if (coef > 0)
				u[i] += coef;
			submodels[i].obj.setLinearCoef(submodels[i].x[j],coef); 
		}
		submodels[i].cplex.solve();
		if (submodels[i].cplex.getStatus() == IloAlgorithm::Feasible || submodels[i].cplex.getStatus() == IloAlgorithm::Optimal)
		{
			IloNumArray subxvals(submodels[i].env);
			submodels[i].cplex.getValues(subxvals, submodels[i].x);
			for (int j = 0; j < knap.nbItems; ++j)
				xvals[i][j] = subxvals[j];
			v[i] = submodels[i].cplex.getObjValue();
			returnval += u[i];
			subxvals.end();
		}
		else
		{
			cout << "Scenario problem is infeasible!";
			exit(0);
		}

	}
	vector<IndexVal> sortvu;
	for (int i = 0; i < knap.nbScens; ++i)
	{
		IndexVal indexval;
		indexval.ind = i;
		indexval.val = v[i]-u[i];
		sortvu.push_back(indexval);
	}
	sort(sortvu.begin(), sortvu.end(), less<IndexVal>());
	for (int i = 0; i < knap.nbScens; ++i)
	{
		if (i < knap.nbScens-mod.p)
			returnval += sortvu[i].val;
		else
		{
			int ind = sortvu[i].ind;
			for (int j = 0; j < knap.nbItems; ++j)
			{
				if (knap.cost[j]*1.0/knap.nbScens+lam[ind][j] > 0)
					xvals[ind][j] = 1;
				else
					xvals[ind][j] = 0;
			}
		}
	}
	return returnval;
}

double V1LD(IloEnv& env, const Knapsack& knap, const Model& mod, vector<subprob>& submodels)
{
	// LD for relaxing nonanticipativity constraints
	vector< vector<double> > lam(knap.nbScens);
	for (int i = 0; i < knap.nbScens; ++i)
		lam[i] = vector<double>(knap.nbItems, 0);
	double theta = -1e5;
	bool flag = 1;
	IloNumVarArray lambda(env, 1+knap.nbScens*knap.nbItems, -50, 50);
	lambda[knap.nbScens*knap.nbItems].setUB(IloInfinity);
	lambda[knap.nbScens*knap.nbItems].setLB(-IloInfinity);
	IloModel ldmodel(env);
	IloExpr ldobj(env);
	ldobj += lambda[knap.nbScens*knap.nbItems];
	ldmodel.add(IloMinimize(env, ldobj));
	ldobj.end();
	for (int j = 0; j < knap.nbItems; ++j)
	{
		IloExpr lhs(env);
		for (int i = 0; i < knap.nbScens; ++i)
			lhs += lambda[i*knap.nbItems+j];
		ldmodel.add(lhs == 0);
		lhs.end();
	}
	IloCplex ldcplex(ldmodel);
	ldcplex.setParam(IloCplex::SimDisplay, 0);
	vector< vector<double> > xvals(knap.nbScens);
	for (int i = 0; i < knap.nbScens; ++i)
		xvals[i] = vector<double>(knap.nbItems);
	// Start subgradient method
	while (flag == 1)
	{
		double result = lr1(knap, mod, lam, xvals, submodels);
		if (theta < result*(1-1e-5))
		{
			// add a cut
			IloExpr lhs(env);
			lhs += lambda[knap.nbItems*knap.nbScens];
			for (int i = 0; i < knap.nbScens; ++i)
			{
				for (int j = 0; j < knap.nbItems; ++j)
					lhs -= xvals[i][j]*lambda[i*knap.nbItems+j];
			}
			double subtract = 0.0;
			for (int i = 0; i < knap.nbScens; ++i)
			{
				for (int j = 0; j < knap.nbItems; ++j)
					subtract += xvals[i][j]*lam[i][j];
			}
			ldmodel.add(lhs >= result-subtract);
			lhs.end();
			// update lam and theta
			ldcplex.solve();
			if (ldcplex.getStatus() == IloAlgorithm::Optimal)
			{
				IloNumArray temp(env);
				ldcplex.getValues(temp, lambda);
				theta = temp[knap.nbItems*knap.nbScens];
				for (int i = 0; i < knap.nbScens; ++i)
				{
					for (int j = 0; j < knap.nbItems; ++j)
						lam[i][j] = temp[i*knap.nbItems+j];
				}
				temp.end();
			}
			else
			{
				cout << "error!" << endl;
				exit(0);
			}
		}
		else 
			flag = 0;
	}	
	ldcplex.end();
	ldmodel.end();
	lambda.end();
	return theta;
}

double V1LP(IloEnv& env, const Knapsack& knap, Model& mod, int option)
{
	double returnval;
	IloNumVarArray u(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
	IloNumVarArray w(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
	IloNumVarArray x(env, knap.nbItems, 0, 1);
	IloNumVarArray z;
	if (option == 0)
		z = IloNumVarArray(env, knap.nbScens, 0, 1);
	if (option == 1)
		z = IloNumVarArray(env, knap.nbScens, 0, 1, ILOINT);
	IloModel model(env);
	// Nonanticipativity
	for (int j = 0; j < knap.nbItems; ++j)
	{
		for (int i = 0; i < knap.nbScens; ++i)
			model.add(u[i*knap.nbItems+j]+w[i*knap.nbItems+j]-x[j] == 0);
	}
	// Au, Du, Dw
	for (int i = 0; i < knap.nbScens; ++i)
	{
		for (int k = 0; k < knap.nbDims; ++k)
		{
			IloExpr lhs(env);
			for (int j = 0; j < knap.nbItems; ++j)
				lhs += u[i*knap.nbItems+j]*knap.weight[i][k*knap.nbItems+j];
			lhs -= knap.capacity[k]*z[i];
			model.add(lhs <= 0);
		}
		for (int j = 0; j < knap.nbItems; ++j)
		{
			model.add(u[i*knap.nbItems+j]-z[i] <= 0);
			model.add(w[i*knap.nbItems+j]+z[i] <= 1);
		}
	}
	// Knapsack
	model.add(IloSum(z) >= (knap.nbScens-mod.p));

	IloExpr obj(env);
	for (int j = 0; j < knap.nbItems; ++j)
		obj += x[j]*knap.cost[j];
	model.add(IloMaximize(env, obj));
	obj.end();
	IloCplex cplex(model);
	if (option == 0)
		cplex.setParam(IloCplex::SimDisplay, 0);
	if (option == 1)
		cplex.setParam(IloCplex::MIPDisplay, 0);
	cplex.solve();
	if (cplex.getStatus() == IloAlgorithm::Feasible || cplex.getStatus() == IloAlgorithm::Optimal)
	{
		returnval = cplex.getObjValue();
		mod.V1IPNnodes = cplex.getNnodes();
	}
	else
	{
		cout << "Infeasible!" << endl;
		exit(0);
	}
	cplex.end();
	model.end();
	u.end();
	w.end();
	x.end();
	z.end();
	return returnval;
}


double quantile(IloEnv& env, const Knapsack& knap, Model& mod)
{
	vector<IndexVal> sortobjvals;
	for (int i = 0; i < knap.nbScens; ++i)
	{
		IloNumVarArray subx(env, knap.nbItems, 0, 1);
		IloModel submodel(env);
		for (int d = 0; d < knap.nbDims; ++d)
		{
			IloExpr lhs(env);			
			for (int j = 0; j < knap.nbItems; ++j)
				lhs += subx[j]*knap.weight[i][d*knap.nbItems+j];
			submodel.add(lhs <= knap.capacity[d]);
			lhs.end();
		}
		IloObjective subobj(env);
		IloExpr obj(env);
		for (int j = 0; j < knap.nbItems; ++j)
			obj += subx[j]*knap.cost[j];
		subobj = IloMaximize(env, obj);
		submodel.add(subobj);
		obj.end();

		IloCplex subcplex(submodel);
		subcplex.setParam(IloCplex::SimDisplay, 0);
		subcplex.solve();
		double subobjval = subcplex.getObjValue();
		IndexVal indexval;
		indexval.ind = i;
		indexval.val = subobjval;
		sortobjvals.push_back(indexval);
		mod.objvals.push_back(subobjval);
		subcplex.end();
		submodel.end();
		subx.end();
	}
	sort(sortobjvals.begin(), sortobjvals.end(), less<IndexVal>());
	return sortobjvals[knap.nbScens-mod.p-1].val;
}

double V2LPiter(IloEnv& env, const Knapsack& knap, const Model& mod, double UB)
{
	// UB could be initialized as the quantile bound
	double returnval;
	double U = UB;
	IloNumVar y(env, -IloInfinity, IloInfinity);
	IloNumVarArray u(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
	IloNumVarArray w(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
	IloNumVarArray x(env, knap.nbItems, 0, 1);
	IloNumVarArray z(env, knap.nbScens, 0, 1);
	IloModel model(env);
	// Nonanticipativity
	for (int j = 0; j < knap.nbItems; ++j)
	{
		for (int i = 0; i < knap.nbScens; ++i)
			model.add(u[i*knap.nbItems+j]+w[i*knap.nbItems+j]-x[j] == 0);
	}
	// objbound
	for (int i = 0; i < knap.nbScens; ++i)
	{
		IloExpr lhs(env);
		for (int j = 0; j < knap.nbItems; ++j)
			lhs += u[i*knap.nbItems+j]*knap.cost[j];
		lhs -= z[i]*U;
		lhs -= y;
		model.add(lhs >= -U);
		lhs.end();
		IloExpr lhs2(env);
		for (int j = 0; j < knap.nbItems; ++j)
			lhs2 += w[i*knap.nbItems+j]*knap.cost[j];
		lhs2 += z[i]*U;
		lhs2 -= y;
		model.add(lhs2 >= 0);
		lhs2.end();
	}
	// Au, Du, Dw
	for (int i = 0; i < knap.nbScens; ++i)
	{
		for (int k = 0; k < knap.nbDims; ++k)
		{
			IloExpr lhs(env);
			for (int j = 0; j < knap.nbItems; ++j)
				lhs += u[i*knap.nbItems+j]*knap.weight[i][k*knap.nbItems+j];
			lhs -= knap.capacity[k]*z[i];
			model.add(lhs <= 0);
		}
		for (int j = 0; j < knap.nbItems; ++j)
		{
			model.add(u[i*knap.nbItems+j]-z[i] <= 0);
			model.add(w[i*knap.nbItems+j]+z[i] <= 1);
		}
	}
	// Knapsack
	model.add(IloSum(z) >= (knap.nbScens-mod.p));

	IloExpr obj(env);
	obj += y;
	model.add(IloMaximize(env, obj));
	obj.end();
	IloCplex cplex(model);
	cplex.setParam(IloCplex::SimDisplay, 0);

	cplex.solve();
	if (cplex.getStatus() == IloAlgorithm::Feasible || cplex.getStatus() == IloAlgorithm::Optimal)
		returnval = cplex.getObjValue();
	else
	{
		cout << "Infeasible!" << endl;
		exit(0);
	}
	cplex.end();
	model.end();
	u.end();
	w.end();
	x.end();
	z.end();
	return returnval;
}

double V2LP(IloEnv& env, const Knapsack& knap, const Model& mod, double UB, int& iterV2)
{
	double U = UB;
	bool flag = 1;
	while (flag == 1)
	{
		iterV2++;
		double objval = V2LPiter(env, knap, mod, U);
		if (fabs((objval-U)*1.0/U) < 1e-5)
			flag = 0;
		else
			U = objval;
	}
	return U;
}

double V3LDiter(IloEnv& env, const Knapsack& knap, const Model& mod, double l)
{
	// UB could be initialized as the quantile bound
	double returnval;
	IloNumVarArray s(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
	IloNumVarArray t(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
	IloNumVarArray u(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
	IloNumVarArray w(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
	IloNumVarArray x(env, knap.nbItems, 0, 1);
	IloNumVarArray z(env, knap.nbScens, 0, 1);
	IloModel model(env);
	// Nonanticipativity
	for (int j = 0; j < knap.nbItems; ++j)
	{
		for (int i = 0; i < knap.nbScens; ++i)
			model.add(u[i*knap.nbItems+j]+w[i*knap.nbItems+j]+s[i*knap.nbItems+j]-t[i*knap.nbItems+j]-x[j] == 0);
	}
	// objbound
	for (int i = 0; i < knap.nbScens; ++i)
	{
		IloExpr lhs(env);
		for (int j = 0; j < knap.nbItems; ++j)
			lhs += u[i*knap.nbItems+j]*knap.cost[j];
		lhs -= z[i]*l;
		model.add(lhs >= 0);
		lhs.end();
		IloExpr lhs2(env);
		for (int j = 0; j < knap.nbItems; ++j)
			lhs2 += w[i*knap.nbItems+j]*knap.cost[j];
		lhs2 += z[i]*l;
		model.add(lhs2 >= l);
		lhs2.end();
	}
	// Au, Du, Dw
	for (int i = 0; i < knap.nbScens; ++i)
	{
		for (int k = 0; k < knap.nbDims; ++k)
		{
			IloExpr lhs(env);
			for (int j = 0; j < knap.nbItems; ++j)
				lhs += u[i*knap.nbItems+j]*knap.weight[i][k*knap.nbItems+j];
			lhs -= knap.capacity[k]*z[i];
			model.add(lhs <= 0);
		}
		for (int j = 0; j < knap.nbItems; ++j)
		{
			model.add(u[i*knap.nbItems+j]-z[i] <= 0);
			model.add(w[i*knap.nbItems+j]+z[i] <= 1);
		}
	}
	// Knapsack
	model.add(IloSum(z) >= (knap.nbScens-mod.p));

	// Check if there are fixed z: used for local cuts
	if (mod.fix_z.size() > 0)
	{
		for (int i = 0; i < knap.nbScens; ++i)
		{
			if (mod.fix_z[i] == 0)
			{
				for (int j = 0; j < knap.nbItems; ++j)
					u[i*knap.nbItems+j].setUB(0);
				z[i].setUB(0);
			}
			if (mod.fix_z[i] == 1)
			{
				for (int j = 0; j < knap.nbItems; ++j)
					w[i*knap.nbItems+j].setUB(0);
				z[i].setLB(1);
			}
		}
	}

	IloExpr obj(env);
	for (int i = 0; i < knap.nbScens; ++i)
	{
		for (int j = 0; j < knap.nbItems; ++j)
			obj += s[i*knap.nbItems+j]+t[i*knap.nbItems+j];
	}
	model.add(IloMinimize(env, obj));
	obj.end();
	IloCplex cplex(model);
	cplex.setParam(IloCplex::SimDisplay, 0);
	cplex.solve();
	if (cplex.getStatus() == IloAlgorithm::Feasible || cplex.getStatus() == IloAlgorithm::Optimal)
		returnval = cplex.getObjValue();
	else
	{
		cout << "Infeasible!" << endl;
		exit(0);
	}
	cplex.end();
	model.end();
	u.end();
	w.end();
	x.end();
	z.end();
	s.end();
	t.end();
	return returnval;
}

double V3LD(IloEnv& env, const Knapsack& knap, const Model& mod, double U, double L, int& iterV3)
{
	// L: a lb provided by a feasible solution
	// U: quantile bound
	double LB = L;
	double UB = U;
	while (fabs((UB-LB)*1.0/UB) > 1e-5)
	{
		iterV3++;
		double l = (LB+UB)*1.0/2;
		double result = V3LDiter(env, knap, mod, l);
		if (result < 1e-5)
			LB = l;
		else
			UB = l;
	}
	return UB;
}

double exact(IloEnv& env, const Knapsack& knap, Model& mod, int option)
{
	IloNumVarArray x(env, knap.nbItems, 0, 1);
	IloNumVarArray z;
	if (option == 0)
		z = IloNumVarArray(env, knap.nbScens, 0, 1);
	if (option == 1)
		z = IloNumVarArray(env, knap.nbScens, 0, 1, ILOINT);
	IloModel model(env);
	IloCplex cplex(model);
	model.add(IloSum(z) >= (knap.nbScens-mod.p));
    // reliability constraint

	vector< vector<double> > bigM;
	
	vector<bool> fixscen(knap.nbScens, 0);

	for (int i = 0; i < knap.nbScens; ++i)
	{
		vector<double> tempBigM(knap.nbDims);
		bigM.push_back(tempBigM);
		bool fflag = 1;
    	for (int k = 0; k < knap.nbDims; k++)
    	{
			double totalweight = 0.0;
        	for (int j = 0; j < knap.nbItems; j++)
				totalweight += knap.weight[i][k*knap.nbItems+j];
			bigM[i][k] = totalweight - knap.capacity[k];
			if (bigM[i][k] > 1e-5)
				fflag = 0;
    	}
		if (fflag == 1)
			fixscen[i] = 1;
	}

	// Add big-M valid inequalities for (x,z) space
	for (int i = 0; i < knap.nbScens; ++i)
	{
		if (fixscen[i] == 1)
			model.add(z[i] >= 1);
		else
		{
			for (int k = 0; k < knap.nbDims; ++k)
			{
				IloExpr lhs(env);
				for (int j = 0; j < knap.nbItems; ++j)
					lhs += x[j]*knap.weight[i][k*knap.nbItems+j];
				lhs += z[i]*bigM[i][k];
				model.add(lhs <= (knap.capacity[k]+bigM[i][k]));
				lhs.end();
			}
		}
	}

    // big-M constraint
    
    IloExpr obj(env);

    for(int j = 0; j < knap.nbItems; j++)
        obj += x[j]*knap.cost[j];
    model.add(IloMaximize(env, obj));
    obj.end();

	if (option == 0)
		cplex.setParam(IloCplex::SimDisplay, 0);
	if (option == 1)
		cplex.setParam(IloCplex::MIPDisplay, 0);
	cplex.solve();
	double returnval = cplex.getObjValue();
	if (option == 1)
		mod.exactNnodes = cplex.getNnodes();
	cplex.end();
	model.end();
	z.end();
	x.end();
	return returnval;
}

double quantile_nogood(IloEnv& env, const Knapsack& knap, const Model& mod, const vector< vector<int> >& exploredxsols)
{
	// calculate the quantile bound including no good cuts
	vector<IndexVal> sortobjvals;
	for (int i = 0; i < knap.nbScens; ++i)
	{
		IloNumVarArray subx(env, knap.nbItems, 0, 1);
		IloModel submodel(env);
		for (int d = 0; d < knap.nbDims; ++d)
		{
			IloExpr lhs(env);			
			for (int j = 0; j < knap.nbItems; ++j)
				lhs += subx[j]*knap.weight[i][d*knap.nbItems+j];
			submodel.add(lhs <= knap.capacity[d]);
			lhs.end();
		}
		// Add no-good cuts for each scenario
		for (int ii = 0; ii < exploredxsols.size(); ++ii)
		{
			IloExpr lhs(env);
			int cardpos = 0;
			for (int jj = 0; jj < knap.nbItems; ++jj)
			{	
				if (exploredxsols[ii][jj] == 0)
					lhs += subx[jj];
				else
				{
					lhs -= subx[jj];
					cardpos++;
				}
			}
			submodel.add(lhs >= 1-cardpos);
			lhs.end();
		}	
				
		IloObjective subobj(env);
		IloExpr obj(env);
		for (int j = 0; j < knap.nbItems; ++j)
			obj += subx[j]*knap.cost[j];
		subobj = IloMaximize(env, obj);
		submodel.add(subobj);
		obj.end();

		IloCplex subcplex(submodel);
		subcplex.setParam(IloCplex::SimDisplay, 0);
		subcplex.solve();
		double subobjval = subcplex.getObjValue();
		IndexVal indexval;
		indexval.ind = i;
		indexval.val = subobjval;
		sortobjvals.push_back(indexval);
		subcplex.end();
		submodel.end();
		subx.end();
	}
	sort(sortobjvals.begin(), sortobjvals.end(), less<IndexVal>());
	return sortobjvals[knap.nbScens-mod.p-1].val;
}

double V2LPiter_nogood(IloEnv& env, const Knapsack& knap, const Model& mod, double UB, const vector< vector<int> >& exploredxsols)
{
	// UB could be initialized as the quantile bound
	double returnval;
	double U = UB;
	IloNumVar y(env, -IloInfinity, IloInfinity);
	IloNumVarArray u(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
	IloNumVarArray w(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
	IloNumVarArray x(env, knap.nbItems, 0, 1);
	IloNumVarArray z(env, knap.nbScens, 0, 1);
	IloModel model(env);
	// Nonanticipativity
	for (int j = 0; j < knap.nbItems; ++j)
	{
		for (int i = 0; i < knap.nbScens; ++i)
			model.add(u[i*knap.nbItems+j]+w[i*knap.nbItems+j]-x[j] == 0);
	}
	// objbound
	for (int i = 0; i < knap.nbScens; ++i)
	{
		IloExpr lhs(env);
		for (int j = 0; j < knap.nbItems; ++j)
			lhs += u[i*knap.nbItems+j]*knap.cost[j];
		lhs -= z[i]*U;
		lhs -= y;
		model.add(lhs >= -U);
		lhs.end();
		IloExpr lhs2(env);
		for (int j = 0; j < knap.nbItems; ++j)
			lhs2 += w[i*knap.nbItems+j]*knap.cost[j];
		lhs2 += z[i]*U;
		lhs2 -= y;
		model.add(lhs2 >= 0);
		lhs2.end();
	}
	// Au, Du, Dw
	for (int i = 0; i < knap.nbScens; ++i)
	{
		for (int k = 0; k < knap.nbDims; ++k)
		{
			IloExpr lhs(env);
			for (int j = 0; j < knap.nbItems; ++j)
				lhs += u[i*knap.nbItems+j]*knap.weight[i][k*knap.nbItems+j];
			lhs -= knap.capacity[k]*z[i];
			model.add(lhs <= 0);
			lhs.end();
		}
		for (int j = 0; j < knap.nbItems; ++j)
		{
			model.add(u[i*knap.nbItems+j]-z[i] <= 0);
			model.add(w[i*knap.nbItems+j]+z[i] <= 1);
		}
	
		// Du, Dw for no-good cuts
		for (int ii = 0; ii < exploredxsols.size(); ++ii)
		{
			IloExpr lhs(env);
			IloExpr lhs2(env);
			int cardpos = 0;
			for (int jj = 0; jj < knap.nbItems; ++jj)
			{	
				if (exploredxsols[ii][jj] == 0)
				{
					lhs += u[i*knap.nbItems+jj];
					lhs2 += w[i*knap.nbItems+jj];
				}
				else
				{
					lhs -= u[i*knap.nbItems+jj];
					lhs2 -= w[i*knap.nbItems+jj];
					cardpos++;
				}
			}
			lhs -= (1-cardpos)*z[i];
			lhs2 += (1-cardpos)*z[i];
			model.add(lhs >= 0);
			model.add(lhs2 >= 1-cardpos);
			lhs.end();
			lhs2.end();
		}	

	}

	// Knapsack
	model.add(IloSum(z) >= (knap.nbScens-mod.p));

	IloExpr obj(env);
	obj += y;
	model.add(IloMaximize(env, obj));
	obj.end();
	IloCplex cplex(model);
	cplex.setParam(IloCplex::SimDisplay, 0);

	cplex.solve();
	if (cplex.getStatus() == IloAlgorithm::Feasible || cplex.getStatus() == IloAlgorithm::Optimal)
		returnval = cplex.getObjValue();
	else
	{
		cout << "Infeasible!" << endl;
		exit(0);
	}
	cplex.end();
	model.end();
	u.end();
	w.end();
	x.end();
	z.end();

	return returnval;
}

double V2LP_nogood(IloEnv& env, const Knapsack& knap, const Model& mod, double UB, const vector< vector<int> >& exploredxsols)
{
	double U = UB;
	bool flag = 1;
	while (flag == 1)
	{
		double objval = V2LPiter_nogood(env, knap, mod, U, exploredxsols);
		cout << "objval = " << objval << endl;
		if (fabs((objval-U)*1.0/U) < 1e-5)
			flag = 0;
		else
			U = objval;
	}
	return U;
}

double Heuristic(const Knapsack& knap, const Model& mod)
{
	// create entireprobmodels
	double lb=0.0;
	double ub=knap.totalcost;
	double value=lb/2.0+ub/2.0;

	entireprob entireprobmodels;
	IloConstraint extra_constraint;

	IloEnv entireprobenv;
	entireprobmodels.env = entireprobenv;
	entireprobmodels.x=IloNumVarArray(entireprobenv, knap.nbItems, 0, 1);
	entireprobmodels.pi=IloNumVarArray(entireprobenv, knap.nbScens, 0, IloInfinity);

	entireprobmodels.model = IloModel(entireprobenv);
	IloExpr subobj(entireprobenv),subobjnew(entireprobenv);

	for (int j = 0; j < knap.nbItems; ++j)
	{
		subobj += entireprobmodels.x[j]*knap.cost[j];
		subobjnew-=entireprobmodels.x[j]*knap.cost[j];
	}


	for (int k = 0; k < knap.nbScens; ++k)
	{
		subobjnew+=entireprobmodels.pi[k];
		for (int d = 0; d < knap.nbDims; ++d)
		{
			IloExpr lhs(entireprobenv);			
			for (int j = 0; j < knap.nbItems; ++j)
				lhs += entireprobmodels.x[j]*knap.weight[k][d*knap.nbItems+j];
			lhs-=entireprobmodels.pi[k];
			entireprobmodels.model.add(lhs <= knap.capacity[d]);
			lhs.end();
		}
	}

	extra_constraint=(subobj>=value);

	entireprobmodels.model.add(extra_constraint);
	entireprobmodels.model.add(IloMinimize(entireprobenv, subobjnew));
	subobjnew.end();

	entireprobmodels.cplex = IloCplex(entireprobmodels.model);
	entireprobmodels.cplex.setParam(IloCplex::SimDisplay, 0);

	// Main loop
	double LB = 0;
	double UB = 1e8;
	while ((UB-LB)*1.0/UB > 1e-4)
	{
		int num_scenario;
		ub=min(UB,ub);
		lb=min(LB,lb);
		int iteration=0;
		while ((ub-lb)*1.0/ub > 1e-4)
		{
			vector<double>  xsols(knap.nbItems, 0);
			double obj = 0.0;
			vector<double>  pis;
			num_scenario=0;
			entireprobmodels.cplex.solve();
			if (entireprobmodels.cplex.getStatus() == IloAlgorithm::Feasible || entireprobmodels.cplex.getStatus() == IloAlgorithm::Optimal)
			{

				IloNumArray xvals(entireprobmodels.env);
				entireprobmodels.cplex.getValues(entireprobmodels.x, xvals);
				for (int j = 0; j < knap.nbItems; ++j)
				{

						xsols[j]=xvals[j];
						obj+=knap.cost[j]*xvals[j];
				}

				xvals.end();

				for (int k = 0; k < knap.nbScens; ++k)
				{			

					double pi1=entireprobmodels.cplex.getValue(entireprobmodels.pi[k]);				

					if (pi1 > 1e-7)
					{
						num_scenario++;
					}					
				}
			}
			else
			{
				cout << "Subproblems "  << " are infeasible!" << endl;
				exit(0);
			}

			entireprobmodels.model.remove(extra_constraint);
			ub=min(UB,ub);
			if(num_scenario>mod.p)
			{
				ub=value;
				value=lb/2.0+ub/2.0;
			}else
			{
				LB=obj;
				lb=max(value,LB);
				value=lb/2.0+ub/2.0;
			}

			extra_constraint=(subobj>=value);
			entireprobmodels.model.add(extra_constraint);

			iteration++;
		}

		UB=min(UB,ub);
	}
	return LB;
}


double V2IP(IloEnv& env, const Knapsack& knap, Model& mod, double UB, int option)
{
	// UB could be initialized as the quantile bound
	double returnval;
	double U = UB;
	Model2 modlocal;
	modlocal.y = IloNumVar(env, -IloInfinity, IloInfinity);
	modlocal.u = IloNumVarArray(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
	modlocal.w = IloNumVarArray(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
	modlocal.x = IloNumVarArray(env, knap.nbItems, 0, 1);
	modlocal.z = IloIntVarArray(env, knap.nbScens, 0, 1);
	IloModel model(env);
	// Nonanticipativity
	for (int j = 0; j < knap.nbItems; ++j)
	{
		for (int i = 0; i < knap.nbScens; ++i)
			model.add(modlocal.u[i*knap.nbItems+j]+modlocal.w[i*knap.nbItems+j]-modlocal.x[j] == 0);
	}
	// objbound
	for (int i = 0; i < knap.nbScens; ++i)
	{
		IloExpr lhs(env);
		for (int j = 0; j < knap.nbItems; ++j)
			lhs += modlocal.u[i*knap.nbItems+j]*knap.cost[j];
		lhs -= modlocal.z[i]*U;
		lhs -= modlocal.y;
		model.add(lhs >= -U);
		lhs.end();
		IloExpr lhs2(env);
		for (int j = 0; j < knap.nbItems; ++j)
			lhs2 += modlocal.w[i*knap.nbItems+j]*knap.cost[j];
		lhs2 += modlocal.z[i]*U;
		lhs2 -= modlocal.y;
		model.add(lhs2 >= 0);
		lhs2.end();
	}
	// Au, Du, Dw
	for (int i = 0; i < knap.nbScens; ++i)
	{
		for (int k = 0; k < knap.nbDims; ++k)
		{
			IloExpr lhs(env);
			for (int j = 0; j < knap.nbItems; ++j)
				lhs += modlocal.u[i*knap.nbItems+j]*knap.weight[i][k*knap.nbItems+j];
			lhs -= knap.capacity[k]*modlocal.z[i];
			model.add(lhs <= 0);
			lhs.end();
		}
		for (int j = 0; j < knap.nbItems; ++j)
		{
			model.add(modlocal.u[i*knap.nbItems+j]-modlocal.z[i] <= 0);
			model.add(modlocal.w[i*knap.nbItems+j]+modlocal.z[i] <= 1);
		}
	}
	// Knapsack
	model.add(IloSum(modlocal.z) >= (knap.nbScens-mod.p));

	IloExpr obj(env);
	obj += modlocal.y;
	model.add(IloMaximize(env, obj));
	obj.end();
	IloCplex cplex(model);
	cplex.setParam(IloCplex::MIPDisplay, 0);
	if (option == 1)
	{
		mod.lastnode = 0;
		cplex.use(new (env) BranchCallback(env));
		cplex.use(new (env) UserCutCallback(env, modlocal, knap, mod));
		// Will use local cuts
	}
	cplex.solve();
	if (cplex.getStatus() == IloAlgorithm::Feasible || cplex.getStatus() == IloAlgorithm::Optimal)
	{
		returnval = cplex.getObjValue();
		if (option == 0)
			mod.V2IPNnodes = cplex.getNnodes();
		if (option == 1)
			mod.V2IPlocalNnodes = cplex.getNnodes();
	}
	else
	{
		cout << "Infeasible!" << endl;
		exit(0);
	}
	cplex.end();
	model.end();
	modlocal.u.end();
	modlocal.w.end();
	modlocal.x.end();
	modlocal.z.end();
	return returnval;
}

double quantile_fixed(const Knapsack& knap, const Model& mod, const IloNumArray& ubs, const IloNumArray& lbs)
{
	vector<IndexVal> sortobjvals;
	int fix_one = 0;
	double min_fix_one = 1e8; // record the smallest objval for scenarios fixed to 1
	for (int j = 0; j < knap.nbScens; ++j)
	{
		if (ubs[j] > 1-1e-5 && lbs[j] > 1-1e-5)
		{
			fix_one++;
			if (mod.objvals[j] < min_fix_one) 
				min_fix_one = mod.objvals[j];
		}
		else
		{
			if (ubs[j] < 1e-5 && lbs[j] < 1e-5)
				continue; // don't consider ones that are fixed to 0
			else
			{
				IndexVal indexval;
				indexval.ind = j;
				indexval.val = mod.objvals[j];
				sortobjvals.push_back(indexval);
			}
		}
	}
	sort(sortobjvals.begin(), sortobjvals.end(), less<IndexVal>());
	double candidate = sortobjvals[knap.nbScens-mod.p-1-fix_one].val;
	if (min_fix_one < candidate)
		return min_fix_one;
	else
		return candidate;
}
