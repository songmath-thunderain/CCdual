/*- mode: C++;
 * Date: Jan 28, 2015
 */

#include <iostream>
#include <ilcplex/ilocplex.h>
#include "ccdual-decomp.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <math.h>
using namespace std;

ILOMIPINFOCALLBACK2(InfoCallback, Model&, mod, IloTimer&, timer) 
{
	if (getNnodes() == 0) 
        mod.rootlb = getBestObjValue(); 
	double U = getBestObjValue();
	double L = getIncumbentObjValue();
	if ((U-L)*1.0/L < 0.005 && mod.fiveflag == 0)
	{
		mod.fiveflag = 1;
		mod.fivetime = timer.getTime();
	}
	if ((U-L)*1.0/L < 0.001 && mod.oneflag == 0)
	{
		mod.oneflag = 1;
		mod.onetime = timer.getTime();
	}
}

class DepthData : public IloCplex::MIPCallbackI::NodeData {
   unsigned int const depth; /**< Depth of node at which this data
                                 *   is stored. */
   public:
      DepthData(unsigned int idepth) : depth(idepth) {}
      unsigned int getDepth() const { return depth; }
};
 
DepthData const rootDepth(0);

struct BranchCallback : public IloCplex::BranchCallbackI {

	const Knapsack& knap;
	Model& mod;

   BranchCallback(IloEnv env, const Knapsack& knapee,Model& modee) : IloCplex::BranchCallbackI(env), knap(knapee), mod(modee) {}
 
   IloCplex::CallbackI *duplicateCallback() const {
         return new (getEnv()) BranchCallback(*this);
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
				vars.end();
				dirs.end();
           	bounds.end();
         }
         else {
            // CPLEX would not create any branch here. Prune the node.
            prune();
         }
   }
};

struct UserCutCallback : public IloCplex::UserCutCallbackI {
	const Knapsack& knap;
	Model& mod;
	IloTimer& timer;
	public:
	UserCutCallback(IloEnv env,const Knapsack& knapee, Model& modee, IloTimer& timee) : IloCplex::UserCutCallbackI(env), knap(knapee), mod(modee), timer(timee) {	
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
		double U = getBestObjValue();
		if (((nNodes > 0 && nNodes != mod.lastnode && (depth < 10 || depth %15==0) ) || nNodes == 0) && (mod.solveoption == 2 || mod.solveoption == 3))
		{
			// at non-root node, just one round of cut generation.
			// also limit to depth<10, or depth is a multiple of 15
			IloNumArray xvals(env);
			getValues(xvals, mod.x);
			for (int j = 0; j < knap.nbItems; ++j)
			{
				if (xvals[j] < 1e-5)
					xvals[j] = 0;
				if (xvals[j] > 1-1e-5)
					xvals[j] = 1;
			}
			IloNumArray zvals(env);
			getValues(zvals, mod.z);
			double yval;
		    if (mod.solveoption	== 2)
				yval = getValue(mod.y);
			else
				yval = 0;
			vector<BendersCoef> temp;
			int option;
			if (depth >= 4)
				option = 1;
			else
				option = 0;
			// Do less intensive version (check only large enough z) after depth 4
			double start = timer.getTime();
			genCuts(env, knap, mod, xvals, zvals, yval, U, temp, option, mod.cutviol);		
			mod.subLPtime += timer.getTime()-start;
			vector<IndexVal> cutlist;
			for (int k = 0; k < temp.size(); ++k)
			{
				IndexVal indexval;
				indexval.ind = k;
				indexval.val = temp[k].normviol;
				cutlist.push_back(indexval);
			}
			sort(cutlist.begin(), cutlist.end(), less<IndexVal>());
			vector<BendersCoef> added;
			vector<double> norm_added;
			int nadded = 0;
			for (int l = 0; l < cutlist.size(); ++l)
			{
				// For all cuts that have been collected, check if they are too parallel to each other and only add diverse cuts
				int k = cutlist[l].ind;
				BendersCoef bcoef = temp[k];
				if (check_parallel(mod, added, bcoef, norm_added) == 1)
				{
					nadded++;
					IloExpr lhs(env);
					//double tempnorm = 0;
					for (int j = 0; j < knap.nbItems; ++j)
					{
						//tempnorm += bcoef.xcoef[j]*bcoef.xcoef[j];
						lhs += mod.x[j]*bcoef.xcoef[j];
					}
					lhs += mod.y*bcoef.ycoef;
					lhs += mod.z[bcoef.zind]*(bcoef.zUcoef*U+bcoef.zrest);
					add(lhs <= (bcoef.rhsUcoef*U+bcoef.rhsrest));
					lhs.end();
					//added.push_back(bcoef);
					//norm_added.push_back(tempnorm);
				}
				if ((depth < 4 && nadded > knap.nbItems)||(depth >= 5 && nadded > 0)) break;
			}
			mod.nUsercuts += nadded;
			xvals.end();
			zvals.end();
		}
		mod.lastnode = getNnodes();
	}
};

struct LazyConstraintCallback : public IloCplex::LazyConstraintCallbackI {
	const Knapsack& knap;
	Model& mod;
	IloTimer& timer;
	public:
	LazyConstraintCallback(IloEnv env,const Knapsack& knapee, Model& modee, IloTimer& timee) : IloCplex::LazyConstraintCallbackI(env), knap(knapee), mod(modee), timer(timee) {	
	}
	IloCplex::CallbackI *duplicateCallback() const {
	   return new (getEnv()) LazyConstraintCallback(*this);
	}

	void main() {
		// start cutcallback
		DepthData const *data = dynamic_cast<DepthData *>(getNodeData());
        if ( !data )
           data = &rootDepth;
		int depth = data->getDepth();
		if (mod.solveoption == 2 || mod.solveoption == 3)
		{
			IloEnv env = getEnv();
			int nNodes = getNnodes();
			IloNumArray xvals(env);
			getValues(xvals, mod.x);
			for (int j = 0; j < knap.nbItems; ++j)
			{
				if (xvals[j] < 1e-5)
					xvals[j] = 0;
				if (xvals[j] > 1-1e-5)
					xvals[j] = 1;
			}
			IloNumArray zvals(env);
			getValues(zvals, mod.z);
			double yval;
		    if (mod.solveoption	== 2)
				yval = getValue(mod.y);
			else
				yval = 0;
			double start = timer.getTime();
			double U = getBestObjValue();
			vector<BendersCoef> temp;
			genCuts(env, knap, mod, xvals, zvals, yval, U, temp, 0, 1e-5);		
			mod.subLPtime += timer.getTime()-start;
			vector<IndexVal> cutlist;
			for (int k = 0; k < temp.size(); ++k)
			{
				IndexVal indexval;
				indexval.ind = k;
				indexval.val = temp[k].normviol;
				cutlist.push_back(indexval);
			}
			sort(cutlist.begin(), cutlist.end(), less<IndexVal>());
			vector<BendersCoef> added;
			vector<double> norm_added;
			int nadded = 0;
			for (int l = 0; l < cutlist.size(); ++l)
			{
				// Going from the mostly violated cut to the least, and check parallelism
				int k = cutlist[l].ind;
				BendersCoef bcoef = temp[k];
				if (check_parallel(mod, added, bcoef, norm_added) == 1)
				{
					nadded++;
					IloExpr lhs(env);
					//double tempnorm = 0;
					for (int j = 0; j < knap.nbItems; ++j)
					{
						//tempnorm += bcoef.xcoef[j]*bcoef.xcoef[j];
						lhs += mod.x[j]*bcoef.xcoef[j];
					}
					lhs += mod.y*bcoef.ycoef;
					lhs += mod.z[bcoef.zind]*(bcoef.zUcoef*U+bcoef.zrest);
					add(lhs <= (bcoef.rhsUcoef*U+bcoef.rhsrest));
					lhs.end();
					//added.push_back(bcoef);
					//norm_added.push_back(tempnorm);
				}
				if (nadded > knap.nbItems) break;
			}
			mod.nLazycuts += nadded;
			xvals.end();
			zvals.end();
			mod.lastnode = getNnodes();
		}
	}
};


int main(int argc, char **argv)
{
	// Run program like this:
	// ./ccdual-decomp instancefile eps option solveoption subproboption cutviol resultfile
	// ./ccdual-decomp instancefile 0.1 1 2 2 0.001 resultfile
	// Solving MIP formulation with a given lower bound L
    IloEnv env;
	IloEnv env2;
    Knapsack knap;
	Model mod;
	mod.env = env2;
	mod.option = atoi(argv[3]);
	// option = 1: solve CCLP by v2LP using the LP relaxation
	// option = 2: solve CCIP by v2LP using the LP relaxation
	// option = 3: same as option 1, but just for set packing instances
	// option = 4: same as option 1, but just for set covering instances
	mod.solveoption = atoi(argv[4]);
	// solveoption = -1: solve big-M formulation
	// solveoption = 1: solve the extended formulation with a lower bound l
	// solveoption = 2: solve the extended formulation with Benders
	// solveoption = 3: solve the first extended formulation with Benders
	mod.subproboption = atoi(argv[5]);
	// subproboption = 1: use naive big-M values, and just use bounds in the z=0 case
	// subproboption = 2: use improved big-M values, and use these as constraints for that scenario in z=0 case
	mod.cutviol = atof(argv[6]);
    knap.weight = NumMatrix(env);
    knap.cost = IloNumArray(env);
    knap.capacity = IloNumArray(env);
	double eps = atof(argv[2]);
    const char* filename = NULL;
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
	if (mod.option == 1 || mod.option == 2)
	{
    	file >> knap.cost >> knap.capacity >> knap.weight;
    	knap.nbItems = knap.cost.getSize();
		knap.nbDims = knap.capacity.getSize();
    	knap.nbScens = knap.weight.getSize();
	}
	if (mod.option == 3)
	{
		// set packing mode
		NumMatrix temp(env);
		NumMatrix temp2(env);
		NumMatrix temp3(env);
		file >> knap.cost >> temp >> temp2 >> temp3;
		knap.nbItems = knap.cost.getSize();
		knap.nbDims = temp2.getSize();
		knap.nbScens = temp3.getSize();
		for (int i = 0; i < knap.nbScens; ++i)
		{
			IloNumArray scen(env, knap.nbDims*knap.nbItems);
			for (int jj = 0; jj < knap.nbDims*knap.nbItems; ++jj)
				scen[jj] = 0;
			for (int d = 0; d < knap.nbDims; ++d)
			{
				for (int j = 0; j < temp2[d].getSize(); ++j)
				{
					int ind = temp2[d][j];
					if (temp3[i][ind] == 1)
						scen[d*knap.nbItems+ind] = 1;
				}
			}
			knap.weight.add(scen);
		}
		knap.capacity = IloNumArray(env, knap.nbDims);
		for (int d = 0; d < knap.nbDims; ++d)
			knap.capacity[d] = 1;
		temp.end();
		temp2.end();
		temp3.end();
	}
	if (mod.option == 4)
	{
		// set covering mode
		NumMatrix temp(env);
		file >> knap.cost >> temp;
		knap.nbItems = knap.cost.getSize();
		knap.nbDims = 1;
		knap.nbScens = temp.getSize();
		for (int j = 0; j < knap.nbItems; ++j)
			knap.cost[j] = -knap.cost[j]*100;
		for (int i = 0; i < knap.nbScens; ++i)
		{
			IloNumArray scen(env, knap.nbItems);
			for (int jj = 0; jj < knap.nbItems; ++jj)
				scen[jj] = -temp[i][jj];
			knap.weight.add(scen);
		}
		knap.capacity = IloNumArray(env, knap.nbDims);
		for (int d = 0; d < knap.nbDims; ++d)
			knap.capacity[d] = -1;
		temp.end();
	}
	
	mod.p = int(eps*knap.nbScens+1e-5);
	mod.lastnode = 0;
	mod.ncuts = 0;
	mod.maxpar = 0.999;
	mod.zthreshold = 1e-3;
	//mod.zthreshold=0.1;
	mod.nUsercuts = 0;
	mod.nLazycuts = 0;
	mod.fiveflag = 0;
	mod.oneflag = 0;
	mod.fivetime = 0;
	mod.onetime = 0;
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

	if (mod.option == 1 || mod.option == 2 || mod.option == 3 || mod.option == 4)
	{
		/*
			Build MIP model based on primal formulation for quantile based bound
		*/
		if (mod.option == 1 || mod.option == 3 || mod.option == 4)
			mod.x = IloNumVarArray(env, knap.nbItems, 0, 1);
		if (mod.option == 2) {
			mod.x = IloNumVarArray(env, knap.nbItems, 0, 1, ILOINT);
		}
    	mod.z = IloNumVarArray(env, knap.nbScens, 0, 1);
	   	mod.y = IloNumVar(env, -IloInfinity, IloInfinity);
		IloModel model(env);
 		IloCplex cplex(model);

    	// This is equivalent to IloCplex cplex(env) + cplex.extract(model)
		if (mod.solveoption == -1) { 
			double buildstart = clock.getTime();
			buildMIP(env, model, knap, mod, 1);
			buildtime = clock.getTime() - buildstart;
			cplex.use(new (env) UserCutCallback(env, knap, mod, clock));
			cplex.use(new (env) LazyConstraintCallback(env, knap, mod, clock));
			cplex.setParam(IloCplex::PreInd, CPX_OFF);
			cout << "build time (includes bigm calcs) = " << buildtime << endl;
		}
			
		if (mod.solveoption == 1 || mod.solveoption == 2 || mod.solveoption == 3)
		{
			double buildstart = clock.getTime();
			buildMIP(env, model, knap, mod, 0);
			if (mod.solveoption == 1)
			{
				mod.u = IloNumVarArray(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
				mod.w = IloNumVarArray(env, knap.nbScens*knap.nbItems, 0, IloInfinity);
			}
			buildtime = clock.getTime() - buildstart;
			cout << "build time (includes bigm calcs) = " << buildtime << endl;
			double Ustart = clock.getTime();
			// Solve the root LP, use V2LP bound as U
			double quantileval = quantile(env, knap, mod);
			/*
			int iterV2 = 0;
			double U = V2LP(env,knap, mod, quantileval, iterV2);
			int iterV2 = 0;
			double U = V3LD(env, knap, mod, quantileval, 0, iterV2);
			*/
			double U;
			if (mod.solveoption == 1)
				U = V2LPiter2(env, knap, mod, quantileval);
			if (mod.solveoption == 2 || mod.solveoption == 3)
			{
				// create the subproblems
				for (int i = 0; i < knap.nbScens; ++i)
				{
					Subprob subprob;
					if (mod.fixscen[i] == 0)
					{
						IloEnv subenv;
						subprob.env = subenv;
						subprob.alpha = IloNumVarArray(subprob.env, knap.nbItems, -IloInfinity, IloInfinity);
						subprob.lambda = IloNumVarArray(subprob.env, 2, 0, IloInfinity);
						if (mod.solveoption == 3) {
							subprob.lambda[0].setUB(0.0);
							subprob.lambda[1].setUB(0.0);
						}
						subprob.gamma1 = IloNumVarArray(subprob.env, knap.nbDims, 0, IloInfinity);
						subprob.gamma2 = IloNumVarArray(subprob.env, knap.nbDims, 0, IloInfinity);
						subprob.beta1 = IloNumVarArray(subprob.env, knap.nbItems, 0, IloInfinity);
						subprob.beta2 = IloNumVarArray(subprob.env, knap.nbItems, 0, IloInfinity);
						//subprob.xi = IloNumVarArray(subprob.env, knap.nbItems, -IloInfinity, IloInfinity);
						subprob.model = IloModel(subprob.env);
						subprob.cplex = IloCplex(subprob.model);
						subprob.obj = IloMaximize(subprob.env);
						subprob.model.add(subprob.obj);
					}
					mod.subprobs.push_back(subprob);
				}
				buildSub(knap, mod);
				// Solve a rootLP, and add all binding constraints into the formulation
				U = rootLP(env, cplex, model, knap, mod, quantileval);
				cout << "u = " << U << endl;
			}
			cout << "root U = " << U << endl;
			mod.U = U;
			mod.Utime = clock.getTime()-Ustart;
			cout << "root U time = " << mod.Utime << endl;
			if (mod.solveoption == 1)
				V2IP(env, model, knap, mod, mod.U);
		}
		/*
		if (mod.option == 2) {
			IloNumArray pri(env);
			pri.add(knap.nbItems,1);
			cplex.setPriorities(mod.x,pri);	
		}
		*/
		cplex.use(InfoCallback(env, mod, clock));

		if (mod.option != 2)
		{
			// Use heuristic as MIPstart
			vector<double> xsol(knap.nbItems, 0);
			double heurstart = clock.getTime();
			double heurobj = Heuristic(knap, mod, xsol);
			mod.heurtime = clock.getTime()-heurstart;
			cout << "heurobj = " << heurobj << endl;
			cout << "heutime  = " << mod.heurtime << endl;
			vector<bool> zvals(knap.nbScens, 0);
			int succ = 0;
			for (int i = 0; i < knap.nbScens; ++i)
			{
				zvals[i] = checkFeas(i, xsol, knap);
				if (zvals[i] == 1)
					succ++;
			}
			if (succ < knap.nbScens-mod.p)
				cout << "Something wrong with heuristic!" << " succ = " << succ << endl;
			else
			{
				IloNumVarArray startVar(env);
				IloNumArray startVal(env);
				if (mod.solveoption == -1 || mod.solveoption == 3)
					MIPstart(startVar, startVal, knap, mod, xsol, zvals, heurobj, 0);
				if (mod.solveoption == 1)
					MIPstart(startVar, startVal, knap, mod, xsol, zvals, heurobj, 1);
				if (mod.solveoption == 2)
					MIPstart(startVar, startVal, knap, mod, xsol, zvals, heurobj, 2);
				try{
					cplex.addMIPStart(startVar,startVal,IloCplex::MIPStartCheckFeas);
				}
				
     			catch (IloException& e) {
					cerr << "Concert exception caught: " << e << endl;
     			}
			}
		}
		if (mod.solveoption == 2 || mod.solveoption == 3)
		{
			// Solve by feasibility cuts that project out u and w
			// construct a feasibility subprob LP
			mod.subLPtime = 0;
			if (mod.solveoption == 2)
			{
				model.add(mod.y <= mod.U);
				IloExpr lhs(env);
				lhs += mod.y;
				for (int j = 0; j < knap.nbItems; ++j)
					lhs -= mod.x[j]*knap.cost[j];
				model.add(lhs <= 0);
				lhs.end();
			}
			/*
			// NEW! Add bigM inequalities
			if (mod.subproboption == 2)
			{
				for (int i = 0; i < knap.nbScens; ++i)
				{
					for (int k = 0; k < knap.nbDims; ++k)
					{
						IloExpr lhs(env);
						for (int j = 0; j < knap.nbItems; ++j)
							lhs += mod.x[j]*knap.weight[i][k*knap.nbItems+j];
						model.add(lhs <= (knap.capacity[k]+mod.bigM[k][i]));
						lhs.end();
					}
				}
			}
			*/
			cplex.use(new (env) BranchCallback(env, knap, mod));
			cplex.use(new (env) UserCutCallback(env, knap, mod, clock));
			cplex.use(new (env) LazyConstraintCallback(env, knap, mod, clock));
			cplex.setParam(IloCplex::PreInd, CPX_OFF);
			// MIP emphasis on optimality
			//cplex.setParam(IloCplex::MIPEmphasis, 2);
		}
		
		cout << "Model size: Rows = " << cplex.getNrows() << ", Cols = " << cplex.getNcols() << endl;
  	  	cplex.setParam(IloCplex::TiLim,3600);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::MIPDisplay, 4);
		double startcplex = clock.getTime();
		cplex.solve();	
		double cplextime = clock.getTime()-startcplex;
		if (mod.fivetime > 1e-5)
			mod.fivetime -= startcplex;
		if (mod.onetime > 1e-5)
			mod.onetime -= startcplex;
		vector<double> xsol(knap.nbItems, 0);
		IloNumArray xsolvals(env);
		cplex.getValues(xsolvals, mod.x);
		IloNumArray zsolvals(env);
		cplex.getValues(zsolvals, mod.z);
		for (int j = 0; j < knap.nbItems; ++j)
			xsol[j] = xsolvals[j];
		vector<bool> zvals(knap.nbScens, 0);
		int succ = 0;
		for (int i = 0; i < knap.nbScens; ++i)
		{
			zvals[i] = checkFeas(i, xsol, knap);
			if (zvals[i] != zsolvals[i])
				cout << "scenario " << i << ": zsol = " << zsolvals[i] << ", actual = " << zvals[i] << endl;
		}
		if (cplex.getStatus() != IloAlgorithm::Feasible && cplex.getStatus() != IloAlgorithm::Optimal)
		{
			cout << "Main problem infeasible!" << endl;
			exit(0);
		}
		double cplexobj = cplex.getObjValue();
		int cplexnode = cplex.getNnodes();
		cout << "cplexobj = " << cplexobj << endl;
		cout << "fivetime = " << mod.fivetime << endl;
		cout << "onetime = " << mod.onetime << endl;
		cout << "cplextime = " << cplextime << endl;
		cout << "cplexnode = " << cplexnode << endl;
		cout << "heurtime = " << mod.heurtime << endl;
		if (mod.solveoption == 1 || mod.solveoption == 2 || mod.solveoption == 3)
		{
			cout << "mod.U = " << mod.U << endl;
			cout << "mod.Utime = " << mod.Utime << endl;
		}
		if (mod.solveoption == 2 || mod.solveoption == 3)
		{
			cout << "subLPtime = " << mod.subLPtime << endl;
			cout << "# user cuts generated = " << mod.nUsercuts << endl;
			cout << "# lazy cuts generated = " << mod.nLazycuts << endl;
			mod.ncuts = mod.nLazycuts + mod.nUsercuts;
			cout << "# cut generated = " << mod.ncuts << endl;
		}
		ofstream out(argv[7], ios::app);
		if (out)  {
				out << setw(24) << argv[1];
				out << ", ";
				out << setw(8) << argv[2];
				out << ", ";
				out << setw(8) << argv[3];
				out << ", ";
				out << setw(8) << argv[4];
				out << ", ";
				out << setw(8) << argv[5];
				out << ", ";
				out << setw(8) << cplextime;
				out << ", ";
				out << setw(8) << cplexnode;
				out << ", ";
				out << setw(8) << cplex.getObjValue();
				out << ", ";
				out << setw(8) << cplex.getBestObjValue();
				out << ", ";
				out << setw(8) << mod.rootlb;
				out << ", ";
				out << setw(8) << mod.heurtime;
				if (mod.solveoption == 1 || mod.solveoption == 2 || mod.solveoption == 3)
				{
					out << ", ";
					out << setw(8) << mod.U;
					out << ", ";
					out << setw(8) << mod.Utime;
				}
				if (mod.solveoption == 2 || mod.solveoption == 3)
				{
					out << ", ";
					out << setw(8) << mod.subLPtime;
					out << ", ";
					out << setw(8) << mod.ncuts;
					out << ", ";
					out << setw(8) << mod.nUsercuts;
					out << ", ";
					out << setw(8) << mod.nLazycuts;
				}
				if (mod.solveoption == -1 || mod.solveoption == 2 || mod.solveoption == 3)
				{
					out << ", ";
					out << setw(8) << mod.fivetime;
					out << ", ";
					out << setw(8) << mod.onetime;
				}
				out << endl;
				out.close();
		}	

		cplex.end();
		model.end();
		mod.x.end();
		mod.z.end();
		if (mod.solveoption == 1)
		{
			mod.y.end();
			mod.u.end();
			mod.w.end();
		}
		if (mod.solveoption == 2 || mod.solveoption == 3)
		{
			for (int i = 0; i < knap.nbScens; ++i)
			{
				if (mod.fixscen[i] == 0)
				{
					mod.subprobs[i].cplex.end();
					mod.subprobs[i].model.end();
					mod.subprobs[i].alpha.end();
					mod.subprobs[i].lambda.end();
					mod.subprobs[i].gamma1.end();
					mod.subprobs[i].gamma2.end();
					mod.subprobs[i].beta1.end();
					mod.subprobs[i].beta2.end();
					mod.subprobs[i].obj.end();
					mod.subprobs[i].env.end();
				}
			}
		}
	}
	env.end();
	env2.end();
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

double calcBigM(const Knapsack& knap, const Model& mod, int row, int scen) {

	// use the given row/scen as objective
	// Solve one-row knapsack LP for every other scenario and matching row
	// Sort and return p largest value as a valid big-m

	vector<IndexVal> optvals;
	
	for (int i = 0; i < knap.nbScens; ++i) {
		double curobj = 0.0;

		if (mod.option != 4) {
			vector<IndexVal> ratiovals;

			for (int j = 0; j < knap.nbItems; j++) {
				IndexVal indexval;
				indexval.ind = j;
				if (knap.weight[scen][row*knap.nbItems+j] < 1e-6) 
					indexval.val = 0.0;
				else
					indexval.val = knap.weight[scen][row*knap.nbItems+j]/knap.weight[i][row*knap.nbItems+j];
				ratiovals.push_back(indexval);
			}
			sort(ratiovals.begin(), ratiovals.end(), less<IndexVal>());
			double curweight = 0.0;
			int curind = 0;
			int curit = ratiovals[0].ind;
			while (curweight + knap.weight[i][row*knap.nbItems+curit] <= knap.capacity[row]) {
				curweight += knap.weight[i][row*knap.nbItems+curit];
				curobj += knap.weight[scen][row*knap.nbItems+curit];
				curind++;
				if (curind == knap.nbItems) break;
				curit = ratiovals[curind].ind;
			}
			if (curind < knap.nbItems) {
				double lastval = (knap.capacity[row] - curweight)/knap.weight[i][row*knap.nbItems+curit]; 	
				if (lastval < -0.000001 || lastval > 1.000001) {
					cout << "Error in continuous knapsack calculation!" << endl;
					cout << "cur it weight = " <<
					knap.weight[i][row*knap.nbItems+curit] << endl;
					cout << "cur it obj = " <<
					knap.weight[scen][row*knap.nbItems+curit] << endl;
					cout << "last val = " << lastval << endl; 
				}
				curobj += knap.weight[scen][row*knap.nbItems+curit]*lastval;		
			}
		}
		else {
			// both obj coef and constraint coefs are negative, find lowest ratio	
			int minind = -1;
			double minratio = 0.0;
			for (int j = 0; j < knap.nbItems; j++) {
				double curratio = knap.weight[scen][row*knap.nbItems+j]/knap.weight[i][row*knap.nbItems+j]; 
				if (minind < 0 || curratio < minratio) {
					minind = j;
					minratio = curratio;
				}
			}
			// solution sets x[minind] = rhs/coef[minind]
			curobj = knap.weight[scen][row*knap.nbItems+minind]*knap.capacity[row]/knap.weight[i][row*knap.nbItems+minind];
		}
		IndexVal objindval;
		objindval.ind = i;
		objindval.val = curobj;
		optvals.push_back(objindval);	
	}
	sort(optvals.begin(), optvals.end(), greater<IndexVal>());
	
	return optvals[mod.p+1].val;

}

void buildMIP(IloEnv& env, IloModel& model, Knapsack& knap, Model& mod, int option)
{
	// option = 0: do nothing but preprocessing
	// option = 1: add bigM inequalities
	model.add(IloSum(mod.z) >= (knap.nbScens-mod.p));
	cout << "building model " << endl;

	for (int k = 0; k < knap.nbDims; ++k)
	{
		vector<double> subM(knap.nbScens, 0);
		mod.bigM.push_back(subM);
	}

	if (mod.subproboption == 2)
	{
		cout << "subproboption = 2" << endl;
		// Strengthened big-M constraint
		vector< vector< vector<double> > > total_pool;
		// Initial big-M constraint
		for (int i = 0; i < knap.nbScens; ++i)
		{	
			if (i%50 == 0) cout << "done with i = " << i << " out of " << knap.nbScens << endl;
			vector< vector<double> > pool;
			for (int k = 0; k < knap.nbDims; ++k)
			{
				vector<double> subpool(knap.nbScens, 0.0);
				pool.push_back(subpool);
			}
			single_sort(knap, mod, i, pool);
			total_pool.push_back(pool);	
		}

		// Add big-M valid inequalities for (x,z) space
		for (int i = 0; i < knap.nbScens; ++i)
		{
			for (int k = 0; k < knap.nbDims; ++k)
				mod.bigM[k][i] = PPlus(knap, mod, i, k, total_pool);
		}
	}
	if (mod.subproboption == 1)
	{
		// calculate naive bigM values
		for (int i = 0; i < knap.nbScens; ++i)
		{
			for (int k = 0; k < knap.nbDims; k++)
			{
				double totalweight = 0.0;
				for (int j = 0; j < knap.nbItems; j++)
					totalweight += knap.weight[i][k*knap.nbItems+j];
				double bigM = totalweight - knap.capacity[k];
				mod.bigM[k][i] = bigM;
			}	
		}
	}
	mod.fixscen = vector<bool>(knap.nbScens, 0);
	for (int i = 0; i < knap.nbScens; ++i)
	{
		bool fflag = 1;
    	for (int k = 0; k < knap.nbDims; k++)
    	{
			double bigM = mod.bigM[k][i]; 
			if (bigM > 1e-5)
			{
				fflag = 0;
				if (option == 1) {
					// Add bigM inequalities
					IloExpr lhs(env);
					for (int j = 0; j < knap.nbItems; ++j)
						lhs += mod.x[j]*knap.weight[i][k*knap.nbItems+j];
					lhs += mod.z[i]*bigM;
					model.add(lhs <= (knap.capacity[k]+bigM));
					lhs.end();
				}
			}
    	}
		if (fflag == 1)
			mod.fixscen[i] = 1;
	}
    // big-M constraint
    for (int i = 0; i < knap.nbScens; i++)
 	{
		if (mod.fixscen[i] == 1)
			model.add(mod.z[i] >= 1);
	}	

	IloExpr obj(env);
	if (mod.solveoption == 1 || mod.solveoption == 2)
	{
		// max y
		obj += mod.y;
	}
	if (mod.solveoption == -1 || mod.solveoption == 3) 
	{
		// The first extended formulation
		for (int j = 0; j < knap.nbItems; ++j)
			obj += mod.x[j]*knap.cost[j];
	}
	model.add(IloMaximize(env, obj));
	obj.end();
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
		subcplex.setParam(IloCplex::Threads, 1);
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
	//for (int i = 0; i < mod.p; i++) {
	//	cout << "scenval[" << i << "] = " << sortobjvals[knap.nbScens-i-1].val << endl;
	//}
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
	cplex.setParam(IloCplex::Threads, 1);
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


double V2LPiter2(IloEnv& env, const Knapsack& knap, const Model& mod, double UB)
{
	// Do iteration inside, so don't have to create/cleanup at every iteration
	// UB could be initialized as the quantile bound
	double Uupdate = 1e8;
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
	IloRangeArray objcons1(env);
	IloRangeArray objcons2(env);
	for (int i = 0; i < knap.nbScens; ++i)
	{
		IloExpr lhs(env);
		for (int j = 0; j < knap.nbItems; ++j)
			lhs += u[i*knap.nbItems+j]*knap.cost[j];
		lhs -= z[i]*U;
		lhs -= y;
		IloRange range(env, -U, lhs, IloInfinity);
		model.add(range);
		objcons1.add(range);
		lhs.end();
		IloExpr lhs2(env);
		for (int j = 0; j < knap.nbItems; ++j)
			lhs2 += w[i*knap.nbItems+j]*knap.cost[j];
		lhs2 += z[i]*U;
		lhs2 -= y;
		IloRange range2(env, 0, lhs2, IloInfinity);
		model.add(range2);
		objcons2.add(range2);
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
	cplex.setParam(IloCplex::Threads, 1);
	bool flag = 1;
	while (flag == 1)
	{
		cplex.solve();
		if (cplex.getStatus() == IloAlgorithm::Feasible || cplex.getStatus() == IloAlgorithm::Optimal)
			Uupdate = cplex.getObjValue();
		else
		{
			cout << "Infeasible!" << endl;
			exit(0);
		}
		if (fabs((Uupdate-U)*1.0/U) <= 1e-5)
		{
			flag = 0;
			break;
		}
		else
		{
			// update the model and resolve
			U = Uupdate;
			for (int i = 0; i < knap.nbScens; ++i)
			{
				objcons1[i].setLinearCoef(z[i], -U);
				objcons2[i].setLinearCoef(z[i], U);
				objcons1[i].setLB(-U);
			}
		}
	}
	
	cplex.end();
	model.end();
	u.end();
	w.end();
	x.end();
	z.end();
	return U;
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

double Heuristic(const Knapsack& knap, const Model& mod, vector<double>& xsols)
{
	// create entireprobmodels
	double lb=0.0;
	double ub=knap.totalcost;
	if (mod.option == 4) {
		lb = knap.totalcost;
		ub = 0.0;
	}
	double value=lb/2.0+ub/2.0;

	vector<double> tempxsols(knap.nbItems, 0);
	entireprob entireprobmodels;
	IloConstraint extra_constraint;

	IloEnv entireprobenv;
	entireprobmodels.env = entireprobenv;
	entireprobmodels.x=IloNumVarArray(entireprobenv, knap.nbItems, 0, 1);
	entireprobmodels.pi=IloNumVarArray(entireprobenv, knap.nbScens, 0, IloInfinity);

	entireprobmodels.model = IloModel(entireprobenv);
	IloExpr subobj(entireprobenv),subobjnew(entireprobenv);

	for (int j = 0; j < knap.nbItems; ++j)
		subobj += entireprobmodels.x[j]*knap.cost[j];

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
	entireprobmodels.cplex.setParam(IloCplex::Threads, 1);

	// Main loop
	double LB = 0;
	double UB = 1e8;
	cout << "Heuristic progress" << endl;
	while ((UB-LB)*1.0/UB > 1e-4)
	{
		int num_scenario;
		ub=min(UB,ub);
		lb=min(LB,lb);
		int iteration=0;
		while ((ub-lb)*1.0/ub > 1e-4)
		{
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

						tempxsols[j]=xvals[j];
						obj+=knap.cost[j]*xvals[j];
				}

				xvals.end();

				for (int k = 0; k < knap.nbScens; ++k)
				{			

					double pi1=entireprobmodels.cplex.getValue(entireprobmodels.pi[k]);				

					if (pi1 > 1e-4)
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
			}
			else
			{
				LB=obj;
				lb=max(value,LB);
				value=lb/2.0+ub/2.0;
				for (int j = 0; j < knap.nbItems; ++j)
					xsols[j]=tempxsols[j];
			}
			cout << "current ub = " << ub << ", currenb lb = " << lb << ", num_scenario = " << num_scenario << endl;

			extra_constraint=(subobj>=value);
			entireprobmodels.model.add(extra_constraint);

			iteration++;
		}

		UB=min(UB,ub);
	}
	subobj.end();
	return LB;
}

void V2IP(IloEnv& env, IloModel& model, const Knapsack& knap, Model& mod, double UB)
{
	// UB could be initialized as the quantile bound
	double U = UB;
	// Nonanticipativity
	for (int j = 0; j < knap.nbItems; ++j)
	{
		for (int i = 0; i < knap.nbScens; ++i)
			model.add(mod.u[i*knap.nbItems+j]+mod.w[i*knap.nbItems+j]-mod.x[j] == 0);
	}
	// objbound
	for (int i = 0; i < knap.nbScens; ++i)
	{
		if (mod.fixscen[i] == 0)
		{
			IloExpr lhs(env);
			for (int j = 0; j < knap.nbItems; ++j)
				lhs += mod.u[i*knap.nbItems+j]*knap.cost[j];
			lhs -= mod.z[i]*U;
			lhs -= mod.y;
			model.add(lhs >= -U);
			lhs.end();
			IloExpr lhs2(env);
			for (int j = 0; j < knap.nbItems; ++j)
				lhs2 += mod.w[i*knap.nbItems+j]*knap.cost[j];
			lhs2 += mod.z[i]*U;
			lhs2 -= mod.y;
			model.add(lhs2 >= 0);
			lhs2.end();
		}
		else
		{
			for (int j = 0; j < knap.nbItems; ++j)
				model.add(mod.w[j] <= 0);
		}
	}
	// Au, Du, Dw
	for (int i = 0; i < knap.nbScens; ++i)
	{
		for (int k = 0; k < knap.nbDims; ++k)
		{
			IloExpr lhs(env);
			for (int j = 0; j < knap.nbItems; ++j)
				lhs += mod.u[i*knap.nbItems+j]*knap.weight[i][k*knap.nbItems+j];
			lhs -= knap.capacity[k]*mod.z[i];
			model.add(lhs <= 0);
			lhs.end();
		}
		for (int j = 0; j < knap.nbItems; ++j)
		{
			model.add(mod.u[i*knap.nbItems+j]-mod.z[i] <= 0);
			model.add(mod.w[i*knap.nbItems+j]+mod.z[i] <= 1);
		}
	}
}

void buildSub(const Knapsack& knap, Model& mod)
{
	for (int i = 0; i < knap.nbScens; ++i)
	{
		if (mod.fixscen[i] == 0)
		{
			for (int j = 0; j < knap.nbItems; ++j)
			{
				// u constraints
				IloExpr u(mod.subprobs[i].env);
				u += mod.subprobs[i].alpha[j];
				u += mod.subprobs[i].lambda[0]*knap.cost[j];
				u -= mod.subprobs[i].beta1[j];
				for (int d = 0; d < knap.nbDims; ++d)
					u -= mod.subprobs[i].gamma1[d]*knap.weight[i][d*knap.nbItems+j];
				mod.subprobs[i].model.add(u <= 0);
				u.end();
				
				// w constraints
				IloExpr w(mod.subprobs[i].env);
				w += mod.subprobs[i].alpha[j];
				w += mod.subprobs[i].lambda[1]*knap.cost[j];
				w -= mod.subprobs[i].beta2[j];
				if (mod.subproboption == 2) {
					for (int d = 0; d < knap.nbDims; ++d)
						w -= mod.subprobs[i].gamma2[d]*knap.weight[i][d*knap.nbItems+j];
				}
				else{
					for (int d = 0; d < knap.nbDims; ++d)
						mod.subprobs[i].model.add(mod.subprobs[i].gamma2[d] <= 0);
				}
				mod.subprobs[i].model.add(w <= 0);
				w.end();
			}

			// rho constraints
			IloExpr rho(mod.subprobs[i].env);
			rho += mod.subprobs[i].lambda[0];
			rho += mod.subprobs[i].lambda[1];
			for (int d = 0; d < knap.nbDims; ++d)
				rho += mod.subprobs[i].gamma1[d];
			if (mod.subproboption == 2) {
				for (int d = 0; d < knap.nbDims; ++d)
					rho += mod.subprobs[i].gamma2[d];
			}
			for (int j = 0; j < knap.nbItems; ++j)
			{
				rho += mod.subprobs[i].beta1[j];
				rho += mod.subprobs[i].beta2[j];
			}
			mod.subprobs[i].model.add(rho <= 1);
			rho.end();

			/*
			// rho1,rho2 constraints
			for (int j = 0; j < knap.nbItems; ++j)
			{
				mod.subprobs[i].model.add(mod.subprobs[i].alpha[j]-mod.subprobs[i].xi[j] <= 0);
				mod.subprobs[i].model.add(-mod.subprobs[i].alpha[j]-mod.subprobs[i].xi[j] <= 0);
			}
			*/
			mod.subprobs[i].cplex.setParam(IloCplex::SimDisplay, 0);
			mod.subprobs[i].cplex.setParam(IloCplex::Threads, 1);
		}
	}
}

void genCuts(IloEnv& env, const Knapsack& knap, Model& mod, const IloNumArray& xvals, const IloNumArray& zvals, double yval, double U, vector<BendersCoef>& bcoef, int option, double cutviol)
{
	// option = 0: do LP separation for all scenarios
	// option = 1: do LP separation just for scenarios with zvals that is big enough
	
	for (int i = 0; i < knap.nbScens; ++i)
	{
		if (mod.fixscen[i] == 0 && (option == 0 || zvals[i] > 1-mod.zthreshold))
		{
			// First check if \hat{x} is feasible or not and only proceed if not.
			bool feasflag = 1;
			if (mod.solveoption == 2)
			{
				double cx = 0;
				for (int j = 0; j < knap.nbItems; ++j)
					cx += xvals[j]*knap.cost[j];
				if (yval > cx*zvals[i]+U*(1-zvals[i])+1e-5)
					feasflag = 0;
				if (yval > cx*(1-zvals[i])+U*zvals[i])
					feasflag = 0;
			}
			if (zvals[i] > 1e-5)
			{
				for (int d = 0; d < knap.nbDims; ++d)
				{
					double lhsval = 0;
					for (int j = 0; j < knap.nbItems; ++j)
						lhsval += xvals[j]*knap.weight[i][d*knap.nbItems+j];
					if (lhsval > knap.capacity[d]+1e-5)
					{
						feasflag = 0;
						break;
					}
				}
			}
			else
			{
				for (int d = 0; d < knap.nbDims; ++d)
				{
					double lhsval = 0;
					for (int j = 0; j < knap.nbItems; ++j)
						lhsval += xvals[j]*knap.weight[i][d*knap.nbItems+j];
					if (lhsval > knap.capacity[d]+mod.bigM[d][i]+1e-5)
					{
						feasflag = 0;
						break;
					}
				}
			}
			if (feasflag == 1)
				continue;
			// update the objective coefficients
			if (mod.solveoption == 2)
			{
				mod.subprobs[i].obj.setLinearCoef(mod.subprobs[i].lambda[0], U*zvals[i]+yval-U);
				mod.subprobs[i].obj.setLinearCoef(mod.subprobs[i].lambda[1], -U*zvals[i]+yval);
			}
			for (int d = 0; d < knap.nbDims; ++d)
				mod.subprobs[i].obj.setLinearCoef(mod.subprobs[i].gamma1[d], -zvals[i]*knap.capacity[d]);
			if (mod.subproboption == 2) {
				for (int d = 0; d < knap.nbDims; ++d)
					mod.subprobs[i].obj.setLinearCoef(mod.subprobs[i].gamma2[d], (zvals[i]-1)*(knap.capacity[d]+mod.bigM[d][i]));
			}
			for (int j = 0; j < knap.nbItems; ++j)
			{
				mod.subprobs[i].obj.setLinearCoef(mod.subprobs[i].alpha[j], xvals[j]);
				mod.subprobs[i].obj.setLinearCoef(mod.subprobs[i].beta1[j], -zvals[i]);
				mod.subprobs[i].obj.setLinearCoef(mod.subprobs[i].beta2[j], zvals[i]-1);
			}
			// then solve
			mod.subprobs[i].cplex.solve();
			if (mod.subprobs[i].cplex.getStatus() == IloAlgorithm::Feasible || mod.subprobs[i].cplex.getStatus() == IloAlgorithm::Optimal)
			{
				double subobj = mod.subprobs[i].cplex.getObjValue();
				if (subobj > cutviol)
				{
					// add a cut!
					IloNumArray alphavals(mod.subprobs[i].env);
					IloNumArray beta1vals(mod.subprobs[i].env);
					IloNumArray beta2vals(mod.subprobs[i].env);
					IloNumArray lambdavals(mod.subprobs[i].env);
					IloNumArray gamma1vals(mod.subprobs[i].env);
					IloNumArray gamma2vals(mod.subprobs[i].env);
					mod.subprobs[i].cplex.getValues(alphavals, mod.subprobs[i].alpha);
					mod.subprobs[i].cplex.getValues(beta1vals, mod.subprobs[i].beta1);
					mod.subprobs[i].cplex.getValues(beta2vals, mod.subprobs[i].beta2);
					mod.subprobs[i].cplex.getValues(lambdavals, mod.subprobs[i].lambda);
					mod.subprobs[i].cplex.getValues(gamma1vals, mod.subprobs[i].gamma1);
					mod.subprobs[i].cplex.getValues(gamma2vals, mod.subprobs[i].gamma2);
					double rhsval = 0;
					double sumbeta1 = 0;
					double sumbeta2 = 0;
					for (int j = 0; j < knap.nbItems; ++j)
					{
						sumbeta1 += beta1vals[j];
						sumbeta2 += beta2vals[j];
					}
					double sumd = 0;
					double sumd2 = 0;
					for (int d = 0; d < knap.nbDims; ++d) {
						sumd += knap.capacity[d]*gamma1vals[d];
						sumd2 += (knap.capacity[d]+mod.bigM[d][i])*gamma2vals[d];
					}
					sumd = sumd - sumd2;
					BendersCoef temp;
					temp.xcoef = vector<double>(knap.nbItems);
					temp.zind = i;
					temp.zUcoef = lambdavals[0]-lambdavals[1];
					temp.zrest = sumbeta2-sumbeta1-sumd;
					temp.rhsUcoef = lambdavals[0];
					temp.rhsrest = sumbeta2+sumd2;
					temp.ycoef = lambdavals[0]+lambdavals[1];
					for (int j = 0; j < knap.nbItems; ++j)
						temp.xcoef[j] = alphavals[j];
					double sumalpha = temp.ycoef*temp.ycoef;
					for (int j = 0; j < knap.nbItems; ++j)
						sumalpha += alphavals[j]*alphavals[j];
					temp.normviol = subobj*1.0/sqrt(sumalpha);
					bcoef.push_back(temp);
					alphavals.end();
					beta1vals.end();
					beta2vals.end();
					lambdavals.end();
					gamma1vals.end();
					gamma2vals.end();
				}
			}
			else
			{
				cout << mod.subprobs[i].cplex.getStatus() << endl;
				cout << "subproblem " << i << " is infeasible!" << endl;
				cout << "xvals = ";
				for (int j = 0; j < knap.nbItems; ++j)
					cout << xvals[j] << " ";
    			cout << endl;	
				exit(0);
			}
		}
	}
}

bool checkFeas(int i, const vector<double>& xsol, const Knapsack& knap)
{
	bool returnflag = 1;
	for (int d = 0; d < knap.nbDims; ++d)
	{
		double lhsval = 0;
		for (int j = 0; j < knap.nbItems; ++j)
			lhsval += xsol[j]*knap.weight[i][d*knap.nbItems+j];
		if (lhsval > knap.capacity[d]+1e-5)
		{
			returnflag = 0;
			break;
		}
	}
	return returnflag;
}

void MIPstart(IloNumVarArray& startVar, IloNumArray& startVal, const Knapsack& knap, const Model& mod, const vector<double>& xsol, const vector<bool>& zvals, double heurobj, int option)
{
	// option = 1: +u,w,y; option = 2: +y
	// add x
	for (int j = 0; j < knap.nbItems; ++j)
	{
		startVar.add(mod.x[j]);
		startVal.add(xsol[j]);
	}
	for (int i = 0; i < knap.nbScens; ++i)
	{
		// add z
		startVar.add(mod.z[i]);
		if (zvals[i] == 1)
			startVal.add(1);
		else
			startVal.add(0);
		if (option == 1)
		{
			// add u
			for (int j = 0; j < knap.nbItems; ++j)
				startVar.add(mod.u[i*knap.nbItems+j]);
			// add w
			for (int j = 0; j < knap.nbItems; ++j)
				startVar.add(mod.w[i*knap.nbItems+j]);

			if (zvals[i] == 1)
			{
				for (int j = 0; j < knap.nbItems; ++j)
					startVal.add(xsol[j]);
				for (int j = 0; j < knap.nbItems; ++j)
					startVal.add(0);
			}
			else
			{
				for (int j = 0; j < knap.nbItems; ++j)
					startVal.add(0);
				for (int j = 0; j < knap.nbItems; ++j)
					startVal.add(xsol[j]);
			}
		}
	}
	if (option == 1 || option == 2)
	{
		startVar.add(mod.y);
		startVal.add(heurobj);
	}
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

	IloExpr obj(env);
	for (int i = 0; i < knap.nbScens; ++i)
	{
		for (int j = 0; j < knap.nbItems; ++j)
			obj += s[i*knap.nbItems+j]+t[i*knap.nbItems+j];
	}
	model.add(IloMinimize(env, obj));
	obj.end();
	IloCplex cplex(model);
	cplex.setParam(IloCplex::Threads, 1);
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

double rootLP(IloEnv& env, IloCplex& cplex, IloModel& model, const Knapsack& knap, Model& mod, double U)
{
	// Update and add Benders cut together here?!
	double UB = U; // initial bound
	IloNumVar y(env, -IloInfinity, IloInfinity);
	IloNumVarArray x(env, knap.nbItems, 0, 1);
	IloNumVarArray z(env, knap.nbScens, 0, 1);
	IloModel lpmodel(env);
	// Save all the Benders cuts generated, in the end, only add active ones, though
	vector<BendersCoef> bcoefall;
	IloRangeArray benders(env);
	// Knapsack
	lpmodel.add(IloSum(z) >= (knap.nbScens-mod.p));
	// Upper bound
	IloExpr obj(env);
	if (mod.solveoption == 2)
	{
		lpmodel.add(y <= UB);
		IloExpr lhs(env);
		lhs += y;
		for (int j = 0; j < knap.nbItems; ++j)
			lhs -= x[j]*knap.cost[j];
		lpmodel.add(lhs <= 0);
		lhs.end();
		obj += y;
	}
	if (mod.solveoption == 3)
	{
		for (int j = 0; j < knap.nbItems; ++j)
			obj += x[j]*knap.cost[j];
	}
	lpmodel.add(IloMaximize(env, obj));
	obj.end();
	IloCplex lpcplex(lpmodel);
	lpcplex.setParam(IloCplex::SimDisplay, 0);
	lpcplex.setParam(IloCplex::Threads, 1);
	lpcplex.setParam(IloCplex::PreInd, CPX_OFF);
	lpcplex.setParam(IloCplex::RootAlg, 2);

	bool iterflag = 1;
	IloNumArray xvals(env);
	IloNumArray zvals(env);
	double yval;
	while (iterflag == 1)
	{
		// Doing Benders
		bool Bendersflag = 1;
		double newUB;
		cout << "Benders root solve iterations" << endl;
		while (Bendersflag == 1)
		{
			lpcplex.solve();
			newUB = lpcplex.getObjValue();
			lpcplex.getValues(xvals, x);
			lpcplex.getValues(zvals, z);
			if (mod.solveoption == 2)
				yval = lpcplex.getValue(y);
			else
				yval = 0;
			vector<BendersCoef> bcoef;
			for (int j = 0; j < knap.nbItems; ++j)
			{
				if (xvals[j] < 1e-5)
					xvals[j] = 0;
				if (xvals[j] > 1-1e-5)
					xvals[j] = 1;
			}
			genCuts(env, knap, mod, xvals, zvals, yval, UB, bcoef, 0, mod.cutviol);
			vector<IndexVal> cutlist;
			for (int k = 0; k < bcoef.size(); ++k)
			{
				IndexVal indexval;
				indexval.ind = k;
				indexval.val = bcoef[k].normviol;
				cutlist.push_back(indexval);
			}
			sort(cutlist.begin(), cutlist.end(), less<IndexVal>());
			vector<BendersCoef> added;
			vector<double> norm_added;
			int nadded = 0;
			for (int l = 0; l < cutlist.size(); ++l)
			{
				int k = cutlist[l].ind;
				BendersCoef curbcoef = bcoef[k];
				if (check_parallel(mod, added, curbcoef, norm_added) == 1)
				{
					bcoefall.push_back(bcoef[k]);
					nadded++;
					IloExpr lhs(env);
					for (int j = 0; j < knap.nbItems; ++j)
						lhs += x[j]*curbcoef.xcoef[j];
					if (mod.solveoption == 2)
						lhs += y*curbcoef.ycoef;
					lhs += z[curbcoef.zind]*(curbcoef.zUcoef*UB+curbcoef.zrest);
					IloRange range(env, -IloInfinity, lhs, (bcoef[k].rhsUcoef*UB+bcoef[k].rhsrest));
					lpmodel.add(range);
					benders.add(range);
					lhs.end();
				}
			}	
			// add Benders cut
			//for (int k = 0; k < bcoef.size(); ++k)
			//{
			//	bcoefall.push_back(bcoef[k]);
			//	IloExpr lhs(env);
			//	for (int j = 0; j < knap.nbItems; ++j)
			//		lhs += x[j]*bcoef[k].xcoef[j];
			//	lhs += y*bcoef[k].ycoef;
			//	lhs += z[bcoef[k].zind]*(bcoef[k].zUcoef*UB+bcoef[k].zrest);
			//	IloRange range(env, -IloInfinity, lhs, (bcoef[k].rhsUcoef*UB+bcoef[k].rhsrest));
			//	double lhsval = 0;
			//	for (int j = 0; j < knap.nbItems; ++j)
			//		lhsval += xvals[j]*bcoef[k].xcoef[j];
			//	lhsval += yval*bcoef[k].ycoef;
			//	lhsval += zvals[bcoef[k].zind]*(bcoef[k].zUcoef*UB+bcoef[k].zrest);
			//	lpmodel.add(range);
			//	benders.add(range);
			//	lhs.end();
			//}
			if (bcoef.size() == 0)
				Bendersflag = 0;

			cout << "current obj = " << newUB << ", current cuts = " << nadded <<endl;
		}
		if (mod.solveoption == 3)
		{
			UB = newUB;
			iterflag = 0;
		}
		else
		{
			if (fabs((newUB-UB)*1.0/UB) <= 1e-5)
				iterflag = 0;
			// update the Benders cut coefficients
			for (int k = 0; k < benders.getSize(); ++k)
			{
				// change coefficient
				benders[k].setLinearCoef(z[bcoefall[k].zind], bcoefall[k].zUcoef*newUB+bcoefall[k].zrest);
				benders[k].setUB(bcoefall[k].rhsUcoef*newUB+bcoefall[k].rhsrest);
			}
			y.setUB(newUB);
			UB = newUB;
		}
	}
	// Check if the stored Benders cuts are binding:
	int initcutadded = 0;
	for (int k = 0; k < bcoefall.size(); ++k)
	{
		double lhsval = 0;
		IloExpr lhs(env);
		for (int j = 0; j < knap.nbItems; ++j)
		{
			lhs += mod.x[j]*bcoefall[k].xcoef[j];
			lhsval += xvals[j]*bcoefall[k].xcoef[j];
		}
		if (mod.solveoption == 2)
		{
			lhs += mod.y*bcoefall[k].ycoef;
			lhsval += yval*bcoefall[k].ycoef;
		}
		lhs += mod.z[bcoefall[k].zind]*(bcoefall[k].zUcoef*UB+bcoefall[k].zrest);
		lhsval += zvals[bcoefall[k].zind]*(bcoefall[k].zUcoef*UB+bcoefall[k].zrest);
		if (fabs(lhsval - bcoefall[k].rhsUcoef*UB-bcoefall[k].rhsrest) < 1e-3)
		{
			initcutadded++;
			model.add(lhs <= bcoefall[k].rhsUcoef*UB+bcoefall[k].rhsrest);
		}
		lhs.end();
	}
	cout << "initcutadded = " << initcutadded << endl;
	cout << "Model size: Rows = " << cplex.getNrows() << ", Cols = " << cplex.getNcols() << endl;
	for (int k = 0; k < benders.getSize(); ++k)
		benders[k].end();
	benders.end();
	xvals.end();
	zvals.end();
	lpcplex.end();
	lpmodel.end();
	x.end();
	z.end();
	return UB;
}

bool check_parallel(const Model& mod, vector<BendersCoef>& added, BendersCoef bcoef, vector<double>& norm_added)
{
	// Return 0 if close-to-parallel cuts already exists, return 1 otherwise
	bool returnflag = 1;
	double normA = 0;
	for (int j = 0; j < bcoef.xcoef.size(); ++j)
		normA += bcoef.xcoef[j]*bcoef.xcoef[j];
	normA += bcoef.ycoef*bcoef.ycoef;
	if (added.size() > 0)
	{
		for (int k = 0; k < added.size(); ++k)
		{
			double para = par(bcoef, added[k], norm_added[k], normA);
			if (para > mod.maxpar)
			{
				returnflag = 0;
				break;
			}
		}
	}
	if (returnflag == 1)
	{
		norm_added.push_back(normA);
		added.push_back(bcoef);
	}
	return returnflag;
}

double par(const BendersCoef& bcoef, const BendersCoef& added, double norm_added, double normA)
{
	double denominator = 0;
	for (int j = 0; j < bcoef.xcoef.size(); ++j)
		denominator += bcoef.xcoef[j]*added.xcoef[j];
	denominator += bcoef.ycoef*added.ycoef;
	double returnval = fabs(denominator)*1.0/sqrt(normA*norm_added);
	return returnval;
}

void single_sort(const Knapsack& knap, const Model& mod, int i, vector< vector<double> >& pool)
{
	// fixed constraint: scenario i, for all dims kk
	// solve for different obj with dim = k from 1 to dim, ii = 1 to nbScens
	for (int k = 0; k < knap.nbDims; ++k)
	{
		for (int ii = 0; ii < knap.nbScens; ++ii)
		{
			// Now calculate the problem with one row LP, which can be done by sorting
			double minval = 1e5;
			//for (int kk = 0; kk < knap.nbDims; ++kk)
			int kk = k;
			if (1)
			{
				vector<IndexVal> cost;
				double returnval = 0.0;
				for (int j = 0; j < knap.nbItems; ++j)
				{
					if (knap.weight[i][kk*knap.nbItems+j] < 1e-5)
						returnval += knap.weight[ii][k*knap.nbItems+j];
					else
					{	
						IndexVal indexval;
						indexval.ind = j;
						indexval.val = knap.weight[ii][k*knap.nbItems+j]*1.0/knap.weight[i][kk*knap.nbItems+j];
						cost.push_back(indexval);
					}
				}
				sort(cost.begin(), cost.end(), less<IndexVal>());
				double currentweight = 0.0;
				int iter = 0;
				while (currentweight < knap.capacity[kk] && iter < cost.size())
				{
					if (currentweight + knap.weight[i][kk*knap.nbItems+cost[iter].ind] <= knap.capacity[kk])
					{
						currentweight += knap.weight[i][kk*knap.nbItems+cost[iter].ind];
						returnval += knap.weight[ii][k*knap.nbItems+cost[iter].ind];
					}
					else
					{
						returnval += knap.weight[ii][k*knap.nbItems+cost[iter].ind]*(knap.capacity[kk]-currentweight)*1.0/knap.weight[i][kk*knap.nbItems+cost[iter].ind];
						break;
					}
					iter++;
				}
				returnval -= knap.capacity[k];
				if (returnval < minval)
					minval = returnval;
			}
			pool[k][ii] = minval;
		}
	}
}


double PPlus(const Knapsack& knap, const Model& mod, int i, int k, vector< vector< vector<double> > >& total_pool)
{
	double returnval = 1e8;
	vector<IndexVal> eta;
	for (int ii = 0; ii < knap.nbScens; ++ii)
	{
		IndexVal indexval;
		indexval.ind = ii;
		indexval.val = total_pool[ii][k][i];
		eta.push_back(indexval);
	}		
	sort(eta.begin(), eta.end(), greater<IndexVal>());
	double subval = eta[mod.p].val;
	if (subval < returnval)
		returnval = subval;
	return returnval;
}



/*
primal subproblem

IloNumVarArray u(env, knap.nbItems, 0, IloInfinity);
IloNumVarArray w(env, knap.nbItems, 0, IloInfinity);
IloNumVar rho(env, 0, IloInfinity);
IloModel primal(env);
for (int j = 0; j < knap.nbItems; ++j)
{
	IloExpr lhsx(env);
	lhsx += u[j];
	lhsx += w[j];
	primal.add(lhsx == xvals[j]);
}
IloExpr lhsu1(env);
for (int j = 0; j < knap.nbItems; ++j)
	lhsu1 += u[j]*knap.cost[j];
lhsu1 += rho;
primal.add(lhsu1 >= U*zvals[i]+yval-U);
lhsu1.end();
IloExpr lhsu2(env);
for (int j = 0; j < knap.nbItems; ++j)
	lhsu2 += w[j]*knap.cost[j];
lhsu2 += rho;
primal.add(lhsu2 >= yval-U*zvals[i]);
lhsu2.end();
for (int d = 0; d < knap.nbDims; ++d)
{
	IloExpr lhsg(env);
	for (int j = 0; j < knap.nbItems; ++j)
		lhsg += u[j]*knap.weight[i][d*knap.nbItems+j];
	lhsg += rho;
	primal.add(lhsg >= -knap.capacity[d]*zvals[i]);
	lhsg.end();
}
for (int j = 0; j < knap.nbItems; ++j)
{
	primal.add(-u[j]+rho >= -zvals[i]);
	primal.add(-w[j]+rho >= zvals[i]-1);
}
IloExpr rhoobj(env);
rhoobj += rho;
primal.add(IloMinimize(env, rhoobj));
IloCplex primalcplex(primal);
primalcplex.solve();
primalcplex.setParam(IloCplex::SimDisplay, 0);
double primalobj = primalcplex.getObjValue();
cout << "primalobj = " << primalobj << endl;
cout << "subobj = " << subobj << ", viol = " << viol << endl;


*/
