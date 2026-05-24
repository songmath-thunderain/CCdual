/*- mode: C++;
 * Date: Jan 27, 2015
 * Do not create so many subproblems
 */

#include <iostream>
#include <ilcplex/ilocplex.h>
#include "ccdual-bounds.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <math.h>
using namespace std;

int main(int argc, char **argv)
{
	// Run program like this:
	// ./ccdual-bounds instancefile eps option resultfile 
	// ./ccdual-bounds instancefile 0.1 1 resultfile
    IloEnv env;
	IloEnv env2;
    Knapsack knap;
	Model mod;
	mod.env = env2;
	mod.cutviol = 1e-5;
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
    file >> knap.cost >> knap.capacity >> knap.weight;
   	knap.nbItems = knap.cost.getSize();
	knap.nbDims = knap.capacity.getSize();
	knap.nbScens = knap.weight.getSize();
	mod.p = int(eps*knap.nbScens+1e-5);
	mod.maxpar = 0.9;
	int option = atoi(argv[3]);
    // option = 1: lp1, lp2, heuristic

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

	double quantileval, heurval, lpval1, lpval2, v1LD1, v1LD2, v2LD1, v2LD2, heurtime, lptime1, lptime2, v1LD1time, v1LD2time, v2LD1time, v2LD2time, bigMtime;

	quantileval = quantile(env, knap, mod);
	cout << "quantileval = " << quantileval << endl;	

	for (int k = 0; k < knap.nbDims; ++k)
	{
		vector<double> subM(knap.nbScens, 0);
		mod.bigM.push_back(subM);
	}

	// Strengthened big-M constraint
	vector< vector< vector<double> > > total_pool;
	double bigMstarttime = clock.getTime();
	// Initial big-M constraint
	for (int i = 0; i < knap.nbScens; ++i)
	{	
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
	bigMtime = clock.getTime()-bigMstarttime;

	if (option == 1)
	{
		// heur, lp1, lp2
		vector<double> xsol(knap.nbItems, 0);
		double heurstart = clock.getTime();
		heurval = Heuristic(knap, mod, xsol);
		heurtime = clock.getTime()-heurstart;
		cout << "heurval = " << heurval << endl;
	
		lp(env, knap, mod, clock, lpval1, lptime1);
		cout << "lpval1 = " << lpval1 << endl;
	}

	if (option == 2 || option == 3)
	{
		// v1LD1 and v2LD1
		// create the subproblems
		mod.fixscen = vector<bool>(knap.nbScens, 0);
		for (int i = 0; i < knap.nbScens; ++i)
		{
			bool fflag = 1;
			for (int k = 0; k < knap.nbDims; k++)
			{
				double bigM = mod.bigM[k][i]; 
				if (bigM > 1e-5)
					fflag = 0;
			}
			if (fflag == 1)
				mod.fixscen[i] = 1;
		}

		if (option == 2)
		{
			rootLP(env, knap, mod, quantileval, 1, 1, clock, v1LD1, v1LD1time);
			cout << "v1LD1 = " << v1LD1 << endl;
		}
		else
		{
			// option == 3
			rootLP(env, knap, mod, quantileval, 2, 1, clock, v2LD1, v2LD1time);
			cout << "v2LD1 = " << v2LD1 << endl;
		}
	}

	// Calculate naive bigM values
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

	if (option == 1)
	{
		// heur, lp1, lp2
		lp(env, knap, mod, clock, lpval2, lptime2);
		cout << "lpval2 = " << lpval2 << endl;
	}

	if (option == 4 || option == 5)
	{
		// v1LD2 and v2LD2
		// create the subproblems

		mod.fixscen = vector<bool>(knap.nbScens, 0);
		for (int i = 0; i < knap.nbScens; ++i)
		{
			bool fflag = 1;
			for (int k = 0; k < knap.nbDims; k++)
			{
				double bigM = mod.bigM[k][i]; 
				if (bigM > 1e-5)
					fflag = 0;
			}
			if (fflag == 1)
				mod.fixscen[i] = 1;
		}

		if (option == 4)
		{
			rootLP(env, knap, mod, quantileval, 1, 0, clock, v1LD2, v1LD2time);
			cout << "v1LD2 = " << v1LD2 << endl;
		}
		else
		{
			// option == 5
			rootLP(env, knap, mod, quantileval, 2, 0, clock, v2LD2, v2LD2time);
			cout << "v2LD2 = " << v2LD2 << endl;
		}
	}
	
	ofstream out(argv[4], ios::app);
	if (out)  {
				out << setw(24) << argv[1];
				out << ", ";
				out << setw(8) << argv[2];
				if (option == 1)
				{
					out << ", ";
					out << setw(8) << quantileval;
					out << ", ";
					out << setw(8) << heurval;
					out << ", ";
					out << setw(8) << lpval1;
					out << ", ";
					out << setw(8) << lpval2;
					out << ", ";
					out << setw(8) << heurtime;
					out << ", ";
					out << setw(8) << lptime1;
					out << ", ";
					out << setw(8) << lptime2;
					out << ", ";
					out << setw(8) << bigMtime;
				}
				if (option == 2)
				{
					out << ", ";
					out << setw(8) << v1LD1;
					out << ", ";
					out << setw(8) << v1LD1time;
				}
				if (option == 3)
				{
					out << ", ";
					out << setw(8) << v2LD1;
					out << ", ";
					out << setw(8) << v2LD1time;
				}
				if (option == 4)
				{
					out << ", ";
					out << setw(8) << v1LD2;
					out << ", ";
					out << setw(8) << v1LD2time;
				}
				if (option == 5)
				{
					out << ", ";
					out << setw(8) << v2LD2;
					out << ", ";
					out << setw(8) << v2LD2time;
				}
				out << endl;
				out.close();
	}
	env.end();
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



double V2LP(IloEnv& env, const Knapsack& knap, const Model& mod, double UB)
{
	double U = UB;
	bool flag = 1;
	while (flag == 1)
	{
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
			extra_constraint=(subobj>=value);
			entireprobmodels.model.add(extra_constraint);

			iteration++;
		}

		UB=min(UB,ub);
	}
	subobj.end();
	return LB;
}

void genCuts(IloEnv& env, const Knapsack& knap, Model& mod, const IloNumArray& xvals, const IloNumArray& zvals, double yval, double U, vector<BendersCoef>& bcoef, int option, bool bigMoption, double threashold)
{
	for (int i = 0; i < knap.nbScens; ++i)
	{
		if (mod.fixscen[i] == 0)
		{
			// First check if \hat{x} is feasible or not and only proceed if not.
			bool feasflag = 1;
			if (option == 2)
			{
				double cx = 0;
				for (int j = 0; j < knap.nbItems; ++j)
					cx += xvals[j]*knap.cost[j];
				if (yval > cx*zvals[i]+U*(1-zvals[i])+1e-5)
					feasflag = 0;
				if (yval > cx*(1-zvals[i])+U*zvals[i]+1e-5)
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
				if (bigMoption == 1)
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
			}
			if (feasflag == 1)
				continue;
			// update the objective coefficients
			buildSub(knap, mod, i);
			if (option == 2)
			{
				mod.subprob.obj.setLinearCoef(mod.subprob.lambda[0], U*zvals[i]+yval-U);
				mod.subprob.obj.setLinearCoef(mod.subprob.lambda[1], -U*zvals[i]+yval);
			}
			for (int d = 0; d < knap.nbDims; ++d)
				mod.subprob.obj.setLinearCoef(mod.subprob.gamma1[d], -zvals[i]*knap.capacity[d]);
			if (bigMoption == 1)
			{
				for (int d = 0; d < knap.nbDims; ++d)
					mod.subprob.obj.setLinearCoef(mod.subprob.gamma2[d], (zvals[i]-1)*(knap.capacity[d]+mod.bigM[d][i]));
			}
			for (int j = 0; j < knap.nbItems; ++j)
			{
				mod.subprob.obj.setLinearCoef(mod.subprob.alpha[j], xvals[j]);
				mod.subprob.obj.setLinearCoef(mod.subprob.beta1[j], -zvals[i]);
				mod.subprob.obj.setLinearCoef(mod.subprob.beta2[j], zvals[i]-1);
			}
			// then solve
			mod.subprob.cplex.solve();
			if (mod.subprob.cplex.getStatus() == IloAlgorithm::Feasible || mod.subprob.cplex.getStatus() == IloAlgorithm::Optimal)
			{
				double subobj = mod.subprob.cplex.getObjValue();
				if (subobj > threashold)
				{
					// add a cut!
					IloNumArray alphavals(mod.subprob.env);
					IloNumArray beta1vals(mod.subprob.env);
					IloNumArray beta2vals(mod.subprob.env);
					IloNumArray lambdavals(mod.subprob.env);
					IloNumArray gamma1vals(mod.subprob.env);
					IloNumArray gamma2vals(mod.subprob.env);
					mod.subprob.cplex.getValues(alphavals, mod.subprob.alpha);
					mod.subprob.cplex.getValues(beta1vals, mod.subprob.beta1);
					mod.subprob.cplex.getValues(beta2vals, mod.subprob.beta2);
					mod.subprob.cplex.getValues(lambdavals, mod.subprob.lambda);
					mod.subprob.cplex.getValues(gamma1vals, mod.subprob.gamma1);
					mod.subprob.cplex.getValues(gamma2vals, mod.subprob.gamma2);
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
					temp.viol = subobj;
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
				cout << mod.subprob.cplex.getStatus() << endl;
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

void rootLP(IloEnv& env, const Knapsack& knap, Model& mod, double U, int option, bool bigMoption, IloTimer& timer, double& val, double& time)
{
	// option = 1: v1LD
	// option = 2: v2LD
	// Update and add Benders cut together here
	double UB = U; // initial bound
	IloNumVar y(env, -IloInfinity, UB);
	IloNumVarArray x(env, knap.nbItems, 0, 1);
	IloNumVarArray z(env, knap.nbScens, 0, 1);
	IloModel lpmodel(env);
	vector<BendersCoef> bcoefall;
	IloRangeArray benders(env);
	// Knapsack
	lpmodel.add(IloSum(z) >= (knap.nbScens-mod.p));
	
	IloExpr obj(env);
	if (option == 2)
	{
		// Upper bound
		IloExpr lhs(env);
		lhs += y;
		for (int j = 0; j < knap.nbItems; ++j)
			lhs -= x[j]*knap.cost[j];
		lpmodel.add(lhs <= 0);
		lhs.end();	
		obj += y;
	}
	if (option == 1)
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

	// create the subproblems
	mod.subprob.env = mod.env;
	mod.subprob.alpha = IloNumVarArray(mod.subprob.env, knap.nbItems, -IloInfinity, IloInfinity);
	mod.subprob.lambda = IloNumVarArray(mod.subprob.env, 2, 0, IloInfinity);
	if (option == 1) {
		mod.subprob.lambda[0].setUB(0.0);
		mod.subprob.lambda[1].setUB(0.0);
	}
	mod.subprob.gamma1 = IloNumVarArray(mod.subprob.env, knap.nbDims, 0, IloInfinity);
	mod.subprob.gamma2 = IloNumVarArray(mod.subprob.env, knap.nbDims, 0, IloInfinity);
	mod.subprob.beta1 = IloNumVarArray(mod.subprob.env, knap.nbItems, 0, IloInfinity);
	mod.subprob.beta2 = IloNumVarArray(mod.subprob.env, knap.nbItems, 0, IloInfinity);
	mod.subprob.model = IloModel(mod.subprob.env);
	mod.subprob.cplex = IloCplex(mod.subprob.model);
	mod.subprob.obj = IloMaximize(mod.subprob.env);
	mod.subprob.model.add(mod.subprob.obj);
	mod.subprob.rangelist = IloRangeArray(mod.subprob.env);
	for (int j = 0; j < knap.nbItems; ++j)
	{
		// u constraints
		IloExpr u(mod.subprob.env);
		u += mod.subprob.alpha[j];
		u += mod.subprob.lambda[0]*knap.cost[j];
		u -= mod.subprob.beta1[j];
		for (int d = 0; d < knap.nbDims; ++d)
			u -= mod.subprob.gamma1[d]*knap.weight[0][d*knap.nbItems+j];
		IloRange range(mod.subprob.env, -IloInfinity, u, 0);
		mod.subprob.rangelist.add(range);
		mod.subprob.model.add(range);
		u.end();
		
		// w constraints
		IloExpr w(mod.subprob.env);
		w += mod.subprob.alpha[j];
		w += mod.subprob.lambda[1]*knap.cost[j];
		w -= mod.subprob.beta2[j];
		if (bigMoption == 1) {
			for (int d = 0; d < knap.nbDims; ++d)
				w -= mod.subprob.gamma2[d]*knap.weight[0][d*knap.nbItems+j];
		}
		else{
			for (int d = 0; d < knap.nbDims; ++d)
				mod.subprob.model.add(mod.subprob.gamma2[d] <= 0);
		}
		IloRange range2(mod.subprob.env, -IloInfinity, w, 0);
		mod.subprob.rangelist.add(range2);
		mod.subprob.model.add(range2);
		w.end();
	}

	// rho constraints
	IloExpr rho(mod.subprob.env);
	rho += mod.subprob.lambda[0];
	rho += mod.subprob.lambda[1];
	for (int d = 0; d < knap.nbDims; ++d)
		rho += mod.subprob.gamma1[d];
	if (bigMoption == 1) {
		for (int d = 0; d < knap.nbDims; ++d)
			rho += mod.subprob.gamma2[d];
	}
	for (int j = 0; j < knap.nbItems; ++j)
	{
		rho += mod.subprob.beta1[j];
		rho += mod.subprob.beta2[j];
	}
	mod.subprob.model.add(rho <= 1);
	mod.subprob.cplex.setParam(IloCplex::SimDisplay, 0);
	mod.subprob.cplex.setParam(IloCplex::Threads, 1);


	bool iterflag = 1;
	IloNumArray xvals(env);
	IloNumArray zvals(env);
	double yval;
	double starttime;
	while (iterflag == 1)
	{
		// Doing Benders
		bool Bendersflag = 1;
		double newUB;
		while (Bendersflag == 1)
		{
			lpcplex.solve();
			newUB = lpcplex.getObjValue();
			cout << "newUB = " << newUB << endl;
			lpcplex.getValues(xvals, x);
			lpcplex.getValues(zvals, z);
			if (option == 2)
				yval = lpcplex.getValue(y);
			else
				yval = 0;
			vector<BendersCoef> bcoef;
			genCuts(env, knap, mod, xvals, zvals, yval, UB, bcoef, option, bigMoption, mod.cutviol);
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
			for (int l = 0; l < cutlist.size(); ++l)
			{
				int k = cutlist[l].ind;
				BendersCoef curbcoef = bcoef[k];
				if (check_parallel(mod, added, curbcoef, norm_added) == 1)
				{
					bcoefall.push_back(bcoef[k]);
					IloExpr lhs(env);
					double tempnorm = 0;
					for (int j = 0; j < knap.nbItems; ++j)
					{
						tempnorm += curbcoef.xcoef[j]*curbcoef.xcoef[j];
						lhs += x[j]*curbcoef.xcoef[j];
					}
					if (option == 2)
						lhs += y*curbcoef.ycoef;
					lhs += z[curbcoef.zind]*(curbcoef.zUcoef*UB+curbcoef.zrest);
					IloRange range(env, -IloInfinity, lhs, (curbcoef.rhsUcoef*UB+curbcoef.rhsrest));
					lpmodel.add(range);
					benders.add(range);
					lhs.end();
					added.push_back(curbcoef);
					norm_added.push_back(tempnorm);
				}
			}	

			if (bcoef.size() == 0)
				Bendersflag = 0;
		}
		if (option == 1)
		{
			UB = newUB;
			iterflag = 0;
		}
		else
		{
			if (fabs((newUB-UB)*1.0/UB) <= 1e-5)
				iterflag = 0;
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

	time = timer.getTime()-starttime;
	val = UB;
	
	for (int k = 0; k < benders.getSize(); ++k)
		benders[k].end();
	benders.end();
	xvals.end();
	zvals.end();
	lpcplex.end();
	lpmodel.end();
	x.end();
	z.end();
	mod.subprob.cplex.end();
	mod.subprob.model.end();
	mod.subprob.obj.end();
	mod.subprob.alpha.end();
	mod.subprob.lambda.end();
	mod.subprob.gamma1.end();
	mod.subprob.gamma2.end();
	mod.subprob.beta1.end();
	mod.subprob.beta2.end();
	mod.subprob.env.end();
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


void lp(IloEnv& env, const Knapsack& knap, const Model& mod, IloTimer& timer, double& lpval, double& lptime)
{
	IloNumVarArray x(env, knap.nbItems, 0, 1);
	IloNumVarArray z(env, knap.nbScens, 0, 1);
	IloModel model(env);
	model.add(IloSum(z) >= knap.nbScens-mod.p);
	IloExpr obj(env);
	for (int j = 0; j < knap.nbItems; ++j)
		obj += x[j]*knap.cost[j];
	model.add(IloMaximize(env, obj));
	obj.end();
	IloCplex cplex(model);
	cplex.setParam(IloCplex::SimDisplay, 0);
	cplex.setParam(IloCplex::Threads, 1);
	vector<bool> fixscen(knap.nbScens, 0);
	for (int i = 0; i < knap.nbScens; ++i)
	{
		bool fflag = 1;
    	for (int k = 0; k < knap.nbDims; k++)
    	{
			double bigM = mod.bigM[k][i]; 
			if (bigM > 1e-5)
			{
				fflag = 0;
				IloExpr lhs(env);
				for (int j = 0; j < knap.nbItems; ++j)
					lhs += x[j]*knap.weight[i][k*knap.nbItems+j];
				lhs += z[i]*bigM;
				model.add(lhs <= (knap.capacity[k]+bigM));
				lhs.end();
			}
    	}
		if (fflag == 1)
			fixscen[i] = 1;
	}
	for (int i = 0; i < knap.nbScens; ++i)
	{
		if (fixscen[i] == 1)
			model.add(z[i] >= 1);
	}
	double start = timer.getTime();
	cplex.solve();
	lptime += timer.getTime()-start;
	if (cplex.getStatus() != IloAlgorithm::Feasible && cplex.getStatus() != IloAlgorithm::Optimal)
	{
		cout << "lp ERROR!" << endl;
		exit(0);
	}
	lpval = cplex.getObjValue();
	cplex.end();
	model.end();
	x.end();
	z.end();
}

double V1LP(IloEnv& env, const Knapsack& knap, const Model& mod)
{
	double returnval;
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
	// Au, Du, Dw
	for (int i = 0; i < knap.nbScens; ++i)
	{
		if (mod.fixscen[i] == 0)
		{
			for (int k = 0; k < knap.nbDims; ++k)
			{
				IloExpr lhs(env);
				for (int j = 0; j < knap.nbItems; ++j)
					lhs += u[i*knap.nbItems+j]*knap.weight[i][k*knap.nbItems+j];
				lhs -= knap.capacity[k]*z[i];
				model.add(lhs <= 0);
				lhs.end();
				IloExpr lhs2(env);
				for (int j = 0; j < knap.nbItems; ++j)
					lhs2 += w[i*knap.nbItems+j]*knap.weight[i][k*knap.nbItems+j];
				lhs2 += (knap.capacity[k]+mod.bigM[k][i])*z[i];
				model.add(lhs2 <= knap.capacity[k]+mod.bigM[k][i]);
				lhs2.end();
			}
			for (int j = 0; j < knap.nbItems; ++j)
			{
				model.add(u[i*knap.nbItems+j]-z[i] <= 0);
				model.add(w[i*knap.nbItems+j]+z[i] <= 1);
			}
		}
		else
			model.add(z[i] >= 1);
	}
	// Knapsack
	model.add(IloSum(z) >= (knap.nbScens-mod.p));

	IloExpr obj(env);
	for (int j = 0; j < knap.nbItems; ++j)
		obj += knap.cost[j]*x[j];
	model.add(IloMaximize(env, obj));
	obj.end();
	IloCplex cplex(model);
	cplex.setParam(IloCplex::SimDisplay, 0);
	cplex.setParam(IloCplex::Threads, 1);
	bool flag = 1;
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

double rootLP2(IloEnv& env, const Knapsack& knap, Model& mod, double U, int option, bool bigMoption)
{
	double UB = U;
	IloNumVar y(env, -IloInfinity, UB);
	IloNumVarArray x(env, knap.nbItems, 0, 1);
	IloNumVarArray z(env, knap.nbScens, 0, 1);
	IloNumVarArray u(env, knap.nbItems*knap.nbScens, 0, 1);
	IloNumVarArray w(env, knap.nbItems*knap.nbScens, 0, 1);
	IloModel lpmodel(env);
	lpmodel.add(IloSum(z) >= (knap.nbScens-mod.p));
	IloExpr obj(env);
	IloRangeArray ranges(env);
	if (option == 2)
	{
		IloExpr lhs(env);
		lhs += y;
		for (int j = 0; j < knap.nbItems; ++j)
			lhs -= x[j]*knap.cost[j];
		lpmodel.add(lhs <= 0);
		lhs.end();
		obj += y;
		for (int i = 0; i < knap.nbScens; ++i)
		{
			IloExpr lhs1(env);
			lhs1 += y;
			lhs1 += z[i]*UB;
			for (int j = 0; j < knap.nbItems; ++j)
				lhs1 -= knap.cost[j]*u[i*knap.nbItems+j];
			IloRange range(env, -IloInfinity, lhs1, UB);
			ranges.add(range);
			lpmodel.add(range);
			lhs1.end();

			IloExpr lhs2(env);
			lhs2 += y;
			lhs2 -= z[i]*UB;
			for (int j = 0; j < knap.nbItems; ++j)
				lhs2 -= knap.cost[j]*w[i*knap.nbItems+j];
			IloRange range2(env, -IloInfinity, lhs2, 0);
			ranges.add(range2);
			lpmodel.add(range2);
			lhs2.end();
		}
	}
	if (option == 1)
	{
		for (int j = 0; j < knap.nbItems; ++j)
			obj += x[j]*knap.cost[j];
	}
	lpmodel.add(IloMaximize(env, obj));
	obj.end();
	for (int i = 0; i < knap.nbScens; ++i)
	{
		for (int d = 0; d < knap.nbDims; ++d)
		{
			if (mod.bigM[d][i] >= 0)
			{
				IloExpr lhs(env);
				for (int j = 0; j < knap.nbItems; ++j)
					lhs += u[i*knap.nbItems+j]*knap.weight[i][d*knap.nbItems+j];
				lhs -= knap.capacity[d]*z[i];
				lpmodel.add(lhs <= 0);
				lhs.end();
			}
			else
			{
				IloExpr lhs(env);
				for (int j = 0; j < knap.nbItems; ++j)
					lhs += u[i*knap.nbItems+j]*knap.weight[i][d*knap.nbItems+j];
				lhs -= (knap.capacity[d]+mod.bigM[d][i])*z[i];
				lpmodel.add(lhs <= 0);
				lhs.end();
			}
			IloExpr lhs2(env);
			for (int j = 0; j < knap.nbItems; ++j)
				lhs2 += w[i*knap.nbItems+j]*knap.weight[i][d*knap.nbItems+j];
			lhs2 += (knap.capacity[d]+mod.bigM[d][i])*z[i];
			lpmodel.add(lhs2 <= (knap.capacity[d]+mod.bigM[d][i]));
			lhs2.end();
		}
		for (int j = 0; j < knap.nbItems; ++j)
		{
			lpmodel.add(u[i*knap.nbItems+j]-z[i] <= 0);
			lpmodel.add(w[i*knap.nbItems+j]+z[i] <= 1);
		}
	}
	for (int j = 0; j < knap.nbItems; ++j)
	{
		for (int i = 0; i < knap.nbScens; ++i)
			lpmodel.add(u[i*knap.nbItems+j]+w[i*knap.nbItems+j]-x[j] == 0);
	}
	IloCplex lpcplex(lpmodel);
	lpcplex.setParam(IloCplex::SimDisplay, 0);
	lpcplex.setParam(IloCplex::Threads, 1);
	lpcplex.setParam(IloCplex::PreInd, CPX_OFF);

	bool iterflag = 1;
	IloNumArray xvals(env);
	IloNumArray zvals(env);
	double yval;
	while (iterflag == 1)
	{
		double newUB;
		lpcplex.solve();
		newUB = lpcplex.getObjValue();
		if (option == 1)
		{
			UB = newUB;
			iterflag = 0;
		}
		if (option == 2)
		{
			if (fabs((newUB-UB)*1.0/UB) <= 1e-5)
				iterflag = 0;
			int iter = 0;
			for (int i = 0; i < knap.nbScens; ++i)
			{
				ranges[iter].setUB(newUB);
				ranges[iter].setLinearCoef(z[i], newUB);
				iter++;
				ranges[iter].setLinearCoef(z[i], -newUB);
				iter++;
			}
			y.setUB(newUB);
			UB = newUB;
		}
	}
	/*
	benders.end();
	xvals.end();
	zvals.end();
	lpcplex.end();
	lpmodel.end();
	x.end();
	z.end();
	*/
	return UB;
}

void buildSub(const Knapsack& knap, Model& mod, int i)
{
	int iter = 0;
	for (int j = 0; j < knap.nbItems; ++j)
	{
		for (int d = 0; d < knap.nbDims; ++d)
			mod.subprob.rangelist[iter].setLinearCoef(mod.subprob.gamma1[d], -knap.weight[i][d*knap.nbItems+j]);
		iter++;
		// w constraints
		if (mod.subproboption == 2) {
			for (int d = 0; d < knap.nbDims; ++d)
				mod.subprob.rangelist[iter].setLinearCoef(mod.subprob.gamma2[d], -knap.weight[i][d*knap.nbItems+j]);
		}
		iter++;
	}
}

