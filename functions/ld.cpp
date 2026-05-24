/*- mode: C++;
  07/02/14
 */

#include <iostream>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <ctime>
#include <algorithm>

using namespace std;

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

double lr2(IloEnv& env, const vector<double>& coef1, const vector<double>& coef2, int nbScens, const vector< vector<double> >& lam, IloNumArray& xvals, int p)
{
	// LR served for V2, solved by an MIP, which could also be solved by bisection
	double returnval;
	IloNumVarArray x(env, 2*nbScens+1, 0, 1);
	IloIntVarArray z(env, nbScens, 0, 1);
	x[2*nbScens].setLB(-IloInfinity);
	x[2*nbScens].setUB(IloInfinity);
	IloModel model(env);
	IloExpr obj(env);
	for (int i = 0; i < nbScens; ++i)
	{
		for (int j = 0; j < 2; ++j)
			obj += lam[i][j]*x[i*2+j];
	}
	obj += x[2*nbScens];
	model.add(IloMaximize(env, obj));
	obj.end();
	int last = x.getSize()-1;
	for (int i = 0; i < nbScens; ++i)
	{
		model.add(2*x[i*2]+3*x[i*2+1]-x[last] >= 0);
		model.add(coef1[i]*x[2*i]+coef2[i]*x[2*i+1]+9*z[i] <= 10);
	}
	model.add(IloSum(z) >= nbScens-p);
	IloCplex cplex(model);
	cplex.setParam(IloCplex::MIPDisplay, 0);
	cplex.solve();
	if (cplex.getStatus() == IloAlgorithm::Optimal)
	{
		cplex.getValues(xvals, x);
		returnval = cplex.getObjValue();
	}
	else
		returnval = -1;
	cplex.end();
	model.end();
	x.end();
	z.end();
	return returnval;
}

double lr1(IloEnv& env, const vector<double>& coef1, const vector<double>& coef2, int nbScens, const vector< vector<double> >& lam, IloNumArray& xvals, int p)
{
	// LR served for V1
	double returnval = 0;
	vector<double> v(nbScens);
	vector<double> u(nbScens);
	xvals = IloNumArray(env, nbScens*2);
	for (int i = 0; i < nbScens; ++i)
	{
		IloNumVarArray x(env, 2, 0, 1);
		IloModel model(env);
		IloExpr obj(env);
		obj += (lam[i][0]+2*1.0/nbScens)*x[0];
		obj += (lam[i][1]+3*1.0/nbScens)*x[1];
		model.add(IloMaximize(env, obj));
		obj.end();
		model.add(coef1[i]*x[0]+coef2[i]*x[1] <= 1);
		IloCplex cplex(model);
		cplex.setParam(IloCplex::SimDisplay, 0);
		cplex.solve();
		if (cplex.getStatus() == IloAlgorithm::Optimal)
		{
			IloNumArray subxval(env);
			cplex.getValues(subxval, x);
			xvals[i*2] = subxval[0];
			xvals[i*2+1] = subxval[1];
			v[i] = cplex.getObjValue();
			subxval.end();
		}
		else
		{
			returnval = -1;
			exit(0);
		}
		u[i] = 0;
		if (lam[i][0]+2*1.0/nbScens > 0)
			u[i] += lam[i][0]+2*1.0/nbScens;
		if (lam[i][1]+3*1.0/nbScens > 0)
			u[i] += lam[i][1]+3*1.0/nbScens;
		x.end();
		model.end();
		cplex.end();
	}
	vector<IndexVal> sortsubobj;
	for (int i = 0; i < nbScens; ++i)
	{
		IndexVal indexval;
		indexval.ind = i;
		indexval.val = v[i]-u[i];
		sortsubobj.push_back(indexval);
	}
	sort(sortsubobj.begin(), sortsubobj.end(), less<IndexVal>());
	for (int i = 0; i < nbScens; ++i)
		returnval += u[i];
	for (int i = 0; i < nbScens-p; ++i)
		returnval += sortsubobj[i].val;
	for (int i = nbScens-p; i < nbScens; ++i)
	{
		int ind = sortsubobj[i].ind;
		if (lam[ind][0]+2*1.0/nbScens > 0)
			xvals[ind*2] = 1;
		else
			xvals[ind*2] = 0;
		if (lam[ind][1]+3*1.0/nbScens > 0)
			xvals[ind*2+1] = 1;
		else
			xvals[ind*2+1] = 0;
	}
	return returnval;
}


double V2(IloEnv& env, vector<double>& coef1, vector<double>& coef2, int nbScens, int p)
{
	// Stronger LD for nonanticipativity
	vector< vector<double> > lam(nbScens);
	for (int i = 0; i < nbScens; ++i)
		lam[i] = vector<double>(2, 0);
	double theta = -1e5;
	bool flag = 1;
	IloNumVarArray lambda(env, nbScens*2+1, -200, 200);
	lambda[nbScens*2].setUB(IloInfinity);
	lambda[nbScens*2].setLB(-IloInfinity);
	IloModel ldmodel(env);
	IloExpr obj(env);
	obj += lambda[nbScens*2];
	ldmodel.add(IloMinimize(env, obj));
	obj.end();
	for (int j = 0; j < 2; ++j)
	{
		IloExpr lhs(env);
		for (int i = 0; i < nbScens; ++i)
			lhs += lambda[i*2+j];
		ldmodel.add(lhs == 0);
		lhs.end();
	}
	IloCplex ldcplex(ldmodel);
	ldcplex.setParam(IloCplex::SimDisplay, 0);
	IloNumArray xvals(env);
	while (flag == 1)
	{
		double result = lr2(env, coef1, coef2, nbScens, lam, xvals, p);
		if (theta < result*(1-1e-5))
		{
			// add a cut
			IloExpr lhs(env);
			lhs += lambda[2*nbScens];
			for (int i = 0; i < nbScens; ++i)
			{
				for (int j = 0; j < 2; ++j)
					lhs -= xvals[i*2+j]*lambda[i*2+j];
			}
			double subtract = 0.0;
			for (int i = 0; i < nbScens; ++i)
			{
				for (int j = 0; j < 2; ++j)
					subtract += xvals[i*2+j]*lam[i][j];
			}
			ldmodel.add(lhs >= result-subtract);
			lhs.end();
			// update lam and theta
			ldcplex.solve();
			if (ldcplex.getStatus() == IloAlgorithm::Optimal)
			{
				IloNumArray temp(env);
				ldcplex.getValues(temp, lambda);
				theta = temp[2*nbScens];
				for (int i = 0; i < nbScens; ++i)
				{
					for (int j = 0; j < 2; ++j)
						lam[i][j] = temp[i*2+j];
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

double V1(IloEnv& env, vector<double>& coef1, vector<double>& coef2, int nbScens, int p)
{
	// LD for relaxing nonanticipativity constraints
	vector< vector<double> > lam(nbScens);
	for (int i = 0; i < nbScens; ++i)
		lam[i] = vector<double>(2, 0);
	double theta = -1e5;
	bool flag = 1;
	IloNumVarArray lambda(env, nbScens*2+1, -50, 50);
	lambda[nbScens*2].setUB(IloInfinity);
	lambda[nbScens*2].setLB(-IloInfinity);
	IloModel ldmodel(env);
	IloExpr ldobj(env);
	ldobj += lambda[nbScens*2];
	ldmodel.add(IloMinimize(env, ldobj));
	ldobj.end();
	for (int j = 0; j < 2; ++j)
	{
		IloExpr lhs(env);
		for (int i = 0; i < nbScens; ++i)
			lhs += lambda[i*2+j];
		ldmodel.add(lhs == 0);
		lhs.end();
	}
	IloCplex ldcplex(ldmodel);
	ldcplex.setParam(IloCplex::SimDisplay, 0);
	IloNumArray xvals(env);
	while (flag == 1)
	{
		double result = lr1(env, coef1, coef2, nbScens, lam, xvals, p);
		if (theta < result*(1-1e-5))
		{
			// add a cut
			IloExpr lhs(env);
			lhs += lambda[2*nbScens];
			for (int i = 0; i < nbScens; ++i)
			{
				for (int j = 0; j < 2; ++j)
					lhs -= xvals[i*2+j]*lambda[i*2+j];
			}
			double subtract = 0.0;
			for (int i = 0; i < nbScens; ++i)
			{
				for (int j = 0; j < 2; ++j)
					subtract += xvals[i*2+j]*lam[i][j];
			}
			ldmodel.add(lhs >= result-subtract);
			lhs.end();
			// update lam and theta
			ldcplex.solve();
			if (ldcplex.getStatus() == IloAlgorithm::Optimal)
			{
				IloNumArray temp(env);
				ldcplex.getValues(temp, lambda);
				theta = temp[2*nbScens];
				for (int i = 0; i < nbScens; ++i)
				{
					for (int j = 0; j < 2; ++j)
						lam[i][j] = temp[i*2+j];
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


double exact(IloEnv& env, vector<double>& coef1, vector<double>& coef2, int nbScens, int p)
{
	double exactobj;
	IloNumVarArray x2(env, 2, 0, 1);
	IloIntVarArray z2(env, nbScens, 0, 1);
	IloModel model2(env);
	IloExpr obj2(env);
	obj2 += 2*x2[0];
	obj2 += 3*x2[1];
	model2.add(IloMaximize(env, obj2));
	obj2.end();
	for (int i = 0; i < nbScens; ++i)
		model2.add(coef1[i]*x2[0]+coef2[i]*x2[1]+9*z2[i] <= 10);
	model2.add(IloSum(z2) >= nbScens-p);
	IloCplex cplex2(model2);
	cplex2.setParam(IloCplex::MIPDisplay, 0);
	cplex2.solve();
	if (cplex2.getStatus() == IloAlgorithm::Optimal)
	{
		exactobj = cplex2.getObjValue();
	}
	else
	  	exactobj = -1;
	cplex2.end();
	model2.end();
	x2.end();
	z2.end();
	return exactobj;
}

double lp(IloEnv& env, vector<double>& coef1, vector<double>& coef2, int nbScens, int p)
{
	double lpobj;
	IloNumVarArray x2(env, 2, 0, 1);
	IloNumVarArray z2(env, nbScens, 0, 1);
	IloModel model2(env);
	IloExpr obj2(env);
	obj2 += 2*x2[0];
	obj2 += 3*x2[1];
	model2.add(IloMaximize(env, obj2));
	obj2.end();
	for (int i = 0; i < nbScens; ++i)
		model2.add(coef1[i]*x2[0]+coef2[i]*x2[1]+9*z2[i] <= 10);
	model2.add(IloSum(z2) >= nbScens-p);
	IloCplex cplex2(model2);
	cplex2.setParam(IloCplex::SimDisplay, 0);
	cplex2.solve();
	if (cplex2.getStatus() == IloAlgorithm::Optimal)
	{
		lpobj = cplex2.getObjValue();
	}
	else
	  	lpobj = -1;
	cplex2.end();
	model2.end();
	x2.end();
	z2.end();
	return lpobj;
}


double primallp(IloEnv& env, vector<double>& coef1, vector<double>& coef2, int nbScens, int p)
{
	double primallpobj;
	IloNumVarArray x(env, 2, 0, 1);
	IloNumVarArray u(env, 2*nbScens, 0, IloInfinity);
	IloNumVarArray w(env, 2*nbScens, 0, IloInfinity);
	IloNumVarArray z(env, nbScens, 0, 1);
	IloModel model(env);
	IloExpr obj(env);
	obj += 2*x[0];
	obj += 3*x[1];
	model.add(IloMaximize(env, obj));
	obj.end();
	for (int i = 0; i < nbScens; ++i)
	{
		IloExpr lhs1(env);
		IloExpr lhs2(env);
		lhs1 += u[2*i];
		lhs1 += w[2*i];
		lhs2 += u[2*i+1];
		lhs2 += w[2*i+1];
		model.add(lhs1 == x[0]);
		model.add(lhs2 == x[1]);
		lhs1.end();
		lhs2.end();
	}
	for (int i = 0; i < nbScens; ++i)
	{
		model.add(coef1[i]*u[2*i]+coef2[i]*u[2*i+1]-z[i] <= 0);
		model.add(w[2*i]+z[i] <= 1);
		model.add(w[2*i+1]+z[i] <= 1);
	}
	model.add(IloSum(z) >= nbScens-p);
	IloCplex cplex(model);
	cplex.setParam(IloCplex::SimDisplay, 0);
	cplex.solve();
	if (cplex.getStatus() == IloAlgorithm::Optimal)
	{
		IloNumArray xvals(env);
		cplex.getValues(xvals, x);
		xvals.end();
		primallpobj = cplex.getObjValue();
	}
	else
	{
		cout << "primallp status = " << cplex.getStatus() << endl;
		exit(0);
	}
	cplex.end();
	model.end();
	x.end();
	z.end();
	return primallpobj;
}

double primallp3(IloEnv& env, vector<double>& coef1, vector<double>& coef2, int nbScens, double UB, int p)
{
	// "primal" formulation W3 for V3
	double primallpobj;
	double U = UB;
	bool fflag = 1;
	while (fflag == 1)
	{
		cout << "U = " << U << endl;
		IloNumVarArray x(env, 3, 0, 1);
		IloNumVarArray u(env, 2*nbScens, 0, IloInfinity);
		IloNumVarArray w(env, 2*nbScens, 0, IloInfinity);
		IloNumVarArray z(env, nbScens, 0, 1);

		x[2].setLB(-IloInfinity);
		x[2].setUB(IloInfinity);
		IloModel model(env);
		IloExpr obj(env);
		obj += x[2];
		model.add(IloMaximize(env, obj));
		obj.end();
		for (int i = 0; i < nbScens; ++i)
		{
			IloExpr lhs1(env);
			IloExpr lhs2(env);
			lhs1 += u[2*i];
			lhs1 += w[2*i];
			lhs2 += u[2*i+1];
			lhs2 += w[2*i+1];
			lhs1 -= x[0];
			lhs2 -= x[1];
			model.add(lhs1 == 0);
			model.add(lhs2 == 0);
			lhs1.end();
			lhs2.end();
		}
		for (int i = 0; i < nbScens; ++i)
		{
			model.add(coef1[i]*u[2*i]+coef2[i]*u[2*i+1]-z[i] <= 0);
			model.add(2*u[2*i]+3*u[2*i+1]-U*z[i]-x[2] >= -U);
			model.add(2*w[2*i]+3*w[2+i+1]+U*z[i]-x[2] >= 0);
			model.add(u[2*i] <= z[i]);
			model.add(u[2*i+1] <= z[i]);
			model.add(w[2*i]+z[i] <= 1);
			model.add(w[2*i+1]+z[i] <= 1);
		}
		model.add(IloSum(z) >= nbScens-p);
		IloCplex cplex(model);
		cplex.setParam(IloCplex::SimDisplay, 0);
		cplex.setParam(IloCplex::PreInd, 0);
		cplex.solve();
		if (cplex.getStatus() == IloAlgorithm::Optimal)
		{
			double objval = cplex.getObjValue();
			if (fabs(objval-U) < 1e-5)
				fflag = 0;
			U = objval;
		}
		else
		{
			cout << "primallp status = " << cplex.getStatus() << endl;
			exit(0);
		}
		cplex.end();
		model.end();
		x.end();
		z.end();
	}
	return U;
}

double quantile(IloEnv& env, vector<double>& coef1, vector<double>& coef2, int nbScens, int p)
{
	vector<IndexVal> tosort;
	for (int i = 0; i < nbScens; ++i)
	{
		IndexVal indexval;
		indexval.ind = i;
		if (2*1.0/coef1[i] > 3*1.0/coef2[i])
		{
			// choose item 1
			if (coef1[i] >= 1)
				indexval.val = 2*1.0/coef1[i];
			else
				indexval.val = 2+(1-coef1[i])*3*1.0/coef2[i];
		}
		else
		{
			// choose item 2
			if (coef2[i] >= 1)
				indexval.val = 3*1.0/coef2[i];
			else
				indexval.val = 3+(1-coef2[i])*2*1.0/coef1[i];
		}
		tosort.push_back(indexval);
	}
	sort(tosort.begin(), tosort.end(), greater<IndexVal>());
	return tosort[p].val;
}


bool sep_mixing(int nbScens, const double& yval, const IloNumArray& zvals, const vector<IndexVal>& subobj, vector<int>& cutind, vector<double>& cutcoef, double& cutrhs, double& violval, int p)
{
	// Separation of mixing: an exact separation
	// Only add the most violated mixing inequality
	vector<IndexVal> sortedz(p+1);
	for (int i = 0; i <= p; i++) 
	{
		sortedz[i].ind = i;
		sortedz[i].val = zvals[subobj[i].ind];
	}
	sort(sortedz.begin(),sortedz.end(), less<IndexVal>()); 
	// Build the best inequality by working from largest zval to smallest
	int maxind = p;
	double lhsval = 0.0;
	double rhsval = subobj[0].val;
	vector<int> subcutind;
	vector<double> subcutcoef;
	for (int i = 0; i <= p; i++)
	{
		int curi = sortedz[i].ind;
		if (curi < maxind)
		{
			subcutind.push_back(subobj[curi].ind);
			subcutcoef.push_back(subobj[maxind].val-subobj[curi].val);
			lhsval += (subobj[maxind].val-subobj[curi].val)*zvals[subobj[curi].ind];
			rhsval += (subobj[maxind].val-subobj[curi].val);
			maxind = curi;
		}
		if (maxind == 0)
			break;
	}
	double viol = yval+lhsval;
	if (viol > rhsval+1e-5)
	{
		cutind = subcutind;
		cutcoef = subcutcoef;
		cutrhs = rhsval;
		violval = viol-rhsval;
		return 1;
	}
	else
		return 0;
}

double mixing(IloEnv& env, vector<double>& coef1, vector<double>& coef2, int nbScens, int p)
{
	double lpobj;
	IloNumVarArray x(env, 2, 0, 1);
	IloNumVarArray z(env, nbScens, 0, 1);
	IloModel model(env);
	IloExpr obj(env);
	obj += 2*x[0];
	obj += 3*x[1];
	model.add(IloMaximize(env, obj));
	obj.end();
	for (int i = 0; i < nbScens; ++i)
		model.add(coef1[i]*x[0]+coef2[i]*x[1]+9*z[i] <= 10);
	model.add(IloSum(z) >= nbScens-p);
	IloCplex cplex(model);
	cplex.setParam(IloCplex::SimDisplay, 0);
	bool flag = 1;
	while (flag == 1)
	{
		flag = 0;
		cplex.solve();
		if (cplex.getStatus() == IloAlgorithm::Optimal || cplex.getStatus() == IloAlgorithm::Feasible)
		{
			lpobj = cplex.getObjValue();
		}
		else
		{	
			cout << "infeasible!" << endl;
	  		exit(0);
		}
		IloNumArray xvals(env);
		cplex.getValues(xvals, x);
		IloNumArray zvals(env);
		cplex.getValues(zvals, z);

		for (int i = 0; i < nbScens; ++i)
		{
			// each row is seen as the coefficient alpha once
			vector<int> cutind;
			vector<double> cutcoef;
			double cutrhs;
			double violval;
			vector<IndexVal> subobj;
			for (int ii = 0; ii < nbScens; ++ii)
			{
				IndexVal indexval;
				indexval.ind = ii;
				if (coef1[i]*1.0/coef1[ii] > coef2[i]*1.0/coef2[ii])
				{
					// choose item 1
					if (coef1[ii] >= 1)
						indexval.val = coef1[i]*1.0/coef1[ii];
					else
						indexval.val = coef1[i]+(1-coef1[ii])*coef2[i]*1.0/coef2[ii];
				}
				else
				{
					// choose item 2
					if (coef2[ii] >= 1)
						indexval.val = coef2[i]*1.0/coef2[ii];
					else
						indexval.val = coef2[i]+(1-coef2[ii])*coef1[i]*1.0/coef1[ii];
				}
				subobj.push_back(indexval);
			}
			sort(subobj.begin(), subobj.end(), less<IndexVal>());

			double yval = coef1[i]*xvals[0]+coef2[i]*xvals[1];
			if (sep_mixing(nbScens, yval, zvals, subobj, cutind, cutcoef, cutrhs, violval, p) == 1)
			{
				// add mixing inequality
				flag = 1;
				IloExpr lhs(env);
				lhs += x[0]*coef1[i];
				lhs += x[1]*coef2[i];

				for (int j = 0; j < cutind.size(); ++j)
				{
					if (fabs(cutcoef[j]) > 1e-5) 
						lhs += z[cutind[j]]*cutcoef[j];
				}
				model.add(lhs <= cutrhs);
				lhs.end();
			}
		}
		xvals.end();
		zvals.end();
	}

	cplex.end();
	model.end();
	x.end();
	z.end();
	return lpobj;
}


int main(int argc, char **argv)
{
	srand (2014);
	int nbScens = atoi(argv[1]);
	double eps = atof(argv[2]);
	int p = int(nbScens*eps);
	cout << "nbScens = " << nbScens << ", p = " << p << endl;
	vector<double> coef1(nbScens);
	vector<double> coef2(nbScens);
	for (int i = 0; i < nbScens; ++i)
	{
		double temp = rand() / double(RAND_MAX);
		coef1[i] = 9*temp + 1;
		coef2[i] = 10-coef1[i];
	}
	IloEnv env;

	// V1
	double v1dual = V1(env, coef1, coef2, nbScens, p);
	cout << "v1 lagrangian dual = " << v1dual << endl;

	// V2
	double v2dual = V2(env, coef1, coef2, nbScens, p);
	cout << "v2 lagrangian dual = " << v2dual << endl;
	
	// exact
	double exactobj = exact(env, coef1, coef2, nbScens, p);
	cout << "exact = " << exactobj << endl;	

	// lp
	double lpobj = lp(env, coef1, coef2, nbScens, p);
	cout << "lp relaxation = " << lpobj << endl;	

	// quantile
	double quantileobj = quantile(env, coef1, coef2, nbScens, p);
	cout << "quantileobj = " << quantileobj << endl;

	// mixing
	//double mixobj = mixing(env, coef1, coef2, nbScens);
	//cout << "mixing relaxation = " << mixobj << endl;	

	// primal lp
	double primallpobj = primallp(env, coef1, coef2, nbScens, p);
	cout << "primallpobj = " << primallpobj << endl;

	// primal lp3
	double primallpobj3 = primallp3(env, coef1, coef2, nbScens, quantileobj, p);
	cout << "primallpobj3 = " << primallpobj3 << endl;

	ofstream out(argv[3], ios::app);
	if (out)  {
			//out << setw(12) << exactobj;
			out << setw(12) << lpobj;
			out << setw(12) << quantileobj;
			//out << setw(12) << v1dual;
			//out << setw(12) << v3dual;
			//out << setw(12) << primallpobj3;
			out << endl;
			out.close();
	}
	env.end();
	return 0;
}



