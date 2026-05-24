option optca = 0.0;
option optcr = 0.0;

set i/1*2/
    j/1*10/
    k/1*20/;

parameters a(i,j,k), maxviol, c(j);

maxviol = 2;

a(i,j,k) = floor(uniform(0.5,1.5));
a(i,'1',k) = 1;
c(j) = 1;
c('1') = card(j) - 1;

binary variables x(j), z(k);
free variable zobj;

* MIP Model
equations obj, cov(i,k), cardin;
obj.. zobj =e= sum(j, c(j)*x(j));
cov(i,k).. sum(j, a(i,j,k)*x(j)) + z(k) =g= 1;
cardin.. sum(k, z(k)) =l= maxviol;
model mipmodel/obj,cov,cardin/;

* Lagrangian model (nonant constraints: sum(k, xl(j,k)) = card(k)*xl(j,'1'))
set cutno/1*100000/;
parameters l(j);
binary variables xl(j,k);
free variable zobjl;
parameter xsol(cutno,j), cutcount;
xsol(cutno,j) = 0;
cutcount = 0;

equations objl, covl(i,k), nonant(j), intcut(cutno,k);
objl.. zobjl =e= (1/card(k))*sum(k, sum(j, c(j)*xl(j,k))) + sum(j, l(j)*(sum(k, xl(j,k)) - card(k)*xl(j,'1')));
covl(i,k).. sum(j, a(i,j,k)*xl(j,k)) + z(k) =g= 1;
nonant(j).. sum(k, xl(j,k)) =e= card(k)*xl(j,'1');
intcut(cutno,k)$(ord(cutno) le cutcount).. sum(j$(xsol(cutno,j) gt 0), (1-xl(j,k))) + sum(j$(xsol(cutno,j) lt 1), xl(j,k)) =g= 1;
model lagmodel/objl,covl,cardin,intcut/;
model fullmodel/objl,covl,cardin,nonant/;

* decomposed lag model (l(j)=0)
binary variables xdl(j);
free variable zobjdl;
parameters adl(i,j);
equations objdl, covdl(i), intcutdl(cutno);
objdl.. zobjdl =e= sum(j, c(j)*xdl(j));
covdl(i).. sum(j, adl(i,j)*xdl(j)) =g= 1;
intcutdl(cutno)$(ord(cutno) le cutcount).. sum(j$(xsol(cutno,j) gt 0), (1-xdl(j))) + sum(j$(xsol(cutno,j) lt 1), xdl(j)) =g= 1;
model dlagmodel/objdl,covdl,intcutdl/;

* knapsack
parameters v(k), u(k), xscen(j,k);
v(k) = 0;
u(k) = 0;
free variable zkobj;
equations kobj, cardin;
kobj.. zkobj =e= (1/card(k))*sum(k, v(k)*(1-z(k)) + u(k)*z(k));
model knapsack/kobj,cardin/;

* Subgradient method
set lagiter/1*1/;
set mainiter/1*50/;
parameter subgrad(j), step, norm, zup, lambda, lb, ub, count, infeas;
parameter mark(k), dis, tol;
alias (k,kk);
l(j) = 0;
zup = c('1'); 
tol = 0.001;

file out/probcov.out/;
put out;
parameter miptime, decomtime;
solve mipmodel minimizing zobj using rmip;
put 'LP bound = ', zobj.l:8:3 /;
solve mipmodel minimizing zobj using mip;
miptime = mipmodel.resusd;
put 'IP bound = ', zobj.l:8:3 /;
lb = 0; 
ub = c('1'); 

cutcount = 0;
decomtime = 0;

*abort$(miptime > 0) "term";

loop(mainiter$(ub-lb gt tol),

	count = 0;
	lambda = 1.5;
	l(j) = 0.0;

$ontext
	loop(lagiter,
		solve lagmodel minmizing zobjl using mip;
		xscen(j,k) = xl.l(j,k);
		subgrad(j) = sum(k, xl.l(j,k)) - card(k)*xl.l(j,'1');
		norm = sqrt(sum(j, subgrad(j)*subgrad(j)));
		lb = max(lb, zobjl.l);
		if (zobjl.l lt lb, 
			count = count + 1;
			else
			count = 0;
		);
		if (count gt 10, lambda = lambda/2; count = 0;);
		step = lambda*(zup - zobj.l)/norm; 
		l(j) = l(j) + step*subgrad(j);
*		put '>> ', ord(lagiter):3:0, ' ', lb:8:3,' ', zobjl.l:8:3, ' ', step:8:5, ' ', norm:8:5, ' ', lambda:8:5, ' ', count:5:0/;
	);
$offtext

*	solve scenario subprobelms
	loop(k,
		adl(i,j) = a(i,j,k);
		solve dlagmodel minimizing zobjdl using mip;
		decomtime = decomtime + dlagmodel.resusd;
		v(k) = zobjdl.l;
		xscen(j,k) = xdl.l(j);
	);

* 	solve knapsack
	u(k) = 0;
	solve knapsack minimizing zkobj using mip;	
	lb = max(lb, zkobj.l); 
	

* 	check each solution mark uniques ones
	mark(k) = 0;
	mark('1') = 1;
	loop(k,
		loop(kk$(ord(kk) gt ord(k)),
			dis = sum(j, abs(xscen(j,k) - xscen(j,kk)));
			if (dis gt 0, mark(kk) = 1; else mark(kk) = 0;); 
		);
	);
	

*	cut off unique solutions
	loop(k,
		if (mark(k) ge 1,
			cutcount = cutcount + 1;
			xsol(cutno,j)$(ord(cutno) eq cutcount) = xscen(j,k);
		);
	);		

* 	check feasibility of unique solutions
	loop(k$(mark(k) gt 0),
		infeas = 0;
		loop(kk,
			if (smin(i, (sum(j, a(i,j,kk)*xscen(j,k)))) lt 1,
				infeas = infeas+1;
			);	
		);
		if (infeas le maxviol,
				ub = min(ub, sum(j, c(j)*xscen(j,k)));
		);
	);

	put  ord(mainiter):3:0, ' ', lb:8:3,' ', ub:8:3, ' ', cutcount:5:0/;

);

put "miptime = ", miptime:8:4 /;
put "dectime = ", decomtime:8:4;
		

