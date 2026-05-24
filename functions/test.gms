$title Test for ccdual

set s / 1*3 /;

parameter p /1/;

$ontext
option seed = 2014;

parameters coef1(s), coef2(s);

parameter temp;

loop(s,
	  	temp = uniform(0,1);
		coef1(s) = 9*temp+1;
		coef2(s) = 10-coef1(s);
);

parameter coef1(s)  /1 1.14186, 2 9.25395, 3 8.10207, 4 7.49229, 5 9.93895, 6 9.85421, 7 7.73071, 8 3.10186, 9 1.92665, 10 5.24343/; 
parameter coef2(s) /1 8.85814, 2 0.746052, 3 1.89793, 4 2.50771, 5 0.0610547, 6 0.145789, 7 2.26929, 8 6.89814, 9 8.07335, 10 4.75657/; 
$offtext

parameter coef1(s)  /1 0.5, 2 9.5, 3 5/; 
parameter coef2(s) /1 9.5, 2 0.5, 3 5/; 


*--------------------------------------------------------------------------
* IP
*--------------------------------------------------------------------------

variables
  intz(s)
  ipobj
;

positive variables x1,x2;
binary variables intz(s);

equations
  knapsackint
  ipobjdef
  intbigM(s) coef1[i]*x1+coef2[i]*x2+9*z[i] <= 10
;


ipobjdef.. ipobj =e= 2*x1+3*x2;

knapsackint.. sum(s, intz(s)) =g= card(s)-p;

intbigM(s).. coef1(s)*x1+coef2(s)*x2+9*intz(s) =l= 10;

x1.up = 1;
x2.up = 1;

model ipmod /ipobjdef, knapsackint, intbigM/;
solve ipmod using mip maximizing ipobj;

parameter exact;
exact = ipobj.l;



*--------------------------------------------------------------------------
* LP
*--------------------------------------------------------------------------

variables
  z(s)
  x1
  x2 
  lpobj
;

equations
  knapsack
  lpobjdef
  bigM(s) coef1[i]*x1+coef2[i]*x2+9*z[i] <= 10
;


lpobjdef.. lpobj =e= 2*x1+3*x2;

knapsack.. sum(s, z(s)) =g= card(s)-p;

bigM(s).. coef1(s)*x1+coef2(s)*x2+9*z(s) =l= 10;

x1.up = 1;
x2.up = 1;
z.up(s) = 1;

model lpmod /lpobjdef, knapsack, bigM/;
solve lpmod using lp maximizing lpobj;

parameter lpobjval;
lpobjval = lpobj.l;

*--------------------------------------------------------------------------
* First V1
*--------------------------------------------------------------------------
variables
  u1(s)
  u2(s) 	
  w1(s)
  w2(s) 	
  v1obj   
;

positive variables u1,u2,w1,w2;

equations
  nac1(s) lhs1 == x[0]
  nac2(s) lhs2 == x[1] 
  scen(s) coef1[i]*u[2*i]+coef2[i]*u[2*i+1]-z[i] <= 0
  bound1(s)  w[2*i]+z[i] <= 1
  bound2(s)  w[2*i+1]+z[i] <= 1
  bound3(s)  u[2*i] <= z[i]
  bound4(s)  u[2*i+1] <= z[i]
  v1objdef
;


v1objdef.. v1obj =e= 2*x1+3*x2;

bound1(s).. w1(s)+z(s) =l= 1;

bound2(s).. w2(s)+z(s) =l= 1;

bound3(s).. u1(s)-z(s) =l= 0;

bound4(s).. u2(s)-z(s) =l= 0;

scen(s).. coef1(s)*u1(s)+coef2(s)*u2(s)-z(s) =l= 0;

nac1(s).. u1(s)+w1(s)-x1 =e= 0;

nac2(s).. u2(s)+w2(s)-x2 =e= 0;



model mod1 /v1objdef, knapsack, bound1, bound2, bound3, bound4, scen, nac1, nac2/;
solve mod1 using lp maximizing v1obj;

parameter v1objval;
v1objval = v1obj.l;

*--------------------------------------------------------------------------
* Quantile
*--------------------------------------------------------------------------
positive variables xscen1, xscen2;
parameters scenobjval(s);
parameters scencoef1, scencoef2;
free variable scenobj;
equations
	scenobjdef
	scencon
;

scenobjdef.. scenobj =e= 2*xscen1+3*xscen2;
scencon.. scencoef1*xscen1+scencoef2*xscen2 =l= 1;

xscen1.up = 1;
xscen2.up = 1;
model scenmod /scenobjdef, scencon/;

loop(s,
		scencoef1 = coef1(s);
		scencoef2 = coef2(s);
		solve scenmod using lp maximizing scenobj;
		scenobjval(s) = scenobj.l;
);

parameter quantile;

alias(s,ss);

parameter count;
loop(s, 
		count = 0;
		loop(ss$(ord(s) <> ord(ss)),
			if (scenobjval(ss) < scenobjval(s),
				count = count+1;
			);
		);
		if (count = p,
			quantile = scenobjval(s);
		);
);

$ontext
*--------------------------------------------------------------------------
* Then V2iterLP
*--------------------------------------------------------------------------
parameter U;

variable y,v2iterlpobj;

equations
  yubound1(s) cu^i >= yz^i
  ywbound2(s) cw^i >= y(1-z^i) 
  v2iterlpobjdef	
;

v2iterlpobjdef.. v2iterlpobj =e= y;
yubound1(s).. 2*u1(s)+3*u2(s)-U*z(s)-y =g= -U;
ywbound2(s).. 2*w1(s)+3*w2(s)+U*z(s)-y =g= 0;

model mod2 /v2iterlpobjdef, knapsack, bound1, bound2, bound3, bound4, scen, nac1, nac2, yubound1, ywbound2/;

U = quantile;
parameter v2iterlpobjval;
parameter flag /1/;

while (flag = 1,
	solve mod2 using lp maximizing v2iterlpobj;
	v2iterlpobjval = v2iterlpobj.l;
	display 'v2iterlpobjval = ', v2iterlpobjval;
	if (abs(v2iterlpobjval-U) < 1e-5,
		flag = 0;
	);
	U = v2iterlpobjval;
);



*--------------------------------------------------------------------------
* Then V2LP
*--------------------------------------------------------------------------

variable v2lpobj;

equations
  yubound(s) cu^i >= yz^i
  ywbound(s) cw^i >= y(1-z^i) 
  v2lpobjdef	
;

v2lpobjdef.. v2lpobj =e= y;
yubound(s).. 2*u1(s)+3*u2(s)-y*z(s) =g= 0;
ywbound(s).. 2*w1(s)+3*w2(s)-y+y*z(s) =g= 0;

model mod3 /v2lpobjdef, knapsack, bound1, bound2, bound3, bound4, scen, nac1, nac2, yubound, ywbound/;
solve mod3 using nlp maximizing v2lpobj;

*--------------------------------------------------------------------------
* Then V3LD
*--------------------------------------------------------------------------

parameter LB /0.3/;

parameter UB;

UB = quantile;

parameter l;

free variable v3ldobj;

positive variables s1(s), t1(s), s2(s), t2(s);

equations
  uboundfix1(s) cu^i >= lz^i
  wboundfix2(s) cw^i >= l(1-z^i) 
  relaxnac1(s)     u^i_w^i+s^i-t^i = x
  relaxnac2(s)  
  v3ldobjdefl	
;

v3ldobjdefl.. v3ldobj =e= sum(s, s1(s)) + sum(s, t1(s)) + sum(s, s2(s)) + sum(s, t2(s));
uboundfix1(s).. 2*u1(s)+3*u2(s)-l*z(s) =g= 0;
wboundfix2(s).. 2*w1(s)+3*w2(s)+l*z(s) =g= l;
relaxnac1(s).. u1(s)+w1(s)+s1(s)-t1(s)-x1 =e= 0;
relaxnac2(s).. u2(s)+w2(s)+s2(s)-t2(s)-x2 =e= 0;

model mod4 /v3ldobjdefl, knapsack, bound1, bound2, bound3, bound4, scen, relaxnac1, relaxnac2, uboundfix1, wboundfix2/;

parameter result;
while (abs(UB-LB)>1e-5,
	l = (LB+UB)*1.0/2;	
	solve mod4 using lp minimizing v3ldobj;
	result = v3ldobj.l;
	if (abs(result) < 1e-5,
		LB = l;
	else
		UB = l;
	);
);
$offtext
display 'exact=', exact;
display 'quantile=', quantile;
display 'lpobjval=', lpobjval;
display 'v1objval=', v1objval;

$ontext
display 'v2iterlpobjval=', U;
display 'v2lpobjval=',v2lpobj.l;
display 'v3ldobjval=',l;
$offtext
