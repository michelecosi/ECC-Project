$title DICE economy homework
$onmulti
$setenv gdxcompress 1
$onrecurse
$eolcom #
              
$onText
Simple one sector nonlinear optimal growth model.
$offText

Set
   t     'extended horizon' / 1*40 /
   tb(t) 'initial period'
   tt(t) 'terminal period';
   
alias(t,ttt);

tb(t) = yes$(ord(t) = 1);
tt(t) = yes$(ord(t) = card(t));
display tb, tt;

Scalar
   dk   'rate of depreciation'      /   .1   /
   el   'exponent on capital'       /   .3   /
   r    'labor force growth rate'   /   .015 /
   eta  'elasticity'                /   1.45 /
   z    'technical progress'        /   .01  /
   rho  'welfare discount'          /   .015 /
   k0   'initial capital'           / 1.3      /
   pop0 'initial population'        / 0.1    /
   tfp0 'total factor productivity' / 2      /
   tstep 'years in a period'        / 5      /;


Parameter
   rr(t)   'discount factor'
   pop(t)  'population'
   tfp(t) 'total factor productivity'
      
   sigma(t)  'emission intensity'  ;
   
rr(t)   = 1/(1 + rho)**( tstep * ( t.val-1) );
pop(t) = min(pop0 * (1 + r) **( (ord(t) - 1) * tstep), pop0*2)  ;
tfp(t) = min(tfp0 * (1 + z * (ord(t) - 1) * tstep),tfp0*2.5) ;

sigma(t) = yes$(ord(t) = 1);
sigma(t) = max(1-0.01*(ord(t)-1)*tstep,0);

Variable
   S(t)  'saving rates'
   Y(t)  'income'
   K(t)  'capital stock'
   U     'utility'
   
   E(t)  'emissions'
   CUMEMI(t) 'cumulative emissions';


Positive variables S,Y,K;

* constraints
S.up(t) = 1;
Y.lo(t) = 1e-5;

* initial and terminal conditions
K.fx(tb) = k0;
S.lo(t)$(ord(t) ge card(t) - 10)  = 0.2;

* initialization
S.l(t) = 0.25;

Equation
   eq_y 'capital stock balance'
   eq_k 'income definition'
   eq_util    'welfare function'
   eq_emi 'emission equation'
   eq_cumemi 'cumulative emissions';

eq_y(t)..    Y(t)   =e= tfp(t) * K(t)**el * pop(t)**(1-el)-0.0003*CUMEMI(t)**2;

eq_k(t+1)..  K(t+1) =e= Y(t) * S(t) * tstep + (1 - dk)**tstep * K(t);

eq_util..    U      =e= sum(t, tstep * rr(t) * pop(t) * (((Y(t) * ( 1 - S(t) ) / pop(t))**(1 - eta) - 1)/( 1 - eta) - 1 ) );

eq_emi(t)..     E(t)   =e= Y(t)*sigma(t);
eq_cumemi(t)..  CUMEMI(t) =e= sum(ttt$(ttt.val le t.val),E(ttt)*tstep);


Model growth / all /;

solve growth maximizing U using nlp;

parameter g(t);

g(t) = div0(Y.l(t) - Y.l(t-1), Y.l(t-1) );

display S.l


execute_unload "data.gdx"