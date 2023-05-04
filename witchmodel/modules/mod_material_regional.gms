*-------------------------------------------------------------------------------
* lit Resources 
*-------------------------------------------------------------------------------
$ifthen %phase%=='conf'

* Definition of the global variable specific to the module 
*$setglobal nolittrade
*-------------------------------------------------------------------------------
$elseif %phase%=='sets'

set fuel /lit/;
set f /lit/;
set extract(f) / lit/;

*-------------------------------------------------------------------------------
$elseif %phase%=='compute_data'
parameter cum_0 /27.94/

parameter trade_scale_lit /1/;
*-------------------------------------------------------------------------------
$elseif %phase%=='include_data'

$gdxin '%datapath%data_mod_lithium'

parameter trade_poly_lit(polydeg,n);
$loaddc trade_poly_lit

$gdxin
*-------------------------------------------------------------------------------
$elseif %phase%=='vars'

* No emissions associated with extraction
Q_EMI_OUT.fx('lit',t,n) = 0;
*-------------------------------------------------------------------------------
$elseif %phase%=='compute_data'
Parameters
 deg0 /27.9416/
 deg4 /7.41426e-16/;

parameter 
  lit_price_coeff /deg0,0,0,0,deg4/;


*-------------------------------------------------------------------------------
$elseif %phase%=='before_solve'
                                                 
*FPRICE.l('lit',t)$(not tfix(t)) = max(valuein(2005,FPRICE.l('lit',tt)), deg0*(wcum('lit',t) - cum_0)+deg4*(wcum('lit',t) - cum_0)**4); #poly((wcum('lit',t) - cum_0), lit_price_coeff));

* Calculate cumulative production (by means of polynomial functions)
cum_prodpp('lit',t,n) = max(0, trade_scale_lit*
                               poly(FPRICE.l('lit',t)/(twh2ej/1000), trade_poly_lit(polydeg,n))
                           );

* Ensure cumulative production in 2005 = 0
cum_prodpp('lit',t,n)$(tfirst(t))= 0;


* This is to avoid negative production (cumulative production declining)
loop((t,tp1)$(pre(t,tp1)),
  cum_prodpp('lit',tp1,n)=max(cum_prodpp('lit',t,n)+1e-5*tlen(t), cum_prodpp('lit',tp1,n) )
);

** Calculate annual production
loop((t,tp1)$(pre(t,tp1)),
  prodpp('lit',t,n) = (cum_prodpp('lit',tp1,n) - cum_prodpp('lit',t,n))/tlen(t)
);

prodpp('lit',t,n)$(year(t) > 2100) = sum(nn,Q_FUEL.l('lit',t,nn)) * valuein(2100, (prodpp('lit',tt,n)/sum(nn,prodpp('lit',tt,nn))));

$ifthen.ssp set nolittrade
Q_OUT.fx('lit',t,n)$(not tfix(t)) = 0;
* fprice defined on the global
$else.ssp
Q_OUT.fx('lit',t,n)$(not tfix(t)) = prodpp('lit',t,n);
$endif.ssp

MCOST_FUEL.fx('lit',t,n)$(not tfix(t)) = FPRICE.l('lit',t); # + p_mkup('lit',t,n);

$elseif %phase%=='gdx_items'

trade_poly_lit

$endif
