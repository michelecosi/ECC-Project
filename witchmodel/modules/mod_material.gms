*-------------------------------------------------------------------------------
* Material Resources
*-------------------------------------------------------------------------------

$ifthen %phase%=='conf'

####à
$setglobal mat lit

* Define default resource availability levels for baseline
$if %mat_scenario%=='global_pool' $setglobal material_extr L
$if %mat_scenario%=='indipentent' $setglobal material_extr M

*------------------------------------------------------------------------------
$elseif %phase%=='sets'

set conf /
'material_extr'.'%material_extr%'
/;

*-------------------------------------------------------------------------------
$elseif %phase%=='include_data'

$gdxin '%datapath%data_mod_material'

parameter trade_poly_mat(polydeg,*,n);
$loaddc trade_poly_mat

parameter trade_scale_mat(*);
$loaddc trade_scale_mat

parameter cmat(*,*) 'mat price function coefficients';
$loaddc cmat

$gdxin

scalar countallomat;
scalar maxitermat;
scalar banstart /2150/;

*------------------------------------------------------------------------------
$elseif %phase%=='vars'

* No emissions associated with extraction
Q_EMI_OUT.fx(%mat%,t,n) = 0;

** FPRICE 'World fuel prices [T$/TWh]'
** Q_FUEL 'Total amount of Energy Sources consumed [TWh]'
** Q_OUT 'Extracted primary energy supplies [TWh]'
** MCOST_FUEL 'Average cost of energy sources [T$/TWh]'

** wcum 'World cumulated quantities'
** cum0   'Cumulative extraction at time zero [millions TWh]'

* a      'Marginal cost of extraction'
* c      'Coefficient on cumulative extraction'
* exp    'Power of cumulative extraction'
* scl    'Scale up parameter from 2002 to 2005 prices'
* fast   'Fraction of tot res at which price grows fast'
* extra  'Extra costs (e.g. conversion and enrichment for uranium)'
* res0   'Start ult.ly rec. resources [Millions TWh]'

*-------------------------------------------------------------------------------
$elseif %phase%=='before_solve'

***?? this sets the price for the first iteration
*** NB: poly gives the value of the polynomial of coefficients "cmat" in the point "wcum(%mat%,t) - cexs(%mat%,'cum0')"
if(ord(siter) eq 1,
FPRICE.l(%mat%,t)$(not tfix(t)) = max(valuein(2005,FPRICE.l(%mat%,tt)), 
                                      poly((wcum(%mat%,t) - cexs(%mat%,'cum0')), cmat(polydeg,'%material_extr%')));  

** setting eq 1
* Calculate cumulative production (by means of polynomial functions)
** cumulative extraction at the regional level 
cum_prodpp(%mat%,t,n) = max(0, trade_scale_mat('%material_extr%') * 
                               poly(FPRICE.l(%mat%,t)/(twh2ej/1000), trade_poly_mat(polydeg,'%material_extr%',N)) 
                           ); 
 
** setting eq 2
* Ensure cumulative production in 2005 = 0 
cum_prodpp(%mat%,t,n)$(tfirst(t)) = 0; 


** setting eq 3
* This is to avoid negative production (cumulative production declining) 
loop((t,tp1)$(pre(t,tp1)), 
  cum_prodpp(%mat%,tp1,n)=max(cum_prodpp(%mat%,t,n)+1e-5*tlen(t), cum_prodpp(%mat%,tp1,n) ) 
); 
 
** Calculate annual production 
loop((t,tp1)$(pre(t,tp1)), 
  prodpp(%mat%,t,n) = (cum_prodpp(%mat%,tp1,n) - cum_prodpp(%mat%,t,n))/tlen(t) 
); 
 

** ?? annual productions has a different behaviour after 2100, production share remains constant
prodpp(%mat%,t,n)$(year(t) > 2100) = sum(nn,Q_FUEL.l(%mat%,t,nn)) * valuein(2100, (prodpp(%mat%,tt,n)/sum(nn,prodpp(%mat%,tt,nn)))); 

** ?? if mat trade is not available, the price is set in a different manner, 


$ifthen.ssp set nomattrade 
Q_OUT.fx(%mat%,t,n)$(not tfix(t)) = 0; 
FPRICE.l(%mat%,t)$(not tfix(t)) = cexs(%mat%,'scl') * (cexs(%mat%,'a') + 
                                  cexs(%mat%,'c') * (wcum(%mat%,t) / (cexs(%mat%,'fast') * cexs(%mat%,'res0') ) )**cexs(%mat%,'exp') ) +
                                  cexs(%mat%,'extra'); 
$else.ssp 
Q_OUT.fx(%mat%,t,n)$(not tfix(t)) = prodpp(%mat%,t,n); 
$endif.ssp 


MCOST_FUEL.fx(%mat%,t,n)$(not tfix(t)) = FPRICE.l(%mat%,t) + p_mkup(%mat%,t,n); 
 
);

external("mat")=no;
internal(%mat%)=yes;
*-------------------------------------------------------------------------------
$elseif %phase%=='before_conv'

if (ord(siter) eq 1,
*price w/o extraction constraints (needed to calculate maximum production)
price_sup(%mat%,t)$(not tfix(t)) = max(valuein(2005,FPRICE.l(%mat%,tt)), 
                                        poly((wcum(%mat%,t) - cexs(%mat%,'cum0')), cmat(polydeg,'%material_extr%')));
price_sup(%mat%,t)$(tfix(t)) = FPRICE.lo(%mat%,t);

** setting eq 1
* Calculate cumulative production (by means of polynomial functions, from price)
cum_prodpp(%mat%,t,n) = max(0, trade_scale_mat('%material_extr%') *
                               poly(price_sup(%mat%,t)/(twh2ej/1000), trade_poly_mat(polydeg,'%material_extr%',N)) );

** setting eq 2
* Ensure cumulative production in 2005 = 0
cum_prodpp(%mat%,t,n)$(tfirst(t)) = 0;

** setting eq 3
* This is to avoid negative production (cumulative production declining)
loop((t,tp1)$(pre(t,tp1)),
  cum_prodpp(%mat%,tp1,n)=max(cum_prodpp(%mat%,t,n)+1e-5*tlen(t), cum_prodpp(%mat%,tp1,n) );
);


*** ?? after the ban, cumprod_max is constrained by the max_prod allowed (by the ban)
*** # where is max prod?
cumprod_max(%mat%,t,n) = cum_prodpp(%mat%,t,n);
loop((t,tp1)$(pre(t,tp1) and year(t) ge banstart),
  cumprod_max(%mat%,tp1,n) = cumprod_max(%mat%,t,n) + tlen(t) * max_prod(%mat%,t,n) );  
);

*** ?? it's the same, only with a iteration order contraint
* Maximum production allowed by constraints (DOES THIS WORK?)
loop((t,tp1)$(pre(t,tp1) and not ord(siter) eq 1),
  cumprod_max(%mat%,tp1,n) = cum_prodpp(%mat%,t,n) + tlen(t) * max_prod(%mat%,t,n) );  
  
*entry conditions for the while loop
allo(%mat%,t)=0;
countallomat=0;
glob_polyfun(%mat%,t) = wcum(%mat%,t) - cexs(%mat%,'cum0');


*** set max iteration for the while loop
#if (external("mat"),
maxitermat=100;
#else 
#maxitermat=1;
#);


** STEP 2 - loop the allocation mechanism to account for production cuts
while (countallomat le maxitermat,
**2.6
***allocation: parameter that quantifies the amount of non economical resources that have to be extracted 
***due to low cost resources being unavailable for production
* Update price accounting for cuts  
  glob_polyfun(%mat%,t) = glob_polyfun(%mat%,t) +  allo(%mat%,t);
  price_sup(%mat%,t)$(not tfix(t) and glob_polyfun(%mat%,t) le 12) =  max(valuein(2005,FPRICE.l(%mat%,tt)), 
                                                                          poly(glob_polyfun(%mat%,t), cmat(polydeg,'%material_extr%')));
  price_sup(%mat%,t)$(not tfix(t) and glob_polyfun(%mat%,t) gt 12) =  max(valuein(2005,FPRICE.l(%mat%,tt)), 
                                                                          poly(12, cmat(polydeg,'%material_extr%')));

** 2.7 (spell error on the thesis)
** setting eq 1_new
* Calculate real cumulative production including bans as min between max cumulative production possible for ss constraints and cum production with new price
  cum_prodpp(%mat%,t,n) = max(0, min(cumprod_max(%mat%,t,n),
                                     trade_scale_mat('%material_extr%') * poly(price_sup(%mat%,t)/(twh2ej/1000),
                                                                           trade_poly_mat(polydeg,'%material_extr%',n) )
                                     )
                             ); 

** setting eq 2
* Ensure cumulative production in 2005 = 0
  cum_prodpp(%mat%,t,n)$(tfirst(t))= 0;

** setting eq 3
* This is to avoid negative production (cumulative production declining)
  loop((t,tp1)$(pre(t,tp1)),
    cum_prodpp(%mat%,tp1,n)=max(cum_prodpp(%mat%,t,n)+1e-5*tlen(t), cum_prodpp(%mat%,tp1,n) );
  );

** 2.8 :the difference between cumulative demand and supply, computed as regional sum of the just calculated regional cumulative production values
  allo(%mat%,t) = wcum(%mat%,t) - cexs(%mat%,'cum0') - sum(n,cum_prodpp(%mat%,t,n)*1e-6);
  allo(%mat%,t)$( allo(%mat%,t) le 0) = 0;
  countallomat=countallomat+1;
);
*** end of while

** Calculate annual production
loop((t,tp1)$(pre(t,tp1)),
  prodpp(%mat%,t,n) = (cum_prodpp(%mat%,tp1,n) - cum_prodpp(%mat%,t,n))/tlen(t)
);

*** ?? annual productions has a different behaviour after 2100, production share remains constant
prodpp(%mat%,t,n)$(year(t) gt 2100) = sum(nn,Q_FUEL.l(%mat%,t,nn)) * valuein(2100, (prodpp(%mat%,tt,n)/sum(nn,prodpp(%mat%,tt,nn))));


*** ?? no trade case, this could be adapted as trade distruption
$ifthen.ssp set nomattrade
Q_OUT.fx(%mat%,t,n)$(not tfix(t)) = Q_FUEL.l(%mat%,t,n);
FPRICE.l(%mat%,t)$(not tfix(t)) = cexs(%mat%,'scl') *  (cexs(%mat%,'a') +  
                                       cexs(%mat%,'c') * (wcum(%mat%,t) / (cexs(%mat%,'fast') * cexs(%mat%,'res0')))**cexs(%mat%,'exp') ) +
                                       cexs(%mat%,'extra'); 
$else.ssp
Q_OUT.fx(%mat%,t,n)$(not tfix(t)) = prodpp(%mat%,t,n);
FPRICE.l(%mat%,t)$(not tfix(t))  = price_sup(%mat%,t); # + sum(ssiter$(ord(ssiter) lt ord(siter) and ord(ssiter) gt 1), delta_price(run,ssiter,%mat%,t));
$endif.ssp

*** ?? the mkup could be transoformed in order to simulate "dazi"

MCOST_FUEL.fx(%mat%,t,n)$(not tfix(t)) = FPRICE.l(%mat%,t) + p_mkup(%mat%,t,n);

*------------------------------------------------------------------------------
$elseif %phase%=='gdx_items'

trade_poly_mat
trade_scale_mat

$endif
