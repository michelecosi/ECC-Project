*-------------------------------------------------------------------------------
* Gas Resources
*-------------------------------------------------------------------------------

$ifthen %phase%=='conf'

* Define default resource availability levels for baseline
$if %baseline%=='ssp1' $setglobal fossil_gas L
$if %baseline%=='ssp2' $setglobal fossil_gas M
$if %baseline%=='ssp3' $setglobal fossil_gas H
$if %baseline%=='ssp4' $setglobal fossil_gas M
$if %baseline%=='ssp5' $setglobal fossil_gas H
$if %baseline%=='ssp5' $setglobal nogastrade
$if %baseline%=='ssp3' $setglobal nogastrade

*------------------------------------------------------------------------------
$elseif %phase%=='sets'

**** usless set
set extract(f) 'Primary energy fuels with an extraction sector' / gas/;

set conf /
'fossil_gas'.'%fossil_gas%'
/;

*-------------------------------------------------------------------------------
$elseif %phase%=='include_data'

$gdxin '%datapath%data_mod_gas'

parameter trade_poly_gas(polydeg,*,n);
$loaddc trade_poly_gas

parameter trade_scale_gas(*);
$loaddc trade_scale_gas

parameter cgas(*,*) 'Gas price function coefficients';
$loaddc cgas

$gdxin

scalar countallogas;
scalar maxitergas;
scalar banstart /2150/;

*------------------------------------------------------------------------------
$elseif %phase%=='vars'

* No emissions associated with extraction
Q_EMI_OUT.fx('gas',t,n) = 0;

** FPRICE 'World fuel prices [T$/TWh]'
** Q_FUEL 'Total amount of Energy Sources consumed [TWh]'
** Q_OUT 'Extracted primary energy supplies [TWh]'
** MCOST_FUEL 'Average cost of energy sources [T$/TWh]'
** Q_IN 'Imported/consumed PES [TWh]';
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
*** NB: poly gives the value of the polynomial of coefficients "cgas" in the point "wcum('gas',t) - cexs('gas','cum0')"
if(ord(siter) eq 1,
FPRICE.l('gas',t)$(not tfix(t)) = max(valuein(2005,FPRICE.l('gas',tt)), 
                                      poly((wcum('gas',t) - cexs('gas','cum0')), cgas(polydeg,'%fossil_gas%')));  

** setting eq 1
* Calculate cumulative production (by means of polynomial functions)
** cumulative extraction at the regional level 
cum_prodpp('gas',t,n) = max(0, trade_scale_gas('%fossil_gas%') * 
                               poly(FPRICE.l('gas',t)/(twh2ej/1000), trade_poly_gas(polydeg,'%fossil_gas%',N)) 
                           ); 
 
** setting eq 2
* Ensure cumulative production in 2005 = 0 
cum_prodpp('gas',t,n)$(tfirst(t)) = 0; 


** setting eq 3
* This is to avoid negative production (cumulative production declining) 
loop((t,tp1)$(pre(t,tp1)), 
  cum_prodpp('gas',tp1,n)=max(cum_prodpp('gas',t,n)+1e-5*tlen(t), cum_prodpp('gas',tp1,n) ) 
); 
 
** Calculate annual production 
loop((t,tp1)$(pre(t,tp1)), 
  prodpp('gas',t,n) = (cum_prodpp('gas',tp1,n) - cum_prodpp('gas',t,n))/tlen(t) 
); 
 

** ?? annual productions has a different behaviour after 2100, production share remains constant
prodpp('gas',t,n)$(year(t) > 2100) = sum(nn,Q_FUEL.l('gas',t,nn)) * valuein(2100, (prodpp('gas',tt,n)/sum(nn,prodpp('gas',tt,nn)))); 

** ?? if gas trade is not available, the price is set in a different manner, 


$ifthen.ssp set nogastrade 
Q_OUT.fx('gas',t,n)$(not tfix(t)) = 0; 
FPRICE.l('gas',t)$(not tfix(t)) = cexs('gas','scl') * (cexs('gas','a') + 
                                  cexs('gas','c') * (wcum('gas',t) / (cexs('gas','fast') * cexs('gas','res0') ) )**cexs('gas','exp') ) +
                                  cexs('gas','extra'); 
$else.ssp 
Q_OUT.fx('gas',t,n)$(not tfix(t)) = prodpp('gas',t,n); 
$endif.ssp 


MCOST_FUEL.fx('gas',t,n)$(not tfix(t)) = FPRICE.l('gas',t) + p_mkup('gas',t,n); 
 
);

external("gas")=no;
internal('gas')=yes;
*-------------------------------------------------------------------------------
$elseif %phase%=='before_conv'

if (ord(siter) eq 1,
*price w/o extraction constraints (needed to calculate maximum production)
price_sup('gas',t)$(not tfix(t)) = max(valuein(2005,FPRICE.l('gas',tt)), 
                                        poly((wcum('gas',t) - cexs('gas','cum0')), cgas(polydeg,'%fossil_gas%')));
price_sup('gas',t)$(tfix(t)) = FPRICE.lo('gas',t);

** setting eq 1
* Calculate cumulative production (by means of polynomial functions, from price)
cum_prodpp('gas',t,n) = max(0, trade_scale_gas('%fossil_gas%') *
                               poly(price_sup('gas',t)/(twh2ej/1000), trade_poly_gas(polydeg,'%fossil_gas%',N)) );

** setting eq 2
* Ensure cumulative production in 2005 = 0
cum_prodpp('gas',t,n)$(tfirst(t)) = 0;

** setting eq 3
* This is to avoid negative production (cumulative production declining)
loop((t,tp1)$(pre(t,tp1)),
  cum_prodpp('gas',tp1,n)=max(cum_prodpp('gas',t,n)+1e-5*tlen(t), cum_prodpp('gas',tp1,n) );
);


*** ?? after the ban, cumprod_max is constrained by the max_prod allowed (by the ban)
*** # where is max prod?
cumprod_max('gas',t,n) = cum_prodpp('gas',t,n);
loop((t,tp1)$(pre(t,tp1) and year(t) ge banstart),
  cumprod_max('gas',tp1,n) = cumprod_max('gas',t,n) + tlen(t) * max_prod('gas',t,n) );  
);

*** ?? it's the same, only with a iteration order contraint
* Maximum production allowed by constraints (DOES THIS WORK?)
loop((t,tp1)$(pre(t,tp1) and not ord(siter) eq 1),
  cumprod_max('gas',tp1,n) = cum_prodpp('gas',t,n) + tlen(t) * max_prod('gas',t,n) );  
  
*entry conditions for the while loop
allo('gas',t)=0;
countallogas=0;
glob_polyfun('gas',t) = wcum('gas',t) - cexs('gas','cum0');


*** set max iteration for the while loop
#if (external("gas"),
maxitergas=100;
#else 
#maxitergas=1;
#);


** STEP 2 - loop the allocation mechanism to account for production cuts
while (countallogas le maxitergas,
**2.6
***allocation: parameter that quantifies the amount of non economical resources that have to be extracted 
***due to low cost resources being unavailable for production
* Update price accounting for cuts  
  glob_polyfun('gas',t) = glob_polyfun('gas',t) +  allo('gas',t);
  price_sup('gas',t)$(not tfix(t) and glob_polyfun('gas',t) le 12) =  max(valuein(2005,FPRICE.l('gas',tt)), 
                                                                          poly(glob_polyfun('gas',t), cgas(polydeg,'%fossil_gas%')));
  price_sup('gas',t)$(not tfix(t) and glob_polyfun('gas',t) gt 12) =  max(valuein(2005,FPRICE.l('gas',tt)), 
                                                                          poly(12, cgas(polydeg,'%fossil_gas%')));

** 2.7 (spell error on the thesis)
** setting eq 1_new
* Calculate real cumulative production including bans as min between max cumulative production possible for ss constraints and cum production with new price
  cum_prodpp('gas',t,n) = max(0, min(cumprod_max('gas',t,n),
                                     trade_scale_gas('%fossil_gas%') * poly(price_sup('gas',t)/(twh2ej/1000),
                                                                           trade_poly_gas(polydeg,'%fossil_gas%',n) )
                                     )
                             ); 

** setting eq 2
* Ensure cumulative production in 2005 = 0
  cum_prodpp('gas',t,n)$(tfirst(t))= 0;

** setting eq 3
* This is to avoid negative production (cumulative production declining)
  loop((t,tp1)$(pre(t,tp1)),
    cum_prodpp('gas',tp1,n)=max(cum_prodpp('gas',t,n)+1e-5*tlen(t), cum_prodpp('gas',tp1,n) );
  );

** 2.8 :the difference between cumulative demand and supply, computed as regional sum of the just calculated regional cumulative production values
  allo('gas',t) = wcum('gas',t) - cexs('gas','cum0') - sum(n,cum_prodpp('gas',t,n)*1e-6);
  allo('gas',t)$( allo('gas',t) le 0) = 0;
  countallogas=countallogas+1;
);
*** end of while

** Calculate annual production
loop((t,tp1)$(pre(t,tp1)),
  prodpp('gas',t,n) = (cum_prodpp('gas',tp1,n) - cum_prodpp('gas',t,n))/tlen(t)
);

*** ?? annual productions has a different behaviour after 2100, production share remains constant
prodpp('gas',t,n)$(year(t) gt 2100) = sum(nn,Q_FUEL.l('gas',t,nn)) * valuein(2100, (prodpp('gas',tt,n)/sum(nn,prodpp('gas',tt,nn))));


*** ?? no trade case, this could be adapted as trade distruption
$ifthen.ssp set nogastrade
Q_OUT.fx('gas',t,n)$(not tfix(t)) = Q_FUEL.l('gas',t,n);
FPRICE.l('gas',t)$(not tfix(t)) = cexs('gas','scl') *  (cexs('gas','a') +  
                                       cexs('gas','c') * (wcum('gas',t) / (cexs('gas','fast') * cexs('gas','res0')))**cexs('gas','exp') ) +
                                       cexs('gas','extra'); 
$else.ssp
Q_OUT.fx('gas',t,n)$(not tfix(t)) = prodpp('gas',t,n);
FPRICE.l('gas',t)$(not tfix(t))  = price_sup('gas',t); # + sum(ssiter$(ord(ssiter) lt ord(siter) and ord(ssiter) gt 1), delta_price(run,ssiter,'gas',t));
$endif.ssp

*** ?? the mkup could be transoformed in order to simulate "dazi"

MCOST_FUEL.fx('gas',t,n)$(not tfix(t)) = FPRICE.l('gas',t) + p_mkup('gas',t,n);

*------------------------------------------------------------------------------
$elseif %phase%=='gdx_items'

trade_poly_gas
trade_scale_gas

$endif
