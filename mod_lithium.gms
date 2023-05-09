*-------------------------------------------------------------------------------       
* Module Lithium
*-------------------------------------------------------------------------------



*-------------------------------------------------------------------------------
$elseif %phase%=='include_data'

* In the phase INCLUDE_DATA you should declare and include exogenous parameters. 
* You can also modify the data loaded in data.gms
* Best practice : - create a .gdx containing those and to loading it 
*                 - this is the only phase where we should have numbers...
$gdxin '%datapath%data_mod_lithium'

parameter trade_poly_lit(polydeg,n);
$loaddc trade_poly_lit

$gdxin


*-------------------------------------------------------------------------------
$elseif %phase%=='compute_data'

* In the phase COMPUTE_DATA you should declare and compute all the parameters 
* that depend on the data loaded in the previous phase. 

**** TO CHECK 
parameter 
    cum_prodpp_l(t,n)
    prodpp(t,n);
parameters
    resgr_m(t) 'Growth rate of ultimately recoverable resources amount'
    res_m(t)   'Ultimately Recoverable Resources [TWh]'
    marg_cost_extr 
    coeff_cum_ext 
    frac_price_fast
    res_0 'Mtons'
    cum_0  'Cumulative extraction at time zero [millions TWh]'
    extra_price
    scale
    exponent
;
scale = 1;
exponent = 4;

marg_cost_extr = 27.94;
coeff_cum_ext = 7.41e-16;
frac_price_fast = 0.85;
res_0 = marg_cost_extr*scale;
extra_price = 0;
cum_0 = res_0;  


resgr_m(t) = 0.006; 
res_m(tfirst) = res_0;
loop((tnofirst(t),tm1)$pre(tm1,t),
    res_m(t) = res_m(tm1)*(1+resgr_m(tm1))**tlen(tm1);
);

* to UPDATE in core_ecnonomy?
wcum('lit',tfirst) = cum_0;

*-------------------------------------------------------------------------------
$elseif %phase%=='vars'

Variables 
    LPRICE(t,n)

    Q_OUT_L(t,n)
    Q_IN_L(t,n)
    Q_LIT(t,n)
    Q_EN_L(t,n)

    MCOST_LIT(t,n)
    COST_LIT(t,n)


*-------------------------------------------------------------------------------
$elseif %phase%=='before_solve'

LPRICE.l(t)$(not tfix(t)) = max(valuein(2005,LPRICE.l(tt)), deg0*(wcum('lit',t) - cum_0)+deg4*(wcum('lit',t) - cum_0)**4); #poly((wcum('lit',t) - cum_0), lit_price_coeff));

* Calculate cumulative production (by means of polynomial functions)
cum_prodpp_l( t,n) = max(0, trade_scale_lit*
                               poly(LPRICE.l(t)/(twh2ej/1000), trade_poly_lit(polydeg,n))
                           );

* Ensure cumulative production in 2005 = 0
cum_prodpp_l(t,n)$(tfirst(t))= 0;


* This is to avoid negative production (cumulative production declining)
loop((t,tp1)$(pre(t,tp1)),
  cum_prodpp_l(tp1,n)=max(cum_prodpp_l(t,n)+1e-5*tlen(t), cum_prodpp_l( tp1,n) )
);

** Calculate annual production
loop((t,tp1)$(pre(t,tp1)),
   prodpp_l(t,n) = (cum_prodpp_l(tp1,n) - cum_prodpp_l(t,n))/tlen(t)
);

 prodpp_l(t,n)$(year(t) > 2100) = sum(nn,Q_L.l(t,nn)) * valuein(2100, ( prodpp_l(tt,n)/sum(nn, prodpp_l(tt,nn))));

$ifthen.ssp set nolittrade
Q_OUT_L.fx(t,n)$(not tfix(t)) = 0;


*** fprice defined on the global (see phase compute_data)
LPRICE.l(t)$(not tfix(t)) = 
                        scale*(marg_cost_extr + coeff_cum_ext*
                         (wcum('lit',t) / (frac_price_fast *res_m(t)))**exponent)+ extra_price;


$else.ssp
Q_OUT_L.fx(t,n)$(not tfix(t)) =  prodpp_l(t,n);
$endif.ssp

MCOST_LIT.fx(t,n)$(not tfix(t)) = LPRICE.l(t); # + p_mkup('lit',t,n);
*-------------------------------------------------------------------------------
$elseif %phase%=='eql'

eqk
eqq_cc
eqq_fen


*-------------------------------------------------------------------------------
$elseif %phase%=='eqs'

eqq_lit(t,n)$..
    Q_LIT(t,n) =e= sum(jfed$(csi(fuel,jfed,t,n)), Q_IN_L(fuel,jfed,t,n));

    
eqcost_pes_..
    COST_LIT(t,n) =e= MCOST_LIT(t,n) * Q_LIT(f,t,n) -
                            (LPRICE.l(t) * Q_OUT_L(t,n));


eqq_en_in_%clt%(jfed,t,n)$(mapn_th('%clt%') and (sum(fuel$csi(fuel,jfed,t,n),1)))..
  Q_EN_L(jfed,t,n) =e= sum(fuel$csi(fuel,jfed,t,n), 
                      csi(fuel,jfed,t,n) * Q_IN_L(fuel,jfed,t,n));
*-------------------------------------------------------------------------------

$elseif %phase%=='gdx_items'
* List the items to be kept in the final gdx 

trade_poly_lit

resgr_m
res_m

$endif
