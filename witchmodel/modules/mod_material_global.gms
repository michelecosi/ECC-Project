*-------------------------------------------------------------------------------       
* Material Resources 
*-------------------------------------------------------------------------------

$ifthen %phase%=='sets'
* material
*set material /lit/;
*set fuel /set.material/;
*set extract(f) /lit/;
set fuel /lit/;

*-------------------------------------------------------------------------------
$elseif %phase%=='compute_data'

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



wcum('lit',tfirst) = cum_0;
*-------------------------------------------------------------------------------
$elseif %phase%=='before_solve'


FPRICE.l('lit',t)$(not tfix(t)) = 
                        scale*(marg_cost_extr + coeff_cum_ext*
                         (wcum('lit',t) / (frac_price_fast *res_m(t)))**exponent)+ extra_price;

MCOST_FUEL.fx('lit',t,n)$(not tfix(t)) = FPRICE.l('lit',t);
*COST_FUEL.fx('lit',t,n)=FPRICE.l('lit',t);
*-------------------------------------------------------------------------------
$elseif %phase%=='gdx_items'

resgr_m
res_m
wcum

$endif
