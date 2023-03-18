$title      ex1 supply and demand

Scalars
scalar_demand /5/
scalar_supply /1.5/
mult_demand /2/;

Variables

DEMAND
SUPPLY
PRICE;

NONNEGATIVE VARIABLES  DEMAND,SUPPLY;

Equations

eq_1 demand
eq_2 supply
eq_3 satisfaction of demand; 


eq_1.. DEMAND =g= scalar_demand - mult_demand*PRICE;
eq_2.. SUPPLY =l= scalar_supply*PRICE;
eq_3.. DEMAND =l= SUPPLY;

MODEL supp_dem /all/;


solve supp_dem minimizing PRICE using lp


execute_unload "results_supp_dem.gdx"
