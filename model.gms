$OnText
French power sector financial modelling for only renewable energies as supply technologies (Offshore and Onshore wind, PV, Hydroelectricity and biogas)
and Battery and PHS (pumped hydro storage) as storage technologies,including primary and secondary reserve requirements for meteo and electricity consumption data of 2016;

Offshore and onshore wind power, Solar power and biogas capacities as well as battery storage and hydrogen (P2G) storage capacity are chosen endogenousely, while hydroelectricity lake and run-of-river and Phumped hydro storage capacities are chosen exogenousely.

Existing capacities by December 2017 are also entered as lower bound of each capacity, and investment cost for existing capacities has been considered zero.

Linear optimisation using one-hour time step with respect to Investment Cost.

By Behrang SHIRIZADEH -  March 2018
$Offtext

*-------------------------------------------------------------------------------
*                                Defining the sets
*-------------------------------------------------------------------------------
sets     h               'all hours'                     /0*8783/
         first(h)        'first hour'
         last(h)         'last hour'
         m               'month'                         /jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec/
         tec             'technology'                    /offshore, onshore, pv, river, lake, biogas, gas, phs, battery, hydrogen, methanation, nuclear/
         gen(tec)        'power plants'                  /offshore, onshore, pv, river, lake, gas, biogas, nuclear/
         vre(tec)        'variable tecs'                 /offshore, onshore, pv/
         ncomb(tec)      'non-combustible generation'    /offshore, onshore, pv, river, lake, phs, battery, nuclear/
         comb(tec)       'combustible generation techs'  /biogas, methanation/
         str(tec)        'storage technologies'          /phs, battery, methanation/
         frr(tec)        'technologies for upward FRR'   /lake, gas, phs, battery/
         nuc             'nuclear elec LCOE scenarios'   /2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500/
;
first(h) = ord(h)=1;
last(h) = ord(h)=card(h);
alias(h,hh);
*-------------------------------------------------------------------------------
*                                Inputs
*-------------------------------------------------------------------------------
$ontext
2016 had 366 days, and the hours of each month is as below:
January from 0 to 743, February from 744 to 1439, March from 1440 to 2183,
April from 2184 to 2903, May from 2904 to 3647, June from 3648 to 4367,
July from 4368 to 5111, August from 5112 to 5855, September from 5856 to 6575,
October from 6576 to 7319, November from 7320 to 8039 and December from 8040 to 8783.
$offtext
parameter month(h)  /0*743 1, 744*1439 2, 1440*2183 3, 2184*2903 4
                    2904*3647 5, 3648*4367 6, 4368*5111 7, 5112*5855 8
                    5856*6575 9, 6576*7319 10, 7320*8039 11, 8040*8783 12/
$Offlisting
parameter load_factor(vre,h) 'Production profiles of VRE'
/
$ondelim
$include  inputs/vre_profiles.csv
$offdelim
/;
parameter demand(h) 'demand profile in each hour in GW'
/
$ondelim
$include inputs/dem_electricity.csv
$offdelim
/;


Parameter lake_inflows(m) 'monthly lake inflows in GWh'
*Resource: RTE - Hourly nationwide electricity generation by sectors in 2016 for France
/
$ondelim
$include  inputs/lake_inflows.csv
$offdelim
/ ;
parameter gene_river(h) 'hourly run of river power generation in GWh'
*Resource: RTE - Hourly nationwide electricity generation by sectors in 2016 for France
/
$ondelim
$include  inputs/run_of_river.csv
$offdelim
/ ;
parameter epsilon(vre) 'additional FRR requirement for variable renewable energies because of forecast errors'
/
$ondelim
$include  inputs/reserve_requirements.csv
$offdelim
/ ;

parameter capa_ex(tec) 'existing capacities of the technologies by December 2017 in GW'
*Resource: RTE
/
$ondelim
$include  inputs/existing_capas.csv
$offdelim
/ ;
$ontext
1) Resource for the costs of power plants : EUR 26950 EN - Joint Research Centre - Institute for Energy and Transport;
"Energy Technology Reference Indicator (ETRI) projections for 2010-2050", 2014, ISBN 978-92-79-44403-6.
2) Resource for the storage costs : FCH JU (fuel cell and hydrogen joint undertaking) and 32 companies and McKinsey & Company;
"commercialization of energy storage in europe", March 2015.
$offtext
parameter capex(tec) 'annualized power capex cost in M€/GW/year'
/
$ondelim
$include  inputs/annuities.csv
$offdelim
/ ;
parameter capex_en(str) 'annualized energy capex cost of storage technologies in M€/GWh/year'
/
$ondelim
$include  inputs/str_annuities.csv
$offdelim
/ ;
parameter fOM(tec) 'annualized fixed operation and maintenance costs M€/GW/year'
/
$ondelim
$include  inputs/fO&M.csv
$offdelim
/ ;
Parameter vOM(tec) 'Variable operation and maintenance costs in M€/GWh'
/
$ondelim
$include  inputs/vO&M.csv
$offdelim
/ ;
parameter nuc_capex(nuc) 'annuity of nuclear electricity for each nuclear price scenario'
/
$ondelim
$include inputs/nuc_annuities.csv
$offdelim
/ ;
parameter nuc_fOM(nuc)   'fixed operation and maintenance cost of nuclear electricity for each nuc scenario'
/
$ondelim
$include inputs/nuc_fOM.csv
$offdelim
/ ;

parameter eta_in(str) 'charging efifciency of storage technologies' /PHS 0.95, battery 0.9, methanation 0.672/;
parameter eta_out(str) 'discharging efficiency of storage technolgoies' /PHS 0.9, battery 0.95, methanation 0.45/;

scalar pump_capa 'pumping capacity in GW' /9.3/;
scalar max_phs 'maximum volume of energy can be stored in PHS reservoir in TWh' /0.18/;
scalar max_biogas 'maxium energy can be generated by biogas in TWh' /15/;
scalar load_uncertainty 'uncertainty coefficient for hourly demand' /0.01/;
scalar delta 'load variation factor'     /0.1/;
scalar ramp_rate 'nuclear power hourly ramp up capacity' /0.15/;
scalar nuc_vOM 'variable operation and maintenance cost of nuclear electricity in €/kWh'  /0.0073/ ;
scalar cf_nuc    'load factor for nuclear power plants'  /0.8/;
scalar transmission_cost 'transmission grid interregional flow reinforcement cost for 18GW of added transmission in M€' /594/;
*-------------------------------------------------------------------------------
*                                Model
*-------------------------------------------------------------------------------
variables        GENE(tec,h)     'hourly energy generation in TWh'
                 CAPA(tec)       'overal yearly installed capacity in GW'
                 STORAGE(str,h)  'hourly electricity input of battery storage GW'
                 STORED(str,h)   'energy stored in each storage technology in GWh'
                 CAPACITY(str)   'energy volume of storage technologies in GWh'
                 RSV(frr,h)      'required upward frequency restoration reserve in GW'
                 COST            'final investment cost in b€'

positive variables GENE(tec,h),CAPA(tec),STORAGE(str,h),STORED(str,h),CAPACITY(str),RSV(frr,h);

equations        gene_vre        'variables renewable profiles generation'
                 gene_capa       'capacity and genration relation for technologies'
                 capa_frr        'capacity needed for the secondary reserve requirements'
                 combustion      'the relationship of combustible technologies'
                 storing         'the definition of stored energy in the storage options'
                 storage_const   'storage in the first hour is equal to the storage in the last hour'
                 lake_res        'constraint on water for lake reservoirs'
                 stored_cap      'maximum energy that is stored in storage units'
                 biogas_const    'maximum energy can be produced by biogas'
                 ramp_up         'the ramp up rate of nuclear power technology'
                 ramp_down       'the ramp down rate of nuclear power technology'
                 nuc_cf          'nuclear power load factor application'
                 reserves        'FRR requirement'
                 adequacy        'supply/demand relation'
                 obj             'the final objective function which is COST';

gene_vre(vre,h)..                GENE(vre,h)             =e=     CAPA(vre)*load_factor(vre,h);
gene_capa(tec,h)..               CAPA(tec)               =g=     GENE(tec,h);
capa_frr(frr,h)..                CAPA(frr)               =g=     GENE(frr,h) + RSV(frr,h);
combustion(h)..                  GENE('gas',h)           =e=     sum(comb,GENE(comb,h));
storing(h,h+1,str)..             STORED(str,h+1)         =e=     STORED(str,h) + STORAGE(str,h)*eta_in(str) - GENE(str,h)/eta_out(str);
storage_const(str,first,last)..  STORED(str,first)       =e=     STORED(str,last);
lake_res(m)..                    lake_inflows(m)         =g=     sum(h$(month(h) = ord(m)),GENE('lake',h))/1000;
stored_cap(str,h)..              STORED(str,h)           =l=     CAPACITY(str);
biogas_const..                   sum(h,GENE('biogas',h)) =l=     max_biogas*1000;
ramp_up(h,h+1)..                 GENE('nuclear',h+1)     =l=     GENE('nuclear',h) + CAPA('nuclear')*ramp_rate ;
ramp_down(h,h+1)..               GENE('nuclear',h+1)     =g=     GENE('nuclear',h) - CAPA('nuclear')*ramp_rate ;
nuc_cf..                         sum(h,GENE('nuclear',h))=l=     CAPA('nuclear')*cf_nuc*8784;
reserves(h)..                    sum(frr, RSV(frr,h))    =e=     sum(vre,epsilon(vre)*CAPA(vre))+ demand(h)*load_uncertainty*(1+delta);
adequacy(h)..                    sum(ncomb,GENE(ncomb,h))+GENE('gas',h)    =g=     demand(h) + sum(str,STORAGE(str,h));
obj..                            COST                    =e=     (sum(tec,(CAPA(tec)-capa_ex(tec))*capex(tec))+ sum(str,CAPACITY(str)*capex_en(str))+sum(tec,(CAPA(tec)*fOM(tec))) +sum((tec,h),GENE(tec,h)*vOM(tec))+transmission_cost)/1000;
*-------------------------------------------------------------------------------
*                                Initial and fixed values
*-------------------------------------------------------------------------------
CAPA.lo(tec) = capa_ex(tec);
CAPA.fx('phs') = pump_capa;
CAPA.fx('river')= capa_ex('river');
CAPA.fx('lake') = 13;
GENE.up('river',h) = gene_river(h)*capa_ex('river');
STORAGE.up('phs',h) = pump_capa;
CAPACITY.fx('phs') = max_phs*1000;
*-------------------------------------------------------------------------------
*                                Model options
*-------------------------------------------------------------------------------
model vre_nuc /all/;
$If exist vre_nuc_p.gdx execute_loadpoint 'vre_nuc_p';
option solvelink=0;
option RESLIM = 1000000;
option lp=CPLEX;
*to save the results in a GDX file for next iteration
option Savepoint=1;
*to replace the results if some already exist
option solveopt = replace;
*to stop the print of the column and row listing
option limcol = 0;
option limrow = 0;
*to stop the echo print of the input
$offlisting
*to stop the complete cross-reference list of symbols
$offsymxref
option SOLPRINT = OFF;
$onecho > cplex.opt
objrng all
rhsrng all
rngrestart rng.inc
$offecho
vre_nuc.optfile=1; vre_nuc.dictfile=2;
*-------------------------------------------------------------------------------
*                                Solve statement
*-------------------------------------------------------------------------------
parameter sumdemand      'the whole demand per year in TWh';
parameter sumgene        'the whole generation per year in TWh';
parameter gene_tec(tec) 'the overall yarly generation of each tec in TWh';
parameter lc 'load curtailment in %';
parameter lcoe 'overall LCOE of the system' ;
parameter nSTORAGE(str,h);

file summary_text /'outputs/summary.txt'/ ;
file hourly_generation /'outputs/hourly_generation.csv' / ;
file summary_csv /'outputs/summary.csv'/ ;
put summary_csv;
summary_csv.pc=5;
put 'scenario', 'cost', 'capa_offshore', 'capa_onshore', 'capa_PV', 'capa_nuclear', 'gene_offshore', 'gene_onshore', 'gene_PV', 'gene_nuclear'/;

*loop starts
loop(test,
         capex('nuclear') = nuc_capex(nuc);
         fOM('nuclear') = nuc_fOM(nuc);
         vOM('nuclear') = nuc_vOM;
         Solve vre_nuc using lp minimizing COST;
         display cost.l;
         display capa.l;
         sumdemand =  sum(h,demand(h))/1000;
         gene_tec(tec) = sum(h,GENE.l(tec,h))/1000;
         sumgene = sum((gen,h),GENE.l(gen,h))/1000 - gene_tec('gas');
         lc = (sumgene-sumdemand)/sumgene;
         lcoe = COST.l*1000/sumgene;
         display sumdemand; display sumgene; display lc; display lcoe;
         display gene_tec;
         display CAPACITY.l;
*-------------------------------------------------------------------------------
         put summary_text;
         put ' the main results for ' test.tl ' as nunclear electricity CAPEX scenario' //
         //
         'I)Overall investment cost is' cost.l 'b€' //
         //
         'II)the Renewable capacity ' //
         'PV              'capa.l('PV')'  GW'//
         'Offshore        'capa.l('offshore')'    GW'//
         'onsore          'capa.l('onshore')'     GW' //
         'run of river    'CAPA.l('river') 'GW' //
         'lake            'CAPA.l('lake') 'GW' //
         'biogas          'CAPA.l('biogas')' GW'//
         'nuclear         'capa.l('nuclear')' GW' //
         //
         'III) The energy generation of Endogenous tecs ' //
         'PV              'gene_tec('pv')'     TWh'//
         'Offshore        'gene_tec('offshore')'       TWh'//
         'Onshore         'gene_tec('onshore')'        TWh'//
         'Biogas          'gene_tec('biogas')'         TWh'//
         'Nuclear         'gene_tec('nuclear')'        TWh'//
         //
         'IV) Installed storage capacity '//
         'Battery installed capacity   ' CAPA.l('battery')' GW'//
         'Hydrogen installed capacity   ' CAPA.l('hydrogen')'  GW'//
         'methanation installed capacity  '  CAPA.l('methanation')'   GW'//
         //
         'V)Needed storage volume' //
         'Battery Storage         'CAPACITY.l("battery")'       TWh' //
         'PHS Storage             'CAPACITY.l("phs")'       TWh'//
         'methane storage         'CAPACITY.l("methanation")'   TWh'//
         //
         ;
*-------------------------------------------------------------------------------
         nSTORAGE(str,h) = 0 - STORAGE.l(str,h);
         put hourly_generation;
         hourly_generation.pc=5;
         put 'scenario', 'hour'; loop(tec, put tec.tl;) put 'demand' , 'ElecStr' , 'Pump' , 'CH4'/ ;
         loop (h,
         put test.tl, h.tl; loop(tec, put gene.l(tec,h);) put demand(h), nSTORAGE('PHS',h), nSTORAGE('battery',h), nSTORAGE('methanation',h)/
         ;);
*-------------------------------------------------------------------------------
         put summary_csv;;
         put test.tl, COST.l, CAPA.l('offshore'), CAPA.l('onshore') , CAPA.l('pv'), CAPA.l('nuclear'), gene_tec('offshore'),gene_tec('onshore'), gene_tec('pv') , gene_tec('nuclear') /
; );
*loop ends
*-------------------------------------------------------------------------------
*                                The End :D
*-------------------------------------------------------------------------------
