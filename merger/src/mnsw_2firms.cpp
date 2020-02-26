// The main program imple/menting algo from MNSW
// Only authority implementing only the AV MMP polcies
// Consumer Value (AV) Policies is not useful in my extension
// Always allowing or blocking does not make sense.
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <algorithm> 
#include <vector>
using namespace std;

// base parameters
struct parameters{
const double intercept = 3; // intercept of demand function A
const double slope = 26;  // slope of demand function B
const double alpha = 2./3.;  // captial share
const double theta = 1.1;  // scale of economy
// Assume uniform distribution for simplicity
const double clb = 3, cub = 6;  // standard cost of investment range
const double init_clb = 6, init_cub = 7;  // greenfield initial investment range
const double philb = 0, phiub = 1;  //  cost of proposing a merger
const double blb = 0, bub = 1; // cost for authority to block a merger
const double beta = 0.8; // discount factor
const double depre = 0.2; //depreciation
const int no_firm = 2;
vector<double> kgrid;  // length is nk
vector<double> investCost, deprTrans, propseProb, 
            mergeProb, propseCost, allowProb;
// parameters for iterations and calculations
const double crit = 1e-1;
const int nk = 21;
const  double klb=0, kub=20, ksteps = 1;
const double augMax=6, augMin=3, greenMax=7, greenMin=6;
const int nsim=12;
} prm;

// results
struct convergence{
vector<double> mc; // marginal cost
vector<double>  vfunc,vfunc_new;
vector<bool> pfuncInvest; // 1-0 vector of the firm investing or not
vector<double> pfuncCond; // conditional expected investment cost function
vector<double> pfuncAuthor; // Probability of authority aproving
vector<double> pfuncPropose; // Prob of merger proposed
vector<double> qstar; // kxk
vector<double> profit;
vector<double> pv1, pv0, pv0old; } res; 

// for post calculations


int main(){    
    srand(120); // set random seed
    init_array();
    
    // VFI 
    double err_vfi = 100.;
    for (auto ki=0; ki<res.profit.size(); ++ki){
        res.pv0[ki] = res.profit[ki]/(1-prm.beta);    
        }
    res.pv1 = res.pv0;

    while (err_vfi> prm.crit){
        res.vfunc = res.vfunc_new;
        investSimulation();
        bellman();
        vector<double> tmp1;
        for (auto i=0;i<res.vfunc.size();++i){
         tmp1.push_back(abs(res.vfunc[i]-res.vfunc_new[i]));
        }
        auto tmp2 = *max_element(tmp1.begin(), tmp1.end());
        err_vfi = tmp2;
    }
    
    
    return 0;
}

void init_array(){  // to fill in arrays
    // calculating some arrays before main loop
    // populate k_grid
    int rowloc;

    for (auto i=prm.klb; i<prm.kub; i+=prm.ksteps){
        prm.kgrid.push_back(i);
    }
    
    for (auto i=0; i<prm.nk; ++i){
    for (auto j=i; j<prm.nk; ++j){  // 
        res.qstar[rowloc] = find_opt_quantity(i, j);
        
    }
    }

    // calculate the marginal cost
    for (auto i=0; i<prm.nk; ++i){
        for (auto j=0; j<prm.nk; ++j){
            rowloc = prm.nk*j+i;
            res.mc[rowloc] = marginal_cost(i, j);
        }
    }

    // calculate 
    for (auto ki=0; ki<res.profit.size(); ++ki){
        prm.mergeProb[ki] = prm.allowProb[ki]*prm.propseProb[ki]; 
    }

}

double marginal_cost(double q , double k){  // based on the 
    double c;
    c = pow(q,((1/((1-prm.alpha)*prm.theta))-1))
    /((1-prm.alpha)*prm.theta*pow(k,(prm.alpha/(1-prm.alpha))));
    return c;
}

double find_opt_quantity(int i, int j){
    double ki = prm.kgrid[i], kj = prm.kgrid[j];
    double step = 0.05;
    double qgrid[prm.nk], pgrid[prm.nk], profit_grid[prm.nk];
    qgrid[0] = 0;
    pgrid[0] = 0;
    double optQuan, initQuan, maxQuan, otherQuan;
    int opt=0; // the optimal  location
    int rowloc = prm.nk*j+i;
    double error;
    
    while (error>prm.crit){ 
        optQuan = qgrid[opt];
        initQuan = optQuan;
        for (auto firm=0; firm<prm.no_firm; ++firm)
        //TODO: optimize by avoiding calculated combos
        // calculate the quantity and profit
        // first, define possible grid points for agents to choose
        profit_grid[0] = 0;
        for (auto i=0; i<prm.nk; ++i){
            qgrid[i] = qgrid[i-1] + step;
            pgrid[i] = prm.intercept - (qgrid[i]+optQuan)/prm.slope;
            profit_grid[i] = qgrid[i]*pgrid[i] - marginal_cost(qgrid[i],prm.kgrid[i]);
        }
        opt = distance(profit_grid, max_element(profit_grid, profit_grid+prm.nk));
    }
    res.qstar[rowloc] = qgrid[opt];
    res.profit[rowloc] = profit_grid[opt];
}

vector<double> prodVal0(){
    double merger_gain, exp_gain;
    int merged_loc, sepi, sepj, other_loc;
    res.pv0.clear();

    for (auto ki=0; ki<res.profit.size(); ++ki){
        sepj = (ki+1)%prm.nk;
        sepi = ki - prm.nk*sepj;
        merged_loc = sepj+sepi; // this the new i position, j=0
        other_loc = sepi*prm.nk + sepj;
        merger_gain = res.pv1[merged_loc]-res.pv1[ki]-res.pv1[other_loc]; 
        exp_gain = prm.mergeProb[ki]*merger_gain/2 
                - prm.propseProb[ki]*prm.propseCost[ki]/2;
        res.pv0[ki] = res.pv1[ki]+exp_gain;
    }
    return res.pv0;
}

vector<double> prodVal1(){
    for (auto ki=0; ki<res.profit.size(); ++ki){
        res.pv1[ki] = res.profit[ki] 
        - prm.investCost[ki]
        + prm.beta*prm.deprTrans[ki]*prm.kgrid[ki];
    }
    return res.pv1;
}

void bellman(){
    int rowLoc;
    double consumerVal, socialVal, producerVal;

    consumerVal = priceNew - priceOld;
    producerVal = valNew - valOld;
    socialVal = consumerVal + producerVal;
}

void investSimulation(){
    vector<double> deprVal;
    vector<vector<double>> simTransProbs, simInvestCost;
    vector<double> transProb_i, investCost_i, expVal, netBenefit; 
    vector<vector<double>> cumCost;
    vector<double> augmentDraws((prm.nk^prm.nsim)), greenDraws(prm.nsim); // of diff len

    // pre calculate depreciation values
    for (auto i=0; i<res.profit.size(); ++i){
        deprVal[i] = prm.beta*prm.deprTrans[i]*res.pv0old[i];
    }

    for (auto i=0; i<12; ++i){  // d0 12 simulations
        transProb_i.clear();
        investCost_i.clear();

        // Random Draws, we need [draws=12] of them
        generate(augmentDraws.begin(), augmentDraws.end(), cusRand(prm.augMin, prm.augMin));
        generate(greenDraws.begin(), greenDraws.end(), cusRand(prm.greenMin, prm.greenMax));

        for (auto i=0;i<(prm.nk-1);++i){
            cumCost.clear(); // cumSize chanege every loop
            for (auto di=0; di<prm.nsim; ++di){
                vector<double> *pCost  = &cumCost[di]; // per draw
                *pCost = 
            }
            

            for (auto j=0;j<prm.nk;++i){
                // Populate the expected values
                for (auto sidx=0; sidx<prm.nk; ++sidx){
                    
                }
                // Calculate the expected benefit of invest
                // find maxloc and 
            }
        }
    } 
}

double cusRand(double min, double max){
    double rr = (double) rand() / RAND_MAX;
    double r = min+ (max-min)*rr;
    return r;
} 