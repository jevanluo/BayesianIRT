//
//  main.c
//  BIRT2PL
//
//  Created by jevan on 10/28/18.
//  Copyright © 2018 jevan. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// set the number of I and P
#define I_items 30
#define P_persons 2000

//structure


int main(int argc, const char * argv[]) {
//clarification
    double Runiform(double rangeLow, double rangeHigh);
    double seq(double rangeLow, double rangeHigh);
    double RNormal(double miu, double sigma, double min, double max);
    double Plogis(double x);
    


//******************************************************
// Data Preparation
//******************************************************
    srand(314159);
    

// set the fixted item and population parameters
    
    //// discrimination
    double a_disc[30];
    for (int i = 0; i < I_items; i++) {
        a_disc[i] = 1 + Runiform(-0.5, 0.5);
    }
    
    //// difficulty
    float diff_range_low  = -3.0;
    float diff_range_high =  3.0;
    float interval = (diff_range_high - diff_range_low) / (30.0 - 1.0);
    double b_diff[30];
    b_diff[29] = diff_range_high;
    int i = 29;
    while (i) {
        b_diff[i - 1] = b_diff[i] - interval;
        --i;
    }
    //// theta
    float theta      = 0;
    float sig_theta = 1.25;
    
    //// generate thetas and the I x P matrix of response probabilities
    double theta_abl[P_persons];
    for (i = 0; i < P_persons; i++) {
        theta_abl[i] = RNormal(theta, sig_theta, -5.0, 5.0);
    }
    
    double P_prob[P_persons][I_items];
    int p;
    double term1,term2;
    for (p = 0; p < P_persons; p++) {
        for (i = 0; i < I_items; i++) {
            term1 = theta_abl[p]*a_disc[i];
            term2 = a_disc[i] * b_diff[i];
            P_prob[p][i] = Plogis(term1-term2);
        }
    }
    
    //// generate the 0/1 responses U as a matrix of Bernoulli draws
    int U[P_persons][I_items];
    for (p = 0; p <P_persons; p++) {
        for(i = 0; i < I_items; i++){
            float B_draw = Runiform(0, 1);
            U[p][i] = (B_draw < P_prob[p][i])?1:0;
            printf("%d,",U[p][i]);
        }
        printf("\n");
    }
    
    
//******************************************************
// Data Preparation
//******************************************************
    
}

double Runiform(double rangeLow, double rangeHigh) {
    double myRand = rand()/(1.0 + RAND_MAX);
    double range = rangeHigh - rangeLow;
    double myRand_scaled = (myRand * range) + rangeLow;
    return myRand_scaled;
}


double Rbeta(double b)
{
    double un;
    un=0.0;
    while(un<=0.0 || un>=1.0) un=rand()*1.0/RAND_MAX;
    return 1.0-exp(1.0/b*log(un));
}

double Normal(double x,double miu,double sigma) //normal density
{
    return 1.0/sqrt(2*M_PI*sigma) * exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
}
double RNormal(double miu, double sigma, double min, double max)
{
    double x;
    double dScope;
    double y;
    do
    {
        x = Runiform(min,max);
        y = Normal(x, miu, sigma);
        dScope = Runiform(0, Normal(miu,miu,sigma));
    }while( dScope > y);
    return x;
}
double Plogis(double x){
    double z;
    z = 1.0 / (1.0 + exp(-(x - 0.0) / 1.0));
    return z;
}
//******************************************************
// MCMC Algorithm Shell
//******************************************************
//blocked.mcmc.update


//******************************************************
// Proposals
//******************************************************





//******************************************************
// Samplers: M-H
//******************************************************



