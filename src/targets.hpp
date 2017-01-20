#pragma once
#include "core.hpp"

namespace bricolage
{

enum ScoringMethod
{
    SCORE_LINEAR = 0,
    SCORE_EXPONENTIAL = 1,
    SCORE_EXPONENTIAL_VEC = 2,
    SCORE_POWER_HAMMING = 3
};

struct cBaseTarget
{
    cBaseTarget(const cWorld_ptr &world, 
                const std::string &name, 
                int_t id,
                ScoringMethod method,
                double strength);

    virtual ~cBaseTarget() {}
    cWorld_ptr world;
    std::string name;
    int_t identifier;
    cRatesVector optimal_rates;
    ScoringMethod scoring_method;
    double strength;
    cRates weighting;

    virtual double assess(const cNetwork &net) const=0;
    void assess_networks(const cNetworkVector &networks, std::vector<double> &scores) const;
    void set_weighting(const cRates &w);
    double score_rates(const cRatesVector &rates) const;


};

struct cDefaultTarget : public cBaseTarget
{
    cDefaultTarget(const cWorld_ptr &world, 
                   const std::string &name, 
                   int_t ident=-1,
                   ScoringMethod method=SCORE_LINEAR, 
                   double strength=1.0);
    double assess(const cNetwork &net) const;
};

struct cNoisyTarget: public cBaseTarget
{
    cNoisyTarget(const cWorld_ptr &world, 
                 const std::string &name, 
                 int_t ident=-1, 
                 ScoringMethod method=SCORE_LINEAR, 
                 double strength=1.0,
                 size_t perturb_count=1,
                 double perturb_prop=1.0,
                 bool e_only=true);
    size_t perturb_count;
    double perturb_prop;
    bool env_only;
    mutable cRatesVector rates_vec;
    double assess(const cNetwork &net) const;
    void calc_perturbation(const cNetwork &net, bool env_only) const;
};

struct cTransNoisyTarget: public cBaseTarget
{
    cTransNoisyTarget(const cWorld_ptr &world, 
                 const std::string &name, 
                 int_t ident=-1, 
                 ScoringMethod method=SCORE_LINEAR, 
                 double strength=1.0);
    mutable cRatesVector rates_vec;
    double assess(const cNetwork &net) const;
};

struct cMultiTarget: public cBaseTarget
{
    cMultiTarget(const cWorld_ptr &world, 
                 const std::string &name, 
                 int_t ident=-1, 
                 ScoringMethod method=SCORE_LINEAR, 
                 double strength=1.0);
                 // const Attractor &pulses);
    double assess(const cNetwork &net) const;
    mutable cRatesVector rates_vec;
    mutable cAttractor pulses;
};

struct cSelectionModel
{
    cSelectionModel(cWorld_ptr &world, bool r=false);
    cWorld_ptr world;
    bool relative;

    bool select(const cRates &scores,
                size_t number, cIndexes &selected) const;
};


} // end namespace bricolage
// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1,/usr/local/include fdm=syntax
