#ifndef HYDROHYBRIDHLL_H
#define HYDROHYBRIDHLL_H

#include "MFP_hydro_riemann.H"
#include "MFP_hydro_hlle.H"
#include "MFP_hllc.H"
#include "MFP_hydro_defs.H"


class HydroHybridHLL : public HydroRiemannSolver
{
public:
    HydroHybridHLL();
    HydroHybridHLL(const sol::table& def);

    HydroHLLC hllc;
    HydroHLLE hlle;

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) override;

    virtual bool requires_shock_detector() const {return true;}

    virtual std::string get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

    Vector<Real> F_hlle, F_hllc;
};

#endif // HYDROHYBRIDHLL_H