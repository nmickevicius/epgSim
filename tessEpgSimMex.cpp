#include "mex.hpp"
#include "mexAdapter.hpp"
#include "EPG.h"
#include <iostream>

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {

    matlab::data::ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) {

        // check input
        checkArguments(outputs, inputs);

        // get some array sizes
        const size_t nrf = inputs[0].getNumberOfElements();
        const size_t nz = inputs[2].getNumberOfElements();
        const size_t numSegments = (size_t)(nrf / 3);

        // get input
        matlab::data::TypedArray<float> fa = std::move(inputs[0]);
        matlab::data::TypedArray<float> phi = std::move(inputs[1]);
        matlab::data::TypedArray<float> sprof = std::move(inputs[2]);
        const float TR = (float)(inputs[3][0]);
        const float TE = (float)(inputs[4][0]);
        const float tau = (float)(inputs[5][0]);
        const float garea = (float)(inputs[6][0]);
        const float tgap = (float)(inputs[7][0]);
        const float T1 = (float)(inputs[8][0]);
        const float T2 = (float)(inputs[9][0]);
        const float D = (float)(inputs[10][0]);
        const int   maxstates = (int)(inputs[11][0]);

        // create output array
        matlab::data::TypedArray<std::complex<float>> F0 = factory.createArray<std::complex<float>>({nrf});
        //matlab::data::TypedArray<std::complex<float>> F0 = factory.createArray<std::complex<float>>({numSegments,3});

        // local defines
        int i, j, k, n, iz, seg, eco, i1, i2;
        int inz = (int)(nz);
        int inrf = (int)(nrf);
        float dt, dk;
        float fgarea, fD;

        // unit conversions for diffusion
        fgarea = 1.0E4 * (garea * 1.0E-3 * 1.0E-3); // convert gradient area from [mT/m * ms] to [G/m * s]
        dk = 4258.0 * 2 * M_PI * fgarea;            // calculate step size in k-space [rad/m]
        fD = D * 1.0E-6;                            // convert ADC from [mm^2/s] to [m^2/s]

        // create standard floating point array for slice profile scaled by flip angle
        float* spfa = new float[inz];

        // create EPG object
        //EPG obj = EPG(inrf, inz);
        EPG obj = EPG(inrf, inz, maxstates); // overwrite obj if maxstates is provided



        //create Omega matrix for EPG calculations
        complex<float>*** Om = new complex<float>**[inz];
        for (i=0; i<inz; i++) {
            Om[i] = new complex<float>*[3];
            for (j=0; j<3; j++)
                Om[i][j] = new complex<float>[obj.getNumStates()];
            Om[i][2][obj.getIdxK0()] = complex<float>(1.0,0.0);
        }

        seg = 0;
        eco = 0;

        // loop over RF pulses
        for (n=0; n<inrf; n++) {

            // get slice profile for current flip angle and convert to radians
            for (i=0; i<inz; i++)
                spfa[i] = fa[n] * sprof[i];

            i1 = obj.getIdxK0() - n - 1;
            i2 = obj.getIdxK0() + n + 1;



            if (maxstates == 0) {

                // apply RF pulse
                obj.rf(Om, spfa, phi[n], i1, i2);

                // model effects between RF pulse and echo time
                if (n%3 == 1) { // F(k=0) pathway
                    obj.relax(Om, T1, T2, TE, i1, i2);
                    obj.diffusion(Om, D, TE, dk, 0.0, i1, i2);
                }
                else if (n%3 == 0) { // F(k=1) pathway (negative crusher)
                    obj.relax(Om, T1, T2, tgap, i1, i2);
                    obj.diffusion(Om, D, tgap, dk, 0.0, i1, i2);
                    obj.relax(Om, T1, T2, tau, i1, i2);
                    obj.diffusion(Om, D, tau, dk, -1.0, i1, i2);
                    obj.dephase(Om, -1, i1, i2);
                    obj.relax(Om, T1, T2, TE-tau-tgap, i1, i2);
                    obj.diffusion(Om, D, TE-tau-tgap, dk, 0.0, i1, i2);
                }
                else if (n%3 == 2) { // F(k=-1) pathway (positive crusher)
                    obj.relax(Om, T1, T2, tgap, i1, i2);
                    obj.diffusion(Om, D, tgap, dk, 0.0, i1, i2);
                    obj.relax(Om, T1, T2, tau, i1, i2);
                    obj.diffusion(Om, D, tau, dk, 1.0, i1, i2);
                    obj.dephase(Om, 1, i1, i2);
                    dt = TE - tau - tgap;
                    obj.relax(Om, T1, T2, dt, i1, i2);
                    obj.diffusion(Om, D, dt, dk, 0.0, i1, i2);
                }

                // sample signal at echo time
                F0[n] = obj.F0(Om);

                // model effects between echo time and end of TR
                if (n%3 == 1) { // F(k=0) pathway (single positive crusher)
                    dt = TR - TE - tgap - tau;
                    obj.relax(Om, T1, T2, dt, i1, i2);
                    obj.diffusion(Om, D, dt, dk, 0.0, i1, i2);
                    obj.relax(Om, T1, T2, tau, i1, i2);
                    obj.diffusion(Om, D, tau, dk, 1.0, i1, i2);
                    obj.dephase(Om,1, i1, i2);
                    obj.relax(Om, T1, T2, tgap, i1, i2);
                    obj.diffusion(Om, D, tgap, dk, 0.0, i1, i2);
                }
                else if (n%3 == 0) { // F(k=1) pathway (two positive crushers)
                    dt = TR - TE - tgap - 2*tau;
                    obj.relax(Om, T1, T2, dt, i1, i2);
                    obj.diffusion(Om, D, dt, dk, 0.0, i1, i2);
                    obj.relax(Om, T1, T2, tau, i1, i2);          // crusher 1
                    obj.diffusion(Om, D, tau, dk, 1.0, i1, i2);  // crusher 1
                    obj.dephase(Om,1, i1, i2);                   // crusher 1
                    obj.relax(Om, T1, T2, tau, i1, i2);          // crusher 2
                    obj.diffusion(Om, D, tau, dk, 1.0, i1, i2);  // crusher 2
                    obj.dephase(Om,1, i1, i2);                   // crusher 2
                    obj.relax(Om, T1, T2, tgap, i1, i2);
                    obj.diffusion(Om, D, tgap, dk, 0.0, i1, i2);
                }
                else if (n%3 == 2) { // F(k=-1) pathway (no crushers)
                    dt = TR - TE;
                    obj.relax(Om, T1, T2, dt, i1, i2);
                    obj.diffusion(Om, D, dt, dk, 0.0, i1, i2);
                }

            }
            else {

                // apply RF pulse
                obj.rf(Om, spfa, phi[n]);

                // model effects between RF pulse and echo time
                if (n%3 == 1) { // F(k=0) pathway
                    obj.relax(Om, T1, T2, TE);
                    obj.diffusion(Om, D, TE, dk, 0.0);
                }
                else if (n%3 == 0) { // F(k=1) pathway (negative crusher)
                    obj.relax(Om, T1, T2, tgap);
                    obj.diffusion(Om, D, tgap, dk, 0.0);
                    obj.relax(Om, T1, T2, tau);
                    obj.diffusion(Om, D, tau, dk, -1.0);
                    obj.dephase(Om, -1);
                    obj.relax(Om, T1, T2, TE-tau-tgap);
                    obj.diffusion(Om, D, TE-tau-tgap, dk, 0.0);
                }
                else if (n%3 == 2) { // F(k=-1) pathway (positive crusher)
                    obj.relax(Om, T1, T2, tgap);
                    obj.diffusion(Om, D, tgap, dk, 0.0);
                    obj.relax(Om, T1, T2, tau);
                    obj.diffusion(Om, D, tau, dk, 1.0);
                    obj.dephase(Om, 1);
                    dt = TE - tau - tgap;
                    obj.relax(Om, T1, T2, dt);
                    obj.diffusion(Om, D, dt, dk, 0.0);
                }

                // sample signal at echo time
                F0[n] = obj.F0(Om);

                // model effects between echo time and end of TR
                if (n%3 == 1) { // F(k=0) pathway (single positive crusher)
                    dt = TR - TE - tgap - tau;
                    obj.relax(Om, T1, T2, dt);
                    obj.diffusion(Om, D, dt, dk, 0.0);
                    obj.relax(Om, T1, T2, tau);
                    obj.diffusion(Om, D, tau, dk, 1.0);
                    obj.dephase(Om,1);
                    obj.relax(Om, T1, T2, tgap);
                    obj.diffusion(Om, D, tgap, dk, 0.0);
                }
                else if (n%3 == 0) { // F(k=1) pathway (two positive crushers)
                    dt = TR - TE - tgap - 2*tau;
                    obj.relax(Om, T1, T2, dt);
                    obj.diffusion(Om, D, dt, dk, 0.0);
                    obj.relax(Om, T1, T2, tau);          // crusher 1
                    obj.diffusion(Om, D, tau, dk, 1.0);  // crusher 1
                    obj.dephase(Om,1);                   // crusher 1
                    obj.relax(Om, T1, T2, tau);          // crusher 2
                    obj.diffusion(Om, D, tau, dk, 1.0);  // crusher 2
                    obj.dephase(Om,1);                   // crusher 2
                    obj.relax(Om, T1, T2, tgap);
                    obj.diffusion(Om, D, tgap, dk, 0.0);
                }
                else if (n%3 == 2) { // F(k=-1) pathway (no crushers)
                    dt = TR - TE;
                    obj.relax(Om, T1, T2, dt);
                    obj.diffusion(Om, D, dt, dk, 0.0);
                }

            }

            eco += 1;
            if (eco >= 3) {
                eco = 0;
                seg += 1;
            }
        }

        // set the output
        outputs[0] = std::move(F0);

        // clean up
        delete [] spfa;

        for (iz=0; iz<inz; iz++) {
            for (i=0; i<3; i++) {
                delete [] Om[iz][i];
            }
            delete [] Om[iz];
        }
        delete [] Om;

    }

    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        ArrayFactory factory;
        if (inputs[0].getNumberOfElements()%3 != 0) {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Number of RF pulses must be divisible by 3 for TESS") }));
        }
    }

};
