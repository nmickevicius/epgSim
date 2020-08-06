#include "EPG.h"
#include <iostream>
#include <fstream>
#include <complex>

using namespace std;

int main()
{

    int i, j, k, n, iz;

    int inrf = 600;
    int inz = 1;
    float flip = 25.0 * M_PI / 180.0;

    float *fa = new float[inrf];
    float *phi = new float[inrf];
    float *sprof = new float[inz];
    float TR = 0.018;
    float TE = 0.006;
    float tau = 0.0028;
    float tgap = 0.0025;
    float garea = 70.0;
    float T1 = 0.8;
    float T2 = 0.08;
    float D = 0.0;

    for (i=0; i<inrf; i++) {
        fa[i] = flip;
        phi[i] = 0.0;
    }
    for (i=0; i<inz; i++)
        sprof[i] = 1.0;
    float *spfa = new float[inz];
    complex<float>* F0 = new complex<float>[inrf];

    float fgarea = 1.0E4 * (garea * 1.0E-3 * 1.0E-3); // convert gradient area from [mT/m * ms] to [G/m * s]
    float dk = 4258.0 * 2 * M_PI * fgarea;            // calculate step size in k-space [rad/m]
    float fD = D * 1.0E-6;                            // convert ADC from [mm^2/s] to [m^2/s]
    float dt;

    // create EPG object
    EPG obj = EPG(inrf, inz);

    //create Omega matrix for EPG calculations
    complex<float>*** Om = new complex<float>**[inz];
    for (i=0; i<inz; i++) {
        Om[i] = new complex<float>*[3];
        for (j=0; j<3; j++)
            Om[i][j] = new complex<float>[obj.getNumStates()];
        Om[i][2][obj.getIdxK0()] = complex<float>(1.0,0.0);
    }

    // loop over RF pulses
    for (n=0; n<inrf; n++) {

        // get slice profile for current flip angle and convert to radians
        for (i=0; i<inz; i++)
            spfa[i] = fa[n] * sprof[i];

        // apply RF pulse
        obj.rf(Om, spfa, phi[n]);

        // F0[n] = obj.F0(Om);
        // obj.relax(Om, T1, T2, TR);
        // obj.dephase(Om,1);

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


    // write data to disk
    std::ofstream outfile("F0.bin");
    //std::fstream fid = std::fstream("F0.bin", std::ios::out | std::ios::binary);
    outfile.write((char*)F0, inrf*sizeof(complex<float>));
    outfile.close();


    delete[] fa;
    delete[] phi;
    delete[] sprof;
    delete[] spfa;
    for (iz=0; iz<inz; iz++) {
        for (i=0; i<3; i++) {
            delete [] Om[iz][i];
        }
        delete [] Om[iz];
    }
    delete [] Om;

}
