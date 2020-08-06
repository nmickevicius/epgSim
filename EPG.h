#ifndef EPG_H
#define EPG_H
#endif

#include <iostream>
#include <complex>

using namespace std;

class EPG {

    int nrf; // maximum number of states
    int nstates;
    int idxF0;
    int nz;   // number of points along slice-select axis
    complex<float>** T; // RF transition state operator

public:

    // constructor
    EPG(int, int, int);

    // apply RF transition operator
    void rf(complex<float>***, float*, float);

    // overloaded apply RF transition operator
    void rf(complex<float>***, float*, float, int, int);

    // apply relaxation and recovery
    void relax(complex<float>***, float, float, float);

    // overloaded relaxation and recovery method
    void relax(complex<float>***, float, float, float, int, int);

    // apply diffusion operator
    void diffusion(complex<float>***, float, float, float, float);

    // overloaded diffusion operator
    void diffusion(complex<float>***, float, float, float, float, int, int);

    // dephasing/crushing/shifting operator
    void dephase(complex<float>***, int);

    // overloaded dephasing operator
    void dephase(complex<float>***, int, int, int);

    // sum F0 signal over slice profile
    complex<float> F0(complex<float>***);

    // get maximum number of states
    int getNumStates();

    // get index to Z(k=0)
    int getIdxZ0();

    // get k=0 index
    int getIdxK0();

    // destructor
    ~EPG();

private:

    // complex-valued matrix multiplication
    void mtimes(complex<float>**, complex<float>**, complex<float>**, int, int, int);

    // overloaded complex-valued matrix multiplication
    void mtimes(complex<float>**, complex<float>**, complex<float>**, int, int, int, int, int);

    // build RF transition operator
    void buildRF(float, float);

    // dephasing with "positive" crusher gradient
    void dephasePlus(complex<float>**, complex<float>**, int);

    // overloaded dephasing with "positive" crusher gradient
    void dephasePlus(complex<float>**, complex<float>**, int, int, int);

    // dephasing with "negative" crusher gradient
    void dephaseMinus(complex<float>**, complex<float>**, int);

    // overloaded dephasing with "negative" crusher gradient
    void dephaseMinus(complex<float>**, complex<float>**, int, int, int);


}; // end of EPG class
