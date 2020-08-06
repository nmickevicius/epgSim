#include "EPG.h"

// overloaded constructor
EPG::EPG(int nrf, int nz, int maxstates)
{
    this->nrf = nrf;
    this->nz = nz;
    if (maxstates > 0) {
        nstates = 2 * maxstates + 1;
        idxF0 = maxstates;
    }
    else {
        nstates = 2*nrf + 1;
        idxF0 = nrf;
    }
    T = new complex<float>*[3];
    for (int i=0; i<3; i++)
        T[i] = new complex<float>[3];
}

// private method to update RF transition matrix
void EPG::buildRF(float alpha, float phi)
{
    float halpha = 0.5*alpha;
    T[0][0] = complex<float> (cos(halpha)*cos(halpha), 0.0);
    T[0][1] = complex<float> (sin(halpha)*sin(halpha)*cos(2.0*phi), sin(halpha)*sin(halpha)*sin(2*phi));
    T[0][2] = complex<float> (sin(alpha)*sin(phi), -1.0*sin(alpha)*cos(phi));
    T[1][0] = complex<float> (sin(halpha)*sin(halpha)*cos(2.0*phi), -1.0*sin(halpha)*sin(halpha)*sin(2.0*phi));
    T[1][1] = complex<float> (cos(halpha)*cos(halpha), 0.0);
    T[1][2] = complex<float> (sin(alpha)*sin(phi), sin(alpha)*cos(phi));
    T[2][0] = complex<float> (-0.5*sin(alpha)*sin(phi), -0.5*sin(alpha)*cos(phi));
    T[2][1] = complex<float> (-0.5*sin(alpha)*sin(phi), 0.5*sin(alpha)*cos(phi));
    T[2][2] = complex<float> (cos(alpha), 0.0);
}

// private method to perform matrix multiplication: C=A*B
void EPG::mtimes(complex<float>** C, complex<float>** A, complex<float>** B, int m, int n, int p)
{
    int i, j, k;
    for (i=0; i<m; i++) {
        for (j=0; j<p; j++) {
            C[i][j] = complex<float> (0.0,0.0);
            for (k=0; k<n; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
    }
}

// overloaded private method to perform matrix multiplication on a subset of columns of second input: C(:,i1:i2) = A*B(:,i1:i2)
void EPG::mtimes(complex<float>** C, complex<float>** A, complex<float>** B, int m, int n, int p, int i1, int i2)
{
    int i, j, k;
    for (i=0; i<m; i++) {
        for (j=i1; j<i2; j++) {
            C[i][j] = complex<float> (0.0,0.0);
            for (k=0; k<n; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
    }
}

// private method for dephasing due to positive crusher
void EPG::dephasePlus(complex<float>** Om, complex<float>** dummy, int nshift)
{
    int i, j;
    for (i=0; i<3; i++) {
        for (j=0; j<nstates; j++)
            dummy[i][j] = Om[i][j];
    }
    for (i=0; i<(nstates-1); i++) { // NJM: added -1 here
        Om[0][i+nshift] = dummy[0][i]; // shift F+ states up
        Om[1][i] = dummy[0][i+nshift]; // shift F- states down
    }
    Om[0][idxF0] = conj(Om[1][idxF0]);
}

// overloaded private method for dephasing due to positive crusher for subset of columns of Om
void EPG::dephasePlus(complex<float>** Om, complex<float>** dummy, int nshift, int i1, int i2)
{
    int i, j;
    for (i=0; i<3; i++) {
        for (j=i1; j<i2; j++)
            dummy[i][j] = Om[i][j];
    }
    for (i=i1; i<i2; i++) {
        Om[0][i+nshift] = dummy[0][i]; // shift F+ states up
        Om[1][i] = dummy[0][i+nshift]; // shift F- states down
    }
    Om[0][idxF0] = conj(Om[1][idxF0]);
}

// private method for dephasing due to negative crusher
void EPG::dephaseMinus(complex<float>** Om, complex<float>** dummy, int nshift)
{
    int i, j;
    for (i=0; i<3; i++) {
        for (j=0; j<nstates; j++)
            dummy[i][j] = Om[i][j];
    }
    for (i=(nstates-1); i>=1; i--) { // NJM: changed >=0 to >=1 here
        Om[0][i-nshift] = dummy[0][i]; // shift F+ states down
        Om[1][i] = dummy[1][i-nshift]; // shift F- states up
    }
    Om[0][idxF0] = conj(Om[1][idxF0]);
}

// private method for dephasing due to negative crusher
void EPG::dephaseMinus(complex<float>** Om, complex<float>** dummy, int nshift, int i1, int i2)
{
    int i, j;
    for (i=0; i<3; i++) {
        for (j=i1; j<i2; j++)
            dummy[i][j] = Om[i][j];
    }
    for (i=i2; i>=i1; i--) {
        Om[0][i-nshift] = dummy[0][i]; // shift F+ states down
        Om[1][i] = dummy[1][i-nshift]; // shift F- states up
    }
    Om[0][idxF0] = conj(Om[1][idxF0]);
}


// public method for applying RF operator
void EPG::rf(complex<float>*** Om, float* zalpha, float phi)
{
    int iz, i, j;
    complex<float>** out = new complex<float>*[3];
    for (i=0; i<3; i++)
        out[i] = new complex<float>[nstates];
    for (iz=0; iz<nz; iz++) {
        buildRF(zalpha[iz], phi);              // build RF operator
        mtimes(out, T, Om[iz], 3, 3, nstates); // apply RF operator
        for (i=0; i<3; i++) {
            for (j=0; j<nstates; j++)
                Om[iz][i][j] = out[i][j];
        }
    }
    for (i=0; i<3; i++)
        delete [] out[i];
    delete [] out;
}

// overloaded public method for applying RF operator to subset of columns of Om
void EPG::rf(complex<float>*** Om, float* zalpha, float phi, int i1, int i2)
{
    int iz, i, j;
    complex<float>** out = new complex<float>*[3];
    for (i=0; i<3; i++)
        out[i] = new complex<float>[nstates];
    for (iz=0; iz<nz; iz++) {
        buildRF(zalpha[iz], phi);              // build RF operator
        mtimes(out, T, Om[iz], 3, 3, nstates, i1, i2); // apply RF operator
        for (i=0; i<3; i++) {
            for (j=0; j<nstates; j++)
                Om[iz][i][j] = out[i][j];
        }
    }
    for (i=0; i<3; i++)
        delete [] out[i];
    delete [] out;
}

// public method for applying relaxation
void EPG::relax(complex<float>*** Om, float T1, float T2, float dt)
{
    complex<float> r1 (exp(-dt/T1),0.0);
    complex<float> r2 (exp(-dt/T2),0.0);
    complex<float> b1 (1.0 - exp(-dt/T1),0.0);
    int iz, row, col, idx;
    for (iz=0; iz<nz; iz++) {
        for (col=0; col<nstates; col++) {
            Om[iz][0][col] *= r2;
            Om[iz][1][col] *= r2;
            Om[iz][2][col] *= r1;
        }
        Om[iz][2][idxF0] += b1;
    }
}

// overloaded public method for applying relaxation to a subset of columns of Om
void EPG::relax(complex<float>*** Om, float T1, float T2, float dt, int i1, int i2)
{
    complex<float> r1 (exp(-dt/T1),0.0);
    complex<float> r2 (exp(-dt/T2),0.0);
    complex<float> b1 (1.0 - exp(-dt/T1),0.0);
    int iz, row, col, idx;
    for (iz=0; iz<nz; iz++) {
        for (col=i1; col<i2; col++) {
            Om[iz][0][col] *= r2;
            Om[iz][1][col] *= r2;
            Om[iz][2][col] *= r1;
        }
        Om[iz][2][idxF0] += b1;
    }
}

// public method for applying diffusion
void EPG::diffusion(complex<float>*** Om, float D, float t, float dk, float gradOnOff)
{

    if (D == 0.0)
        return;

    int iz, ik, j, idx;
    float k, kvalz, kvalfp, kvalfm;
    float bvalz, bvalfp, bvalfm;
    complex<float> expbzD, expbfpD, expbfmD;

    // loop over points in slice profile
    for (iz=0; iz<nz; iz++) {
        for (ik=0; ik<nstates; ik++) {

            k = (float)(ik - idxF0);

            kvalz = k * dk;
            kvalfp = (k + 0.5*gradOnOff)*dk;
            kvalfm = (-1.0*k + 0.5*gradOnOff)*dk;

            bvalz = kvalz * kvalz * t;
            bvalfp = (kvalfp*kvalfp + gradOnOff*dk*dk/12.0) * t;
            bvalfm = (kvalfm*kvalfm + gradOnOff*dk*dk/12.0) * t;

            expbzD = complex<float> (exp(-bvalz * D),0.0);
            expbfpD = complex<float> (exp(-bvalfp * D),0.0);
            expbfmD = complex<float> (exp(-bvalfm * D),0.0);

            Om[iz][0][ik] *= expbfpD;
            Om[iz][1][ik] *= expbfmD;
            Om[iz][2][ik] *= expbzD;

        }
    }

}

// overloaded public method for applying diffusion to subset of columns of Om
void EPG::diffusion(complex<float>*** Om, float D, float t, float dk, float gradOnOff, int i1, int i2)
{

    if (D == 0.0)
        return;

    int iz, ik, j, idx;
    float k, kvalz, kvalfp, kvalfm;
    float bvalz, bvalfp, bvalfm;
    complex<float> expbzD, expbfpD, expbfmD;

    // loop over points in slice profile
    for (iz=0; iz<nz; iz++) {
        for (ik=i1; ik<i2; ik++) {

            k = (float)(ik - idxF0);

            kvalz = k * dk;
            kvalfp = (k + 0.5*gradOnOff)*dk;
            kvalfm = (-1.0*k + 0.5*gradOnOff)*dk;

            bvalz = kvalz * kvalz * t;
            bvalfp = (kvalfp*kvalfp + gradOnOff*dk*dk/12.0) * t;
            bvalfm = (kvalfm*kvalfm + gradOnOff*dk*dk/12.0) * t;

            expbzD = complex<float> (exp(-bvalz * D),0.0);
            expbfpD = complex<float> (exp(-bvalfp * D),0.0);
            expbfmD = complex<float> (exp(-bvalfm * D),0.0);

            Om[iz][0][ik] *= expbfpD;
            Om[iz][1][ik] *= expbfmD;
            Om[iz][2][ik] *= expbzD;

        }
    }

}


// public method for applying dephasing operator
void EPG::dephase(complex<float>*** Om, int nshift)
{
    int iz, i;
    complex<float>**dummy = new complex<float>*[3];
    for (i=0; i<3; i++)
        dummy[i] = new complex<float>[nstates];
    for (iz=0; iz<nz; iz++) {
        if (nshift > 1)
            dephasePlus(Om[iz], dummy, abs(nshift));
        else
            dephaseMinus(Om[iz], dummy, abs(nshift));
    }
    for (i=0; i<3; i++)
        delete [] dummy[i];
    delete [] dummy;
}

// overloaded public method for applying dephasing operator applied to subset of columns of Om
void EPG::dephase(complex<float>*** Om, int nshift, int i1, int i2)
{
    if (i1 <= 0)
        i1 = 1;
    if (i2 >= nstates)
        i2 = nstates - 1;
    int iz, i;
    complex<float>**dummy = new complex<float>*[3];
    for (i=0; i<3; i++)
        dummy[i] = new complex<float>[nstates];
    for (iz=0; iz<nz; iz++) {
        if (nshift > 1)
            dephasePlus(Om[iz], dummy, abs(nshift), i1, i2);
        else
            dephaseMinus(Om[iz], dummy, abs(nshift), i1, i2);
    }
    for (i=0; i<3; i++)
        delete [] dummy[i];
    delete [] dummy;
}

// public method to get measurable F+(k=0) signal
complex<float> EPG::F0(complex<float>*** Om)
{
    int iz;
    complex<float> sig (0.0,0.0);
    for (iz=0; iz<nz; iz++)
        sig += Om[iz][0][idxF0];
    return sig;
} 

// public method to get number of states
int EPG::getNumStates()
{
    return nstates;
}

// public method to get index to Z(k=0)
int EPG::getIdxZ0()
{
    return 2*nstates + idxF0;
}

int EPG::getIdxK0()
{
    return idxF0;
}

// destructor
EPG::~EPG()
{
    for (int i=0; i<3; i++)
        delete [] T[i];
    delete [] T;
}
