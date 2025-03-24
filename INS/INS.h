#pragma once
#ifndef INS_H
#define INS_H

#include "structures.h"
#include <vector>
#include "matrix.h"

class INS {
private:
    double dt;
    position pos;
    velocity vel;
    accelerometer accel;
    gyroscope gyro;
    Matrix<double> Cbn;
    Matrix<double> roc;
    Matrix<double> alphaib;
    double magAlpha;
    Matrix<double> alphaSkew;
    Matrix<double> wien;
    Matrix<double> radicurv;
    Matrix<double> owenn;
    Matrix<double> I;
    Matrix<double> Cbtbn;
    Matrix<double> Cnewold;
    Matrix<double> AvCbn;
    Matrix<double> fibb;
    Matrix<double> fibn;
    Matrix<double> wenn;
    Matrix<double> oldCbn;
    double g0;
    Matrix<double> g;
    double omega_ie;
    double R0;
    double ecc;
    double RP;
    double f;
    double mu;

public:
    INS(double deltaT, position poS, velocity veL, accelerometer accelM, gyroscope gyroS, Matrix<double> CbodyNed);

    // Getter Functions
    Matrix<double> getAlphaib() const;
    double getMagAlpha() const;
    Matrix<double> getSkew() const;
    Matrix<double> getWien() const;
    Matrix<double> getOwenn() const;
    Matrix<double> getCbtbn() const;
    Matrix<double> getAvCbn() const;
    Matrix<double> getfibn() const;
    velocity getVel() const;
    position getPos() const;
    Matrix<double> getCbn() const;

    // Main Methods
    Matrix<double> radiofcurv(double lat);
    Matrix<double> gravityNed(double lat, double alt);
    void Preliminaries();
    void specTrans();
    void update();
};

#endif

