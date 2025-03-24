#include "INS.h"
#include <cmath>

// Constructor
INS::INS(double deltaT, position poS, velocity veL, accelerometer accelM, gyroscope gyroS, Matrix<double> CbodyNed)
    : dt(deltaT), pos(poS), vel(veL), accel(accelM), gyro(gyroS), Cbn(CbodyNed),
    omega_ie(7.292115e-5), R0(6378137), ecc(0.0818191908425), RP(6356752.31425),
    f(1 / 298.257223563), mu(3.986004418E14) {}

// Getter implementations
Matrix<double> INS::getAlphaib() const { return alphaib; }
double INS::getMagAlpha() const { return magAlpha; }
Matrix<double> INS::getSkew() const { return alphaSkew; }
Matrix<double> INS::getWien() const { return wien; }
Matrix<double> INS::getOwenn() const { return owenn; }
Matrix<double> INS::getCbtbn() const { return Cbtbn; }
Matrix<double> INS::getAvCbn() const { return AvCbn; }
Matrix<double> INS::getfibn() const { return fibn; }
velocity INS::getVel() const { return vel; }
position INS::getPos() const { return pos; }
Matrix<double> INS::getCbn() const { return Cbn; }

// radiofcurv function
Matrix<double> INS::radiofcurv(double lat) {
    double temp = 1 - ((ecc * sin(lat)) * (ecc * sin(lat)));
    double RN = R0 * (1 - ecc * ecc) / pow(temp, 1.5);
    double RE = R0 / sqrt(temp);
    roc = Matrix<double>(2, 1);
    roc[0][0] = RN;
    roc[1][0] = RE;
    return roc;
}

// gravityNed function
Matrix<double> INS::gravityNed(double lat, double alt) {
    g = Matrix<double>(3, 1);
    g0 = 9.7803253359 * (1 + 0.001931853 * pow(sin(lat), 2)) / sqrt(1 - pow(ecc, 2) * pow(sin(lat), 2));
    g[0][0] = -8.08E-9 * alt * sin(2 * lat);
    g[1][0] = 0.0;
    g[2][0] = g0 * (1 - 2 / R0 * (1 + f * (1 - 2 * pow(sin(lat), 2)) + (pow(omega_ie, 2) * pow(R0, 2) * RP / mu)) * alt + 3 * pow(alt, 2) / pow(R0, 2));
    return g;
}

// Preliminaries function
void INS::Preliminaries() {
	std::vector<std::vector<double>> alpha = {
	{gyro.omegax * dt},
	{gyro.omegay * dt},
	{gyro.omegaz * dt}
	}; //alpha = omega*dt;

	alphaib = Matrix<double>(alpha); //Convert from vector to matrix
	Matrix<double> temp = (Matrix<double>::multiply(Matrix<double>::transpose(alphaib), alphaib));
	magAlpha = sqrt(temp[0][0]); // Magneto of alpha matrix
	alphaSkew = Matrix<double>::SkewSymetric(alphaib); //skew symetric of alpha matrix

	std::vector<std::vector<double>> nedRot = {
		{cos(pos.latidude)},
		{0.0},
		{-sin(pos.latidude)}
	};

	Matrix<double> nedrot(nedRot);
	wien = nedrot * omega_ie;
	radicurv = radiofcurv(pos.latidude);
	owenn = Matrix<double>(3, 1);
	owenn[0][0] = vel.vEast / (radicurv[1][0] + pos.altitude);
	owenn[1][0] = -vel.vNorth / (radicurv[0][0] + pos.altitude);
	owenn[2][0] = -vel.vEast * tan(pos.latidude) / (radicurv[1][0] + pos.altitude);
}

// specTrans function
void INS::specTrans() {
	if (magAlpha > 1E-8) {
		I = Matrix<double>(3, 3);
		I[0][0] = 1; I[1][1] = 1; I[2][2] = 1;

		Cbtbn = I + (1 - cos(magAlpha)) / pow(magAlpha, 2) * alphaSkew + (1 / pow(magAlpha, 2)) * (1 - (sin(magAlpha) / magAlpha)) * (alphaSkew * alphaSkew);
		AvCbn = Cbn * Cbtbn - 0.5 * Matrix<double>::SkewSymetric(wien + owenn) * Cbn * dt;
	}
	else {
		AvCbn = Cbn - 0.5 * Matrix<double>::SkewSymetric(wien + owenn) * Cbn * dt;
	}
	fibb = Matrix<double>(3, 1);
	fibb[0][0] = accel.fx; fibb[1][0] = accel.fy; fibb[2][0] = accel.fz;
	fibn = AvCbn * fibb;
}

// update function
void INS::update() {
	//Velocity Update
	Matrix<double> updV(3, 1);
	updV[0][0] = vel.vNorth;  updV[1][0] = vel.vEast; updV[2][0] = vel.vDown;
	Matrix<double> oldV(3, 1);
	oldV[0][0] = vel.vNorth;  oldV[1][0] = vel.vEast;  oldV[2][0] = vel.vDown;
	updV = updV + dt * (fibn + gravityNed(pos.latidude, pos.altitude) - Matrix<double>::SkewSymetric(owenn + 2 * wien) * updV);
	vel.vNorth = updV[0][0];
	vel.vEast = updV[1][0];
	vel.vDown = updV[2][0];

	//Position Update (Curvilinear)
	Matrix<double> oldPos(3, 1);
	oldPos[0][0] = pos.latidude; oldPos[1][0] = pos.longitude; oldPos[2][0] = pos.altitude;
	pos.altitude = oldPos[2][0] - 0.5 * dt * (oldV[2][0] + vel.vDown);
	pos.latidude = oldPos[0][0] + 0.5 * dt * ((oldV[0][0] / (radicurv[0][0] + oldPos[2][0])) + (vel.vNorth / (radicurv[0][0] + pos.altitude)));
	Matrix<double> newRadiCurv = radiofcurv(pos.latidude);
	pos.longitude = oldPos[1][0] + 0.5 * dt * ((oldV[1][0] / ((radicurv[1][0] + oldPos[2][0]) * cos(oldPos[0][0]))) + (vel.vEast / ((newRadiCurv[1][0] + pos.altitude) * cos(pos.latidude))));

	//Attitude Update
	wenn = Matrix<double>(3, 1);
	wenn[0][0] = vel.vEast / (newRadiCurv[1][0] + pos.altitude);
	wenn[1][0] = -vel.vNorth / (newRadiCurv[0][0] + pos.altitude);
	wenn[2][0] = -(vel.vEast * tan(pos.latidude)) / (newRadiCurv[1][0] + pos.altitude);

	if (magAlpha > 1.E-8) {
		Cnewold = I + (sin(magAlpha) / magAlpha) * alphaSkew + ((1 - cos(magAlpha)) / pow(magAlpha, 2)) * (alphaSkew * alphaSkew);
	}
	else {
		Cnewold = I + alphaSkew;
	}

	Cbn = (I - Matrix<double>::SkewSymetric(wien + 0.5 * owenn + 0.5 * wenn) * dt) * Cbn * Cnewold;
}