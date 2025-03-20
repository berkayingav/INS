#pragma once
#ifndef init_H
#define init_H

struct velocity {
	double vNorth;
	double vEast;
	double vDown;
};

struct position {
	double latidude;
	double longitude;
	double altitude;
};

struct accelerometer {
	double fx;
	double fy;
	double fz;
};

struct gyroscope {
	double omegax;
	double omegay;
	double omegaz;
};

#endif