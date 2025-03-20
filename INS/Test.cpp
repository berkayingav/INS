#include "gtest/gtest.h"
#include "matrix.h"
#include "structures.h"
#include "INS.h"
#include <iomanip>
using namespace std;

std::vector<std::vector<int>> matrix3x1 = {
    {1},
    {2},
    {3}
}; 

std::vector<std::vector<int>> matrix1x3 = {
    {1,2,3}
};

std::vector<std::vector<int>> matrix3x3 = {
    {1,2,3},
    {2,4,6},
    {3,6,9}
};

std::vector<std::vector<int>> matrix3x4 = {
    {2,4,3},
    {2,5,6},
    {3,6,7}
};

std::vector<std::vector<int>> matrix3x4Add = {
    {3,6,6},
    {4,9,12},
    {6,12,16}
};

std::vector<std::vector<int>> matrix3x4Subt = {
    {1,2,0},
    {0,1,0},
    {0,0,-2}
};

std::vector<std::vector<int>> matrix3x4Trans = {
    {1,2,3},
    {2,4,6},
    {3,6,9}
};

std::vector<std::vector<int>> matrixSkew = {
    {0,-3,2},
    {3,0,-1},
    {-2,1,0}
};




// INSProject.exe --gtest_filter=MathTest.Addition     --> bash'ta bu þekilde çalýþtýrabilirsin C:\Users\berka\source\repos\INSProject\x64\Debug 
// yoluna gittikten sonra cd ile 

TEST(MatrixTest, MatrixAdd) {
    Matrix<int> mat3x3(matrix3x3);
    Matrix<int> mat3x4(matrix3x4);
    Matrix<int> mat3x4Add(matrix3x4Add);
    Matrix<int> result = Matrix<int>::add(mat3x3, mat3x4);

    EXPECT_EQ(result,mat3x4Add);
}


TEST(MatrixTest, MatrixSubtract) {
    Matrix<int> mat3x3(matrix3x3);
    Matrix<int> mat3x4(matrix3x4);
    Matrix<int> mat3x4Subt(matrix3x4Subt);
    Matrix<int> result = Matrix<int>::subtract(mat3x4, mat3x3);

    EXPECT_EQ(result, mat3x4Subt);
}

TEST(MatrixTest, MatrixMultiply) {
    Matrix<int> mat3x1(matrix3x1);
    Matrix<int> mat1x3(matrix1x3);
    Matrix<int> expected(matrix3x3);
    Matrix<int> mat3x4(matrix3x4);
    Matrix<int> result = Matrix<int>::multiply(mat3x1, mat1x3);

    EXPECT_EQ(result, expected);
}

TEST(MatrixTest, MatrixTranspose) {
    Matrix<int> mat3x3(matrix3x3);
    Matrix<int> mat3x4Trans(matrix3x4Trans);
    Matrix<int> result = Matrix<int>::transpose(mat3x3);

    EXPECT_EQ(result, mat3x4Trans);
}

TEST(MatrixTest, SkewSymetric) {
    Matrix<int> input(matrix3x1);
    Matrix<int> output(matrixSkew);
    Matrix<int> result = Matrix<int>::SkewSymetric(input);

    EXPECT_EQ(result, output);

}

//Operator Overloading Test
TEST(MatrixTest, matscal) {
    std::vector<std::vector<int>> matscal = {
    {1,2,1},
    {1,1,-1},
    {3,2,3}
    };
    std::vector<std::vector<int>> matscalres = {
    {5,10,5},
    {5,5,-5},
    {15,10,15}};
    double number = 5;
    Matrix<int> matscall(matscal);
    Matrix<int> expected(matscalres);
    Matrix<int> input = matscall*5;
    EXPECT_EQ(input, expected);

}

//INS Functions Test

double deltaT = 0.1;  // 

position pos = { 0.6973, 0.5699, 3000.0 };  // latitude, longitude, altitude for radian
velocity vel = { 50.0, 30, -3 };    // vNorth, vEast, vDown
accelerometer accel = { 1.5, -1.2, 10.2 }; // fx, fy, fz
gyroscope gyro = { 0.05, -0.08, 0.1 }; // omegax, omegay, omegaz

std::vector<std::vector<double>> CbnData = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}
};

Matrix<double> Cbn(CbnData);

INS ins(deltaT, pos, vel, accel, gyro, Cbn);


TEST(INSTest, gravityNed) {
    Matrix<double> g;
    g = ins.gravityNed(pos.latidude,pos.altitude);
    Matrix<double> result(3, 1); result[0][0] = -2.3864e-05; result[1][0] = 0;   result[2][0] = 9.7924;

    for (size_t i = 0; i < result.getRows(); ++i) {
        for (size_t j = 0; j < result.getCols(); ++j) {
            EXPECT_NEAR(g[i][j], result[i][j], 1e-4);
        }
    }

}

TEST(INSTest, radiofcurv) {
    Matrix<double> radiNE;
    radiNE = ins.radiofcurv(pos.latidude);
    Matrix<double> result(2, 1); result[0][0] = 6.3618e+06;  result[1][0] =6.3870e+06;
    for (size_t i = 0; i < result.getRows(); ++i) {
        for (size_t j = 0; j < result.getCols(); ++j) {
            EXPECT_NEAR(radiNE[i][j], result[i][j], 1e+3);
        }
    }
}

TEST(INSTest, preliminaries) {
    ins.Preliminaries();
    auto testAlpha = ins.getAlphaib();
    auto testMagAlpha = ins.getMagAlpha();
    auto testSkew = ins.getSkew();
    auto testWien = ins.getWien();
    auto testOwenn = ins.getOwenn();
    Matrix<double> alphaResult(3, 1); alphaResult[0][0] = 0.0050; alphaResult[1][0] = -0.0080; alphaResult[2][0] = 0.0100;
    double magalphaResult = 0.0137;
    Matrix<double> skewResult(3, 3); skewResult[0][0] = 0; skewResult[0][1] = -0.0100; skewResult[0][2] = -0.0080; 
    skewResult[1][0] = 0.0100; skewResult[1][1] = 0; skewResult[1][2] = -0.0050;
    skewResult[2][0] = 0.0080; skewResult[2][1] = 0.0050; skewResult[2][2] = 0;
    Matrix<double> wienResult(3, 1); wienResult[0][0] = 5.5902e-05; wienResult[1][0] = 0; wienResult[2][0] = -4.6824e-05;
    Matrix<double> owennResult(3, 1); owennResult[0][0] = 4.6949e-06; owennResult[1][0] = -7.8558e-06; owennResult[2][0] = -3.9325e-06;


    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_NEAR(testSkew[i][j], skewResult[i][j], 1e-3);
        }
    }

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 1; ++j) {
            EXPECT_NEAR(testAlpha[i][j], alphaResult[i][j], 1e-4);
            EXPECT_NEAR(testWien[i][j], wienResult[i][j], 1e-5);
            EXPECT_NEAR(testOwenn[i][j], owennResult[i][j], 1e-7);
        }
    }

    EXPECT_NEAR(testMagAlpha, magalphaResult,1e-4);
    
}


 
TEST(INSTest, spectrans) {
    ins.Preliminaries();
    ins.specTrans();
    //auto testCbtbn = ins.getCbtbn();
    auto testAvCbn = ins.getAvCbn();
    auto testfibn = ins.getfibn();
    Matrix<double> AvCbnResult(3, 3);AvCbnResult[0][0] = 1; AvCbnResult[0][1] = -0.0050; AvCbnResult[0][2] = -0.0040;
    AvCbnResult[1][0] = 0.0050; AvCbnResult[1][1] = 1; AvCbnResult[1][2] = -0.0025;
    AvCbnResult[2][0] = 0.0040; AvCbnResult[2][1] = 0.0025; AvCbnResult[2][2] = 1;
    Matrix<double> fibnResult(3, 1); fibnResult[0][0] = 1.4653; fibnResult[1][0] = -1.2181; fibnResult[2][0] = 10.2029;

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_NEAR(testAvCbn[i][j], AvCbnResult[i][j], 1e-4);
        }
    }

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 1; ++j) {
            EXPECT_NEAR(testfibn[i][j], fibnResult[i][j], 1e-4);
        }
    }
}

TEST(INSTest, update) {
    ins.Preliminaries();
    ins.specTrans();
    ins.update();
    auto testVel = ins.getVel();
    auto testPos = ins.getPos();
    auto testCbn = ins.getCbn();
    position pos = { 0.6973, 0.5699, 3000.20004 };  // latitude, longitude, altitude for radian
    velocity vel = { 50.1462, 29.8786, -1.0009 };    // vNorth, vEast, vDown
    Matrix<double> CbnResult(3, 3); CbnResult[0][0] = 0.9999; CbnResult[0][1] = -0.0100; CbnResult[0][2] = -0.0080;
    CbnResult[1][0] = 0.0100; CbnResult[1][1] = 0.9999; CbnResult[1][2] = -0.0050;
    CbnResult[2][0] = 0.0080; CbnResult[2][1] = 0.0050; CbnResult[2][2] = 1;

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_NEAR(testCbn[i][j], CbnResult[i][j], 1e-4);
        }
    }

    EXPECT_NEAR(testVel.vNorth, vel.vNorth, 1e-3);
    EXPECT_NEAR(testVel.vEast, vel.vEast, 1e-4);
    EXPECT_NEAR(testVel.vDown, vel.vDown, 1e-4);
    EXPECT_NEAR(testPos.altitude, pos.altitude, 1e-4);
    EXPECT_NEAR(testPos.latidude, pos.latidude, 1e-4);
    EXPECT_NEAR(testPos.longitude, pos.longitude, 1e-4);
}



 int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    //::testing::GTEST_FLAG(filter) = "MathTest.Addition";  Belli bir test için

    return RUN_ALL_TESTS();
} 


