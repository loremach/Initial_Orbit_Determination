#include "..\include\NutAngles.h"
#include "..\include\utils.h"

void NutAngles(double Mjd_TT, double& dpsi, double& deps){
    double T  = (Mjd_TT-(Const::MJD_J2000))/36525.0;
    double T2 = T*T;
    double T3 = T2*T;
    double rev = 360.0*3600.0; // arcsec/revolution

    double C[106][9] = {
    // l  l' F  D Om    dpsi    *T     deps     *T
        {0, 0, 0, 0, 1,-1719960,-1742,  920250,   89 },   //   1
        {0, 0, 0, 0, 2,   20620,    2,   -8950,    5 },   //   2
        {-2, 0, 2, 0, 1,     460,    0,    -240,    0},   //   3
        {2, 0,-2, 0, 0,     110,    0,       0,    0 },   //   4
        {-2, 0, 2, 0, 2,     -30,    0,      10,    0},   //   5
        {1,-1, 0,-1, 0,     -30,    0,       0,    0 },   //   6
        {0,-2, 2,-2, 1,     -20,    0,      10,    0 },   //   7
        {2, 0,-2, 0, 1,      10,    0,       0,    0 },   //   8
        {0, 0, 2,-2, 2, -131870,  -16,   57360,  -31 },   //   9
        {0, 1, 0, 0, 0,   14260,  -34,     540,   -1 },   //  10
        {0, 1, 2,-2, 2,   -5170,   12,    2240,   -6 },   //  11
        {0,-1, 2,-2, 2,    2170,   -5,    -950,    3 },   //  12
        {0, 0, 2,-2, 1,    1290,    1,    -700,    0 },   //  13
        {2, 0, 0,-2, 0,     480,    0,      10,    0 },   //  14
        {0, 0, 2,-2, 0,    -220,    0,       0,    0 },   //  15
        {0, 2, 0, 0, 0,     170,   -1,       0,    0 },   //  16
        {0, 1, 0, 0, 1,    -150,    0,      90,    0 },   //  17
        {0, 2, 2,-2, 2,    -160,    1,      70,    0 },   //  18
        {0,-1, 0, 0, 1,    -120,    0,      60,    0 },   //  19
        {-2, 0, 0, 2, 1,     -60,    0,      30,    0},   //  20
        {0,-1, 2,-2, 1,     -50,    0,      30,    0 },   //  21
        {2, 0, 0,-2, 1,      40,    0,     -20,    0 },   //  22
        {0, 1, 2,-2, 1,      40,    0,     -20,    0 },   //  23
        {1, 0, 0,-1, 0,     -40,    0,       0,    0 },   //  24
        {2, 1, 0,-2, 0,      10,    0,       0,    0 },   //  25
        {0, 0,-2, 2, 1,      10,    0,       0,    0 },   //  26
        {0, 1,-2, 2, 0,     -10,    0,       0,    0 },   //  27
        {0, 1, 0, 0, 2,      10,    0,       0,    0 },   //  28
        {-1, 0, 0, 1, 1,      10,    0,       0,    0},   //  29
        {0, 1, 2,-2, 0,     -10,    0,       0,    0 },   //  30
        {0, 0, 2, 0, 2,  -22740,   -2,    9770,   -5 },   //  31
        {1, 0, 0, 0, 0,    7120,    1,     -70,    0 },   //  32
        {0, 0, 2, 0, 1,   -3860,   -4,    2000,    0 },   //  33
        {1, 0, 2, 0, 2,   -3010,    0,    1290,   -1 },   //  34
        {1, 0, 0,-2, 0,   -1580,    0,     -10,    0 },   //  35
        {-1, 0, 2, 0, 2,    1230,    0,    -530,    0},   //  36
        {0, 0, 0, 2, 0,     630,    0,     -20,    0 },   //  37
        {1, 0, 0, 0, 1,     630,    1,    -330,    0 },   //  38
        {-1, 0, 0, 0, 1,    -580,   -1,     320,    0},   //  39
        {-1, 0, 2, 2, 2,    -590,    0,     260,    0},   //  40
        {1, 0, 2, 0, 1,    -510,    0,     270,    0 },   //  41
        {0, 0, 2, 2, 2,    -380,    0,     160,    0 },   //  42
        {2, 0, 0, 0, 0,     290,    0,     -10,    0 },   //  43
        {1, 0, 2,-2, 2,     290,    0,    -120,    0 },   //  44
        {2, 0, 2, 0, 2,    -310,    0,     130,    0 },   //  45
        {0, 0, 2, 0, 0,     260,    0,     -10,    0 },   //  46
        {-1, 0, 2, 0, 1,     210,    0,    -100,    0},   //  47
        {-1, 0, 0, 2, 1,     160,    0,     -80,    0},   //  48
        {1, 0, 0,-2, 1,    -130,    0,      70,    0 },   //  49
        {-1, 0, 2, 2, 1,    -100,    0,      50,    0},   //  50
        {1, 1, 0,-2, 0,     -70,    0,       0,    0 },   //  51
        {0, 1, 2, 0, 2,      70,    0,     -30,    0 },   //  52
        {0,-1, 2, 0, 2,     -70,    0,      30,    0 },   //  53
        {1, 0, 2, 2, 2,     -80,    0,      30,    0 },   //  54
        {1, 0, 0, 2, 0,      60,    0,       0,    0 },   //  55
        {2, 0, 2,-2, 2,      60,    0,     -30,    0 },   //  56
        {0, 0, 0, 2, 1,     -60,    0,      30,    0 },   //  57
        {0, 0, 2, 2, 1,     -70,    0,      30,    0 },   //  58
        {1, 0, 2,-2, 1,      60,    0,     -30,    0 },   //  59
        {0, 0, 0,-2, 1,     -50,    0,      30,    0 },   //  60
        {1,-1, 0, 0, 0,      50,    0,       0,    0 },   //  61
        {2, 0, 2, 0, 1,     -50,    0,      30,    0 },   //  62
        {0, 1, 0,-2, 0,     -40,    0,       0,    0 },   //  63
        {1, 0,-2, 0, 0,      40,    0,       0,    0 },   //  64
        {0, 0, 0, 1, 0,     -40,    0,       0,    0 },   //  65
        {1, 1, 0, 0, 0,     -30,    0,       0,    0 },   //  66
        {1, 0, 2, 0, 0,      30,    0,       0,    0 },   //  67
        {1,-1, 2, 0, 2,     -30,    0,      10,    0 },   //  68
        {-1,-1, 2, 2, 2,     -30,    0,      10,    0},   //  69
        {-2, 0, 0, 0, 1,     -20,    0,      10,    0},   //  70
        {3, 0, 2, 0, 2,     -30,    0,      10,    0 },   //  71
        {0,-1, 2, 2, 2,     -30,    0,      10,    0 },   //  72
        {1, 1, 2, 0, 2,      20,    0,     -10,    0 },   //  73
        {-1, 0, 2,-2, 1,     -20,    0,      10,    0},   //  74
        {2, 0, 0, 0, 1,      20,    0,     -10,    0 },   //  75
        {1, 0, 0, 0, 2,     -20,    0,      10,    0 },   //  76
        {3, 0, 0, 0, 0,      20,    0,       0,    0 },   //  77
        {0, 0, 2, 1, 2,      20,    0,     -10,    0 },   //  78
        {-1, 0, 0, 0, 2,      10,    0,     -10,    0},  //  79
        {1, 0, 0,-4, 0,     -10,    0,       0,    0 },   //  80
        {-2, 0, 2, 2, 2,      10,    0,     -10,    0},   //  81
        {-1, 0, 2, 4, 2,     -20,    0,      10,    0},   //  82
        {2, 0, 0,-4, 0,     -10,    0,       0,    0 },   //  83
        {1, 1, 2,-2, 2,      10,    0,     -10,    0 },   //  84
        {1, 0, 2, 2, 1,     -10,    0,      10,    0 },   //  85
        {-2, 0, 2, 4, 2,     -10,    0,      10,    0},   //  86
        {-1, 0, 4, 0, 2,      10,    0,       0,    0},  //  87
        {1,-1, 0,-2, 0,      10,    0,       0,    0 },   //  88
        {2, 0, 2,-2, 1,      10,    0,     -10,    0 },   //  89
        {2, 0, 2, 2, 2,     -10,    0,       0,    0 },   //  90
        {1, 0, 0, 2, 1,     -10,    0,       0,    0 },   //  91
        {0, 0, 4,-2, 2,      10,    0,       0,    0 },   //  92
        {3, 0, 2,-2, 2,      10,    0,       0,    0 },   //  93
        {1, 0, 2,-2, 0,     -10,    0,       0,    0 },   //  94
        {0, 1, 2, 0, 1,      10,    0,       0,    0 },   //  95
        {-1,-1, 0, 2, 1,      10,    0,       0,    0},   //  96
        {0, 0,-2, 0, 1,     -10,    0,       0,    0 },   //  97
        {0, 0, 2,-1, 2,     -10,    0,       0,    0 },   //  98
        {0, 1, 0, 2, 0,     -10,    0,       0,    0 },   //  99
        {1, 0,-2,-2, 0,     -10,    0,       0,    0 },   // 100
        {0,-1, 2, 0, 1,     -10,    0,       0,    0 },   // 101
        {1, 1, 0,-2, 1,     -10,    0,       0,    0 },   // 102
        {1, 0,-2, 2, 0,     -10,    0,       0,    0 },   // 103
        {2, 0, 0, 2, 0,      10,    0,       0,    0 },   // 104
        {0, 0, 2, 4, 2,     -10,    0,       0,    0 },   // 105
        {0, 1, 0, 1, 0,      10,    0,       0,    0 }    // 106
    };






    // Mean arguments of luni-solar motion

    //   l   mean anomaly of the Moon
    //   l'  mean anomaly of the Sun
    //   F   mean argument of latitude
    //   D   mean longitude elongation of the Moon from the Sun 
    //   Om  mean longitude of the ascending node

    double l  = custom_mod ( 485866.733 + (1325.0*rev +  715922.633)*T  + 31.310*T2 + 0.064*T3, rev );
    double lp = custom_mod ( 1287099.804 + (  99.0*rev + 1292581.224)*T  -  0.577*T2 - 0.012*T3, rev );
    double F  = custom_mod ( 335778.877 + (1342.0*rev +  295263.137)*T - 13.257*T2 + 0.011*T3, rev );
    double D  = custom_mod ( 1072261.307 + (1236.0*rev + 1105601.328)*T  -  6.891*T2 + 0.019*T3, rev );
    double Om = custom_mod ( 450160.280 - (   5.0*rev +  482890.539)*T +  7.455*T2 + 0.008*T3, rev );
    // Nutation in longitude and obliquity [rad]

    dpsi = 0.0;
    deps = 0.0;

    for (int i=0; i<106; i++){
        double arg  =  ( C[i][0]*l+C[i][1]*lp+C[i][2]*F+C[i][3]*D+C[i][4]*Om )/(Const::Arcs);
        dpsi = dpsi + ( C[i][5]+C[i][6]*T ) * sin(arg);
        deps = deps + ( C[i][7]+C[i][8]*T ) * cos(arg);
    }

    dpsi = 1.0e-5 * dpsi/(Const::Arcs);
    deps = 1.0e-5 * deps/(Const::Arcs);
}


