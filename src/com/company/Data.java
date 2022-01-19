package com.company;

public class Data {
    static final int data_points = 3;                                       // punkty całkowania
    static final double const_TO = 100.0;                                    // temaperatura początkowa
    static final double[] conts_temp_1 = {1200.0, 1200.0, 1200.0, 1200.0};  // maksymalna temperatura końcowa
    static final double const_K = 25.0;                                     // anizotropowy współczynnik przewodzenia ciepła
    static final double const_RO = 7800;                                    // gęstość
    static final double const_C = 700.0;                                    // ciepło właściwe
    static final double const_ALPHA = 300.0;                                 // współczynnik konwekcyjnej wymiany ciepła
    static final double delta_T = 50.0;                                      // zmiana czasu
    static final double const_time = 500.0;
    static final double const_grid_B = 0.1;                               // szerokość siatki
    static final double const_grid_H = 0.1;                               // wysokość siatki
    static final double const_grid_NB = 4.0;                               // ilość węzłów na szerokość (oś x)
    static final double const_grid_NH = 4.0;                               // ilość węzłów na wysokość (oś y)

}
