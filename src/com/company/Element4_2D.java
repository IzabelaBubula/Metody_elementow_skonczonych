package com.company;

public class Element4_2D {
    public double[][] d_n_d_ksi = new double[Data.data_points * Data.data_points][4];
    public double[][] d_n_d_eta = new double[Data.data_points * Data.data_points][4];
    public double[][] gauss;

    Element4_2D(){
        Gauss g1 = new Gauss();
        gauss = g1.createGauss(Data.data_points);
        int c = 0;
        // dn/dksi~eta~ = +/- 0.25 * (1 +/- eta~ksi~)
        for(int i = 0; i < Data.data_points; i ++){
            for(int j = 0; j < Data.data_points; j ++){
                this.d_n_d_ksi[c][0] = - (1.0/4.0) * (1.0 - this.gauss[0][i]);
                this.d_n_d_ksi[c][1] =  (1.0/4.0) * (1.0 - this.gauss[0][i]);
                this.d_n_d_ksi[c][2] =  (1.0/4.0) * (1.0 + this.gauss[0][i]);
                this.d_n_d_ksi[c][3] = - (1.0/4.0) * (1.0 + this.gauss[0][i]);

                    this.d_n_d_eta[c][0] = -(1.0 / 4.0) * (1.0 - this.gauss[0][j]);
                    this.d_n_d_eta[c][1] = -(1.0 / 4.0) * (1.0 + this.gauss[0][j]);
                    this.d_n_d_eta[c][2] = (1.0 / 4.0) * (1.0 + this.gauss[0][j]);
                    this.d_n_d_eta[c][3] = (1.0 / 4.0) * (1.0 - this.gauss[0][j]);
                c++;
            }
        }
    }
}
