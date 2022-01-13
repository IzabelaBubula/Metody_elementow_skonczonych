package com.company;

public class Element4_2D {
    public double[][] dn_dksi = new double[Data.data_points * Data.data_points][4];
    public double[][] dn_deta = new double[Data.data_points * Data.data_points][4];
    public double[][] gauss;

    Element4_2D(){
        Gauss g1 = new Gauss();
        gauss = g1.createGauss(Data.data_points);
        int c = 0;
        // dn/dksi~eta~ = +/- 0.25 * (1 +/- eta~ksi~)
        for(int i = 0; i < Data.data_points; i ++){
            for(int j = 0; j < Data.data_points; j ++){
                this.dn_dksi[c][0] = - (1.0/4.0) * (1.0 - this.gauss[0][i]);
                this.dn_dksi[c][1] =  (1.0/4.0) * (1.0 - this.gauss[0][i]);
                this.dn_dksi[c][2] =  (1.0/4.0) * (1.0 + this.gauss[0][i]);
                this.dn_dksi[c][3] = - (1.0/4.0) * (1.0 + this.gauss[0][i]);

                if(c == 2) { // ze względu na złe układanie się 3 i 4 wiersza zastaosowano warunek
                    c = 3;
                    this.dn_deta[c][0] = -(1.0 / 4.0) * (1.0 - this.gauss[0][j]);
                    this.dn_deta[c][1] = -(1.0 / 4.0) * (1.0 + this.gauss[0][j]);
                    this.dn_deta[c][2] = (1.0 / 4.0) * (1.0 + this.gauss[0][j]);
                    this.dn_deta[c][3] = (1.0 / 4.0) * (1.0 - this.gauss[0][j]);
                    c = 2;
                } else if (c == 3){
                    c = 2;
                    this.dn_deta[c][0] = -(1.0 / 4.0) * (1.0 - this.gauss[0][j]);
                    this.dn_deta[c][1] = -(1.0 / 4.0) * (1.0 + this.gauss[0][j]);
                    this.dn_deta[c][2] = (1.0 / 4.0) * (1.0 + this.gauss[0][j]);
                    this.dn_deta[c][3] = (1.0 / 4.0) * (1.0 - this.gauss[0][j]);
                    c = 3;
                } else {
                    this.dn_deta[c][0] = -(1.0 / 4.0) * (1.0 - this.gauss[0][j]);
                    this.dn_deta[c][1] = -(1.0 / 4.0) * (1.0 + this.gauss[0][j]);
                    this.dn_deta[c][2] = (1.0 / 4.0) * (1.0 + this.gauss[0][j]);
                    this.dn_deta[c][3] = (1.0 / 4.0) * (1.0 - this.gauss[0][j]);
                }
                c++;
            }
        }
    }
}
