package com.company;

public class Element {
    public double[][] h = new double[4][4];
    public double[][] hbc = new double[4][4];
    public double[] p = new double[4];
    public double[][] c = new double[4][4];
    public Jakobian[] jakobiany;
    private double[][] dn_x;
    public double[][] dn_y;
    public int[] ID;

    Element(int[] ID){
        jakobiany = new Jakobian[Data.data_points * Data.data_points];
        this.ID = ID;
        for(int i =0; i < 4; i++){
            for(int j = 0; j < 4; j ++){
                this.hbc[i][j] = 0;
                this.c[i][j] = 0;
                this.h[i][j] = 0;
            }
            this.p[i] = 0;
        }
    }

    public void countH(Element4_2D element42D){
        double[] dndx = new double[4];
        double[] dndy = new double[4];

        for(int i = 0; i < Math.pow(Data.data_points, 2); i++){

            //wyliczanie dn/dx dla każdej funkcji kształtu
            //jakobian inverted[0] * dN/dksi~eta~ + (jakobian inverted[1] * Dn/deta~ksi~)
            dndx[0] = (jakobiany[i].jakobian_inverted[0][0] * element42D.d_n_d_ksi[i][0]) + (jakobiany[i].jakobian_inverted[0][1] * element42D.d_n_d_eta[i][0]);
            dndx[1] = (jakobiany[i].jakobian_inverted[0][0] * element42D.d_n_d_ksi[i][1]) + (jakobiany[i].jakobian_inverted[0][1] * element42D.d_n_d_eta[i][1]);
            dndx[2] = (jakobiany[i].jakobian_inverted[0][0] * element42D.d_n_d_ksi[i][2]) + (jakobiany[i].jakobian_inverted[0][1] * element42D.d_n_d_eta[i][2]);
            dndx[3] = (jakobiany[i].jakobian_inverted[0][0] * element42D.d_n_d_ksi[i][3]) + (jakobiany[i].jakobian_inverted[0][1] * element42D.d_n_d_eta[i][3]);

            //wyliczanie dn/dy dla każdej funkcji kształtu
            dndy[0] = (jakobiany[i].jakobian_inverted[1][0] * element42D.d_n_d_ksi[i][0]) + (jakobiany[i].jakobian_inverted[1][1] * element42D.d_n_d_eta[i][0]);
            dndy[1] = (jakobiany[i].jakobian_inverted[1][0] * element42D.d_n_d_ksi[i][1]) + (jakobiany[i].jakobian_inverted[1][1] * element42D.d_n_d_eta[i][1]);
            dndy[2] = (jakobiany[i].jakobian_inverted[1][0] * element42D.d_n_d_ksi[i][2]) + (jakobiany[i].jakobian_inverted[1][1] * element42D.d_n_d_eta[i][2]);
            dndy[3] = (jakobiany[i].jakobian_inverted[1][0] * element42D.d_n_d_ksi[i][3]) + (jakobiany[i].jakobian_inverted[1][1] * element42D.d_n_d_eta[i][3]);

            //Show.showdndxdy(dndx, dndy);
            for(int k = 0; k < 4; k ++){
                for(int j = 0; j < 4; j ++){
                    this.h[k][j] += ((dndx[k] * dndx[j]) + (dndy[k] * dndy[j])) * (element42D.gauss[1][i % Data.data_points] * element42D.gauss[1][(int)(i / Data.data_points)])* Data.const_K * this.jakobiany[i].det_j;

//                this.H[k][j] += this.h[k][j];
                }
            }

        }
//        System.out.println("\n\nH matrix: ");
//        for(int i = 0; i < 4; i ++){
//            for(int j = 0; j < 4; j ++){
//                System.out.print(H[i][j] + " ");
//            }
//            System.out.println();
//        }
//        System.out.println("\n\n");
    }


    //warunek brzegowy dla elem skończonego, HBC i P
    public void countHbcP(Element4_2D element42D, Grid g1){
        int i = 0, iTwo = 0;
        //po węzłach
        for(i = 0; i < 4; i ++){
            if(i + 1 > 3)
                iTwo = 0;
            else
                iTwo = i + 1;

            Node node = g1.nodes.get(ID[i]);
            Node nodeTwo = g1.nodes.get(ID[iTwo]);
//            assert node instanceof Node;
//            assert nodeTwo instanceof Node;

            //sprawdzenie czy brzeg
            //jeżeli któryś z node nie jest na brzegu to pomijamy tą iteracje
            //szukamy sytuacji gdzy oba nody są na ścianie zewnętrznej
            if(node.bc == 0 || nodeTwo.bc == 0)
                continue;

            //dł odc w kartezjańskim
            double odcinek = Math.sqrt((Math.pow((node.x - nodeTwo.x), 2)) + (Math.pow((node.y - nodeTwo.y), 2))) / 2.0;
            //System.out.println(odcinek);

            //całmkowanie po powierzchni, oblicznaie Hbc i P dla powierzchni
            double ksi, eta;

            for(int j = 0; j < Data.data_points; j++){

                //deklarowanie tablicy funkcji kształtu
                double[] shape_functions = new double[] {0.0, 0.0, 0.0, 0.0};

                //funkcja kształtu N = 1/4 * (1 - Ksi) * (1 - ETa)
                if(i == 1 || i == 3) {
                    ksi = 1;
                    eta = element42D.gauss[0][j];
                } else{
                    ksi = element42D.gauss[0][j];
                    eta = 1;
                }

                //dla Dane.data_points == 3, ksi || eta == 0 na środku ściany
                if(Data.data_points == 3){
                    if((i == 0 || i == 2) && j == (int) Math.floor(Data.data_points / 2)) {
                        ksi = 0;
                        eta = eta;
                    } else if ((i == 1 || i == 3) && j == (int)(Math.floor(Data.data_points / 2))){
                        ksi = ksi;
                        eta = 0;
                    }
                }

                //odwrócenie znaków na ścianie, rozwiązeie siatki karezjańskiej
                if(i >= 2)
                    ksi = ksi * (-1);
                else
                    ksi = ksi;
                if(i == 3 || i == 0)
                    eta = eta * (-1);
                else
                    eta = eta;

                //wyliczanie funkcji kształu dla Hbc i P(dwie funkcje dla 2 pkt na brzegu)
                if(i == 0 ) shape_functions[i] = (1.0/4.0) * (1.0 + ((-1.0) * ksi)) * (1.0 + ((-1.0) * eta));
                else if(i == 1) shape_functions[i] = (1.0/4.0) * (1.0 + ((1.0) * ksi)) * (1.0 + ((-1.0) * eta));
                else if(i == 2) shape_functions[i] = (1.0/4.0) * (1.0 + ((1.0) * ksi)) * (1.0 + ((1.0) * eta));
                else if(i == 3) shape_functions[i] = (1.0/4.0) * (1.0 + ((-1.0) * ksi)) * (1.0 + ((1.0) * eta));

                if(iTwo == 0) shape_functions[iTwo] = (1.0/4.0) * (1.0 + ((-1.0) * ksi)) * (1.0 + ((-1.0) * eta));
                else if(iTwo == 1) shape_functions[iTwo] = (1.0/4.0) * (1.0 + ((1.0) * ksi)) * (1.0 + ((-1.0) * eta));
                else if(iTwo == 2) shape_functions[iTwo] = (1.0/4.0) * (1.0 + ((1.0) * ksi)) * (1.0 + ((1.0) * eta));
                else if(iTwo == 3) shape_functions[iTwo] = (1.0/4.0) * (1.0 + ((-1.0) * ksi)) * (1.0 + ((1.0) * eta));

                for(int n = 0; n < shape_functions.length; n++){
                    this.hbc[0][n] += shape_functions[0] * shape_functions[n] * Data.const_ALPHA * element42D.gauss[1][j] * odcinek;
                    this.hbc[1][n] += shape_functions[1] * shape_functions[n] * Data.const_ALPHA * element42D.gauss[1][j] * odcinek;
                    this.hbc[2][n] += shape_functions[2] * shape_functions[n] * Data.const_ALPHA * element42D.gauss[1][j] * odcinek;
                    this.hbc[3][n] += shape_functions[3] * shape_functions[n] * Data.const_ALPHA * element42D.gauss[1][j] * odcinek;
                }

                for(int k = 0; k < 4; k ++){
                    this.p[k] += shape_functions[k] *Data.const_ALPHA * odcinek * Data.conts_temp_1[j] * element42D.gauss[1][j];
                }

            }
        }
    }

    public void countC(Element4_2D element42D){
        int temp = 0;
        for(int k = 0; k < Data.data_points; k++){
            double gauss1 = element42D.gauss[0][k]; //współczynnik gauss-lagrange'a

            for(int m = 0; m < Data.data_points; m ++){
                double gauss2 = element42D.gauss[0][m];
                double[] tab1 = new double[4]; //funkcje kształtu

                tab1[0] = (1.0/4.0) * (1.0 - (gauss1)) * (1.0 - (gauss2));
                tab1[1] = (1.0/4.0) * (1.0 + (gauss1)) * (1.0 - (gauss2));
                tab1[2] = (1.0/4.0) * (1.0 + (gauss1)) * (1.0 + (gauss2));
                tab1[3] = (1.0/4.0) * (1.0 - (gauss1)) * (1.0 + (gauss2));

                for (int i = 0; i < 4; i ++){
                    for (int j = 0; j < 4; j ++){
                        this.c[i][j] += tab1[i] * tab1[j] * Data.const_RO * Data.const_C * element42D.gauss[1][k] * element42D.gauss[1][m] * this.jakobiany[temp].det_j;
                    }
                }
                temp++;
            }
        }
    }

}
