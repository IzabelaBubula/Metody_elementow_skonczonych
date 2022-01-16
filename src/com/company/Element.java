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
    public double H[][] = new double[4][4];

    Element(int[] ID){
        jakobiany = new Jakobian[Data.data_points * Data.data_points];
        this.ID = ID;
        for(int i =0; i < 4; i++){
            for(int j = 0; j < 4; j ++){
                this.hbc[i][j] = 0;
                this.c[i][j] = 0;
                this.h[i][j] = 0;
                this.H[i][j] = 0.0;
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
            dndx[0] = (jakobiany[i].jakobian_inverted[0][0] * element42D.dn_dksi[i][0]) + (jakobiany[i].jakobian_inverted[0][1] * element42D.dn_deta[i][0]);
            dndx[1] = (jakobiany[i].jakobian_inverted[0][0] * element42D.dn_dksi[i][1]) + (jakobiany[i].jakobian_inverted[0][1] * element42D.dn_deta[i][1]);
            dndx[2] = (jakobiany[i].jakobian_inverted[0][0] * element42D.dn_dksi[i][2]) + (jakobiany[i].jakobian_inverted[0][1] * element42D.dn_deta[i][2]);
            dndx[3] = (jakobiany[i].jakobian_inverted[0][0] * element42D.dn_dksi[i][3]) + (jakobiany[i].jakobian_inverted[0][1] * element42D.dn_deta[i][3]);

            //wyliczanie dn/dy dla każdej funkcji kształtu
            dndy[0] = (jakobiany[i].jakobian_inverted[1][0] * element42D.dn_dksi[i][0]) + (jakobiany[i].jakobian_inverted[1][1] * element42D.dn_deta[i][0]);
            dndy[1] = (jakobiany[i].jakobian_inverted[1][0] * element42D.dn_dksi[i][1]) + (jakobiany[i].jakobian_inverted[1][1] * element42D.dn_deta[i][1]);
            dndy[2] = (jakobiany[i].jakobian_inverted[1][0] * element42D.dn_dksi[i][2]) + (jakobiany[i].jakobian_inverted[1][1] * element42D.dn_deta[i][2]);
            dndy[3] = (jakobiany[i].jakobian_inverted[1][0] * element42D.dn_dksi[i][3]) + (jakobiany[i].jakobian_inverted[1][1] * element42D.dn_deta[i][3]);

            //Show.showdndxdy(dndx, dndy);
            //użwyany waga * waga
//            System.out.println("Waga 1: " + element42D.gauss[1][i%Data.data_points] + " waga 2: " + element42D.gauss[1][(int)(i/Data.data_points)] + " iloczyn: " + element42D.gauss[1][i % Data.data_points] * element42D.gauss[1][(int)(i / Data.data_points)]);
            this.generateHMatrix((element42D.gauss[1][i % Data.data_points] * element42D.gauss[1][(int)(i / Data.data_points)]),
                    i, dndx, dndy);

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

    private void generateHMatrix(double GW, int pt, double[] dndx, double[] dndy) { // gaussWeight, element, dn_x, dn_y
        for(int i = 0; i < 4; i ++){
            for(int j = 0; j < 4; j ++){
                this.h[i][j] += ((dndx[i] * dndx[j]) + (dndy[i] * dndy[j])) * GW * Data.const_K * this.jakobiany[pt].det_j;

//                this.H[i][j] += this.h[i][j];
            }
        }
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
                double[] functions = new double[] {0.0, 0.0, 0.0, 0.0};

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
                if(i == 0 ) functions[i] = (1.0/4.0) * (1.0 + ((-1.0) * ksi)) * (1.0 + ((-1.0) * eta));
                else if(i == 1) functions[i] = (1.0/4.0) * (1.0 + ((1.0) * ksi)) * (1.0 + ((-1.0) * eta));
                else if(i == 2) functions[i] = (1.0/4.0) * (1.0 + ((1.0) * ksi)) * (1.0 + ((1.0) * eta));
                else if(i == 3) functions[i] = (1.0/4.0) * (1.0 + ((-1.0) * ksi)) * (1.0 + ((1.0) * eta));

                if(iTwo == 0) functions[iTwo] = (1.0/4.0) * (1.0 + ((-1.0) * ksi)) * (1.0 + ((-1.0) * eta));
                else if(iTwo == 1) functions[iTwo] = (1.0/4.0) * (1.0 + ((1.0) * ksi)) * (1.0 + ((-1.0) * eta));
                else if(iTwo == 2) functions[iTwo] = (1.0/4.0) * (1.0 + ((1.0) * ksi)) * (1.0 + ((1.0) * eta));
                else if(iTwo == 3) functions[iTwo] = (1.0/4.0) * (1.0 + ((-1.0) * ksi)) * (1.0 + ((1.0) * eta));

                // Liczenie macierzy Hbc, tutaj jest wzór integral_S = alfa({N} * {N}_trans)*dS
                // dS to det_j to długość odcinka / 2
                countHBC(element42D.gauss[1][j], functions, odcinek);

                // liczenie wektora P: integral_S =  alfa * {N} * t_0 * dS
                // tak jak poprzednio dS to det_j czyli długość odcinka / 2
                countP(j, element42D.gauss[1][j], functions, odcinek);
            }
        }
    }

    // Macierz C - czyli pojemność cieplna elementu, liczona jest tak jak macierz H dla funkcji kształtu
    // każdego elementu ze wzoru integral_V = RO * C * ({N}*{N}_trans)dV
    // przy wyliczeniu tego det_j jest to wyznacznik z obiektu jakobian, który jest liczony wyżej
    public void countC(Element4_2D element42D){
        int c = 0;
        for(int i = 0; i < Data.data_points; i ++){
            double pc1 = element42D.gauss[0][i]; //współczynnik gauss-lagrange'a

            for(int j = 0; j < Data.data_points; j ++){
                double pc2 = element42D.gauss[0][j];
                double[] shape = new double[4]; //funkcje kształtu

                shape[0] = (1.0/4.0) * (1.0 - (pc1)) * (1.0 - (pc2));
                shape[1] = (1.0/4.0) * (1.0 + (pc1)) * (1.0 - (pc2));
                shape[2] = (1.0/4.0) * (1.0 + (pc1)) * (1.0 + (pc2));
                shape[3] = (1.0/4.0) * (1.0 - (pc1)) * (1.0 + (pc2));

                generateC(element42D.gauss[1][i] * element42D.gauss[1][j], c, shape);

                c++;
            }
        }
    }

    private void generateC(double GW, int pt, double[] shapes){
        for (int i = 0; i < 4; i ++){
            for (int j = 0; j < 4; j ++){
                this.c[i][j] += shapes[i] * shapes[j] * Data.const_RO * Data.const_C * GW * this.jakobiany[pt].det_j;
            }
        }
    }

    private void countP(int wall, double weigth, double[] shapes, double det_jacobian){
        for(int i = 0; i < 4; i ++){
            this.p[i] += shapes[i] *Data.const_ALPHA * det_jacobian * Data.conts_temp_1[wall] * weigth;
        }
    }

    private void countHBC(double weight, double[] shapes, double det_jacobian){
        for(int i = 0; i < shapes.length; i++){
            this.hbc[0][i] += shapes[0] * shapes[i] * Data.const_ALPHA * weight * det_jacobian;
            this.hbc[1][i] += shapes[1] * shapes[i] * Data.const_ALPHA * weight * det_jacobian;
            this.hbc[2][i] += shapes[2] * shapes[i] * Data.const_ALPHA * weight * det_jacobian;
            this.hbc[3][i] += shapes[3] * shapes[i] * Data.const_ALPHA * weight * det_jacobian;
        }
    }

}
