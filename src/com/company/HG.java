package com.company;

public class HG {
    public static double[][] H;
    public static double[][] C;
    public static double[][] HOnly;
    public static double[][] expand;
    public static double[] P;

    public static void initialization(int node) {
        //inicjalizacja tablic
        H = new double[node][node]; //przetrzymuje H+C/delta_T <- krok czasowy
        C = new double[node][node];
        P = new double[node];
        HOnly = new double[node][node]; //przetrzymuje czyste H
        for(int i = 0; i < node; i++){ //zerowanie macierzy i wektora P
            for(int j = 0; j < node; j++){
                H[i][j] = 0.0;
                C[i][j] = 0.0;
                HOnly[i][j] = 0.0;
            }
            P[i] = 0.0;
        }
    }

    //macierz rozszerzona to H z P tzn H[H1, H2, H3] i P{P1} = [H1, H2, H3, P1]
    public static void expandMatrix(){
        //tworzenie nowej tablicy rozszerzoenj o jedną kolumnę
        expand = new double[H.length][H[0].length + 1];

        //przepisywanie wartości z macierzy H do expand
        for(int i = 0; i < H.length; i ++){
            for(int j = 0; j < H[0].length; j ++){
                expand[i][j] = H[i][j];
            }
        }
        //dopisywanie wartości P do ostatniej kolumny maceirzy H
        for(int i = 0; i < H.length; i++){
            expand[i][H.length] = P[i];
        }
    }

    /*
    do agregacji H wykorzystujemy:
        - lokalne macierze H każdego elementu
        - Hbc(część warunku brzegowego)
        - C (pojemność cieplne)

     do agregacji P wykorzystujemy:
        - P lokalne
        - C
     */
    public static void agragate(Element e1, Grid g1){
        //agregacja H
        for(int i = 0; i < 4; i ++){
            for(int j = 0; j < 4; j ++){
                HOnly[e1.ID[i]][e1.ID[j]] += e1.h[i][j];
                C[e1.ID[i]][e1.ID[j]] += e1.c[i][j];
                H[e1.ID[i]][e1.ID[j]] += e1.h[i][j] + e1.hbc[i][j] + (e1.c[i][j]/Data.delta_T);

            }
        }
        //agregacja P
        for(int i = 0; i < 4; i ++){
            double c = 0.0;

            //suma wierszy macierzy C dla kroku czasowero * t_0
            for(int j = 0; j < 4; j ++){
                c += e1.c[i][j] / Data.delta_T * g1.nodes.get(e1.ID[j]).t0;
            }
            P[e1.ID[i]] += e1.p[i] + c;
        }
    }

}
