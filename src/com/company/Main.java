package com.company;

import java.util.Arrays;

public class Main {

    public static void main(String[] args) {
	    // generowanie siatki
        Grid g1 = Grid.generate(Data.const_grid_H, Data.const_grid_B, Data.const_grid_NH, Data.const_grid_NB);
        //Show.showElements(g1);
//        Show.showAllNodes(g1);
//        Show.showElements(g1);

        //liczenie temperatury w kroku czasowym delta_T
        for(int i = (int)Data.delta_T; i <= Data.const_time; i+= Data.delta_T){

            //inicjalizacja H_global do agregacji i równań
            HG.initialization(g1.nN);

            //generowanie H i P glob oraz dla układu niestacjonarnego macierz C
            Element4_2D element42D = new Element4_2D();
            for(int k = 0; k < g1.nE; k++) {
                Element e1 = g1.elements.get(k);

                //każdy punkt całkowania 2D ma swój jakobian
                //jakobian: stosunek pół układu globalnego k lokalnego
                Jakobian.generateJakobian(k, element42D, g1);

                //generowanie macierzy H dla każdego eleementu
                e1.countH(element42D);

                //generowanie Hbc i P
                e1.countHbcP(element42D, g1);

                //generowanie C (pojemność cieplna) potrzebne do obliczenia temp w kroku delta_t
                e1.countC(element42D);

                //agregacja
                HG.agragate(e1, g1);
            }

            //łączenie H i P
            HG.expandMatrix();
//            if(i == (int)Data.delta_T){
//                System.out.println("\nH: ");
//                Show.showExpandedMatrix(HG.HOnly);
//                System.out.println("\nC: ");
//                Show.showExpandedMatrix(HG.C);
//                System.out.println("\nHC: ");
//                Show.showExpandedMatrix(HG.H);
//                System.out.println("\nP:");
//                Show.showP(HG.P);
//            }

            //liczenie wektora temp macierzy expand
            int n = HG.expand.length;
            double[] temperature = Gauss.count(HG.expand);

            //przypisanie temperatur do węzłó
            for(int j = 0; j < n; j ++){
                g1.nodes.get(j).t0 = temperature[j];
            }

            //wypisanie min i max temp dla każdej iteracji
            System.out.println("\nIn " + i + ": min: " + Arrays.stream(temperature).min().getAsDouble() + ", max: " + Arrays.stream(temperature).max().getAsDouble());

        }
    }
}
