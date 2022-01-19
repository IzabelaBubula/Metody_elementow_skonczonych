package com.company;

public class Jakobian {
    public static double[][] jakobian;
    static double[][] jakobian_inverted;
    public static double det_j = 0.0;


    Jakobian(){
        jakobian = new double[2][2];
        jakobian_inverted = new double[2][2];
        for(int i = 0; i < 2; i ++){
            for(int j = 0; j < 2; j ++){
                jakobian[i][j] = 0.0;
                jakobian_inverted[i][j] = 0.0;
            }
        }
    }

    public static void generateJakobian(int id, Element4_2D element42d, Grid g1){
        for(int i = 0; i < Math.pow(Data.data_points, 2); i++){
            Element e1 = g1.elements.get(id);
            e1.jakobiany[i] = jakobianCount(i, id, element42d, g1);

        }
    }

    public static Jakobian jakobianCount(int nr, int id, Element4_2D element42d, Grid g1){
        Jakobian j1 = new Jakobian();
        Element e1 = g1.elements.get(id);
        //Wzór N1x1 + N2x2 + N3x3 + N4x4, N <- funkcja kształtu
        j1.jakobian[0][0] =
                element42d.d_n_d_ksi[nr][0] * g1.nodes.get(e1.ID[0]).x
                + element42d.d_n_d_ksi[nr][1] * g1.nodes.get(e1.ID[1]).x
                + element42d.d_n_d_ksi[nr][2] * g1.nodes.get(e1.ID[2]).x
                + element42d.d_n_d_ksi[nr][3] * g1.nodes.get(e1.ID[3]).x;

        j1.jakobian[1][0] =
                element42d.d_n_d_ksi[nr][0] * g1.nodes.get(e1.ID[0]).y
                + element42d.d_n_d_ksi[nr][1] * g1.nodes.get(e1.ID[1]).y
                + element42d.d_n_d_ksi[nr][2] * g1.nodes.get(e1.ID[2]).y
                + element42d.d_n_d_ksi[nr][3] * g1.nodes.get(e1.ID[3]).y;

        j1.jakobian[0][1] =
                element42d.d_n_d_eta[nr][0] * g1.nodes.get(e1.ID[0]).x
                + element42d.d_n_d_eta[nr][1] * g1.nodes.get(e1.ID[1]).x
                + element42d.d_n_d_eta[nr][2] * g1.nodes.get(e1.ID[2]).x
                + element42d.d_n_d_eta[nr][3] * g1.nodes.get(e1.ID[3]).x;

        j1.jakobian[1][1] =
                element42d.d_n_d_eta[nr][0] * g1.nodes.get(e1.ID[0]).y
                + element42d.d_n_d_eta[nr][1] * g1.nodes.get(e1.ID[1]).y
                + element42d.d_n_d_eta[nr][2] * g1.nodes.get(e1.ID[2]).y
                + element42d.d_n_d_eta[nr][3] * g1.nodes.get(e1.ID[3]).y;

        j1.det_j = (jakobian[0][0] * jakobian[1][1]) - (jakobian[1][0] * jakobian[0][1]);

        j1.jakobian_inverted[0][0] = jakobian[1][1] / det_j;
        j1.jakobian_inverted[1][1] = jakobian[0][0] / det_j;
        j1.jakobian_inverted[1][0] = - jakobian[1][0] / det_j;
        j1.jakobian_inverted[0][1] = - jakobian[0][1] /det_j;

        return j1;
    }
}
