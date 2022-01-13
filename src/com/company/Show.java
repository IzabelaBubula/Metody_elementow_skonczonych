package com.company;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

public class Show {

    public static void showElements(Grid g1){
        for(int i = 0 ; i < g1.nE; i++){
            System.out.println("Element nr : " + (i + 1));
            System.out.println("----- Node 1: " + g1.elements.get(i).ID[0]);
            System.out.println("----- Node 2: " + g1.elements.get(i).ID[1]);
            System.out.println("----- Node 3: " + g1.elements.get(i).ID[2]);
            System.out.println("----- Node 4: " + g1.elements.get(i).ID[3]);
            System.out.println();
        }
    }

    public static void showNode(Grid g1, int nr){
        System.out.println("Node nr : " + nr);
        System.out.println("X: " + g1.nodes.get(nr).x + ", y: " + g1.nodes.get(nr).y);
    }

    public static void showAllNodes(Grid g1){
        for(int i = 0; i < (Data.const_grid_NB * Data.const_grid_NH); i ++){
            System.out.println("Node nr : " + i);
            System.out.println("X: " + g1.nodes.get(i).x + ", y: " + g1.nodes.get(i).y);
        }
    }

    public static void showdndxdy(double[] dndx, double[] dndy){
        System.out.println("Dn/Dx: ");
        for(int i = 0; i < 4; i ++){
            System.out.print(dndx[i] + " ");
        }
        System.out.println("Dn/Dy: ");
        for(int i = 0; i < 4; i ++){
            System.out.print(dndy[i] + " ");
        }
    }

    public static void showExpandedMatrix(double[][] matrix){
        DecimalFormat numbers = new DecimalFormat("$#.00");
        for(int i = 0 ; i < matrix.length; i ++){
            for(int j = 0; j < matrix[0].length; j ++){
                System.out.printf("%.3f", matrix[i][j]);
                System.out.print(" | ");
//                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    public static void showP(double[] P){
        for(int i = 0; i < P.length; i++){
            System.out.printf("%.3f", P[i]);
            System.out.print(" | ");
        }
    }
}
