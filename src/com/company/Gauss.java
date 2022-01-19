package com.company;

import java.util.ArrayList;

public class Gauss {
    public double[][] gaussTable;

    public double[][] createGauss(int i){
        //tabela współczynników i kwardratur Gauss-Lagrange'a
        if(i == 2){
            gaussTable = new double[2][2];
            gaussTable[0][0] = -1.0/Math.sqrt(3.0);
            gaussTable[0][1] = 1.0/Math.sqrt(3.0);
            gaussTable[1][0] = 1.0;
            gaussTable[1][1] = 1.0;
        } else if(i == 3){
            gaussTable = new double[2][3];
            gaussTable[0][0] = -Math.sqrt(3.0/5.0);
            gaussTable[0][1] = 0.0;
            gaussTable[0][2] = Math.sqrt(3.0/5.0);
            gaussTable[1][0] = 5.0/9.0;
            gaussTable[1][1] = 8.0/9.0;
            gaussTable[1][2] = 5.0/9.0;
        }
        return gaussTable;
    }

    public static double[] count(double[][] expanded){
        int temp = expanded.length;
        double[] x = new double[temp];

        for(int i = 0; i <= temp - 1; i ++){
            for(int j = i + 1; j < temp; j ++){
                double p = expanded[j][i]/expanded[i][i];
                for(int k = i; k <= temp; k ++){
                    expanded[j][k] -= (p * expanded[i][k]);
                }
            }
        }
        x[temp - 1] = expanded[temp - 1][temp] / expanded[temp - 1][temp - 1];
        for(int i = temp - 2; i >= 0; i--){
            double s = 0;
            for(int j = i + 1; j < temp; j ++){
                s += expanded[i][j] * x[j];
                x[i] = (expanded[i][temp]-s) / expanded[i][i];
            }
        }
        return x;
    }
}
