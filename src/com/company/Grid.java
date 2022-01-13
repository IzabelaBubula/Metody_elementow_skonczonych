package com.company;

import java.util.ArrayList;

public class Grid {
    public double H;
    public double B;
    public double nB;
    public double nH;
    public int nN;
    public int nE;
    public ArrayList<Node> nodes;
    public ArrayList<Element> elements;

    public Grid(double h, double b, double nB, double nH) {
        this.H = h;
        this.B = b;
        this.nB = nB;
        this.nH = nH;
        this.nN = (int) (this.nH * this.nB);
        this.nE = (int)((this.nH - 1) * (this.nB - 1));
    }

    public static Grid generate(double h, double b, double nH, double nB){
        Grid g1 = new Grid(h, b, nB, nH);
        ArrayList<Node> nodes = new ArrayList<>();
        ArrayList<Element> elements = new ArrayList<>();

        double x = 0;
        double dx = (g1.B / (g1.nB - 1));
        double dy = (g1.H / (g1.nH - 1));

        //do sprawdzenia czy nw danym node występuje warunek brzegowy(czy jest na ścainie)
        for(int i = 0; i < g1.nN; i ++){
            if(i < (g1.nH + (x * g1.nH)))
                x = x;
            else
                x = x+1;

            double x_tmp = x * dx;
            double y = (i % g1.nH) * dy;

            double bc;
            if(x_tmp == 0.0 || y == 0.0 || x_tmp == g1.B || y == g1.H)
                bc = 1;
            else
                bc = 0;

//            nodes.set(i, new Node(x_tmp, y, Data.const_TO, bc));
            nodes.add(new Node(x_tmp, y, Data.const_TO, bc));
        }

        //dopisywanie nodów do elementu
        //warunek przeniecienia nodu do kolejnej kolumny jeżeli i > nH - 1
        double e = 0;
        for (int i = 0; i < g1.nE; i ++){
            if(i > 0 && (i % (g1.nH - 1.0) == 0)) // modulo reszta z dzielenia
                e += 1;
            else
                e = e;

            int[] temp_t = new int[4];
            temp_t[0] = (int)(i + e);
            temp_t[1] = (int)(i + e + g1.nH);
            temp_t[2] = (int)(i + e + g1.nH + 1);
            temp_t[3] = (int)(i + e + 1);
//            elements.set(i, new Element(temp_t));
            elements.add(new Element(temp_t));
        }

        g1.nodes = nodes;
        g1.elements = elements;
//        nodes.removeAll(nodes);
//        elements.removeAll(elements);
        return g1;
    }

}
