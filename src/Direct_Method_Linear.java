/**
 * Created by user on 6/5/2016.
 */

import java.io.FileWriter;
import java.io.IOException;

public class Direct_Method_Linear {

    private static int REACTANTS = 100;
    private static int REACTIONS = 99;
    private static double SIMULATION_TIME = 0.1;
    private static int nc = 10;
    private static double EPSILON = 0.2;
    private static int INFINITY = 99999999;

    public static int minint(int[] a) {
        int min = a[0];

        for (int i = 1; i < a.length; ++i)
            if (a[i] < min)
                min = a[i];

        return min;
    }

    public static double maximum(double[] a) {
        double max = a[0];

        for (int i = 1; i < a.length; ++i)
            if (a[i] > max)
                max = a[i];

        return max;
    }

    public static double minimum(double[] a) {
        double min = a[0];

        for (int i = 1; i < a.length; ++i)
            if (a[i] < min)
                min = a[i];

        return min;
    }

    public static double summation(double[] a) {

        double sum = 0;

        for (int i = 0; i < a.length; ++i)
            sum += a[i];

        return sum;
    }

    public static int reaction_selector(double[] a, double a0, double r2) {
        int reaction = -1;
        double sum_au = 0.0;
        double ra0 = r2 * a0;

        for (int i = 0; i < REACTIONS; ++i) {
            sum_au += a[i];
            if (ra0 < sum_au) {
                reaction = i;
                break;
            }
        }
        return reaction;

    }

    public static void main(String[] args) throws IOException {
        int[][] SMatrix = new int[REACTIONS][REACTANTS];

        for (int j = 0; j < REACTIONS; ++j)
            for (int i = 0; i < REACTANTS - 1; ++i) {
                if (i == j) {
                    SMatrix[j][i] = -1;
                    SMatrix[j][i + 1] = 1;
                }
            }



		/*
		 * while (reaction < DECOMPOSITION_REACTIONS) { for (int n = 1; n <=
		 * Math.floor(REACTANTS / 2); ++n) { for (int m = n; m <= REACTANTS - n;
		 * ++m) { --SMatrix[reaction][n - 1]; --SMatrix[reaction][m - 1];
		 * ++SMatrix[reaction][n + m - 1]; ++reaction; } } }
		 *
		 * System.out.println("Reaction = " + reaction);
		 *
		 * while (reaction < REACTIONS) { for (int n = 1; n <=
		 * Math.floor(REACTANTS / 2); ++n) { for (int m = n; m <= REACTANTS - n;
		 * ++m) { ++SMatrix[reaction][n - 1]; ++SMatrix[reaction][m - 1];
		 * --SMatrix[reaction][n + m - 1]; ++reaction; } } }
		 */


        FileWriter writer = new FileWriter("DirectMethod.csv");

        double t = 0.0, tao = 0.0;
        int iteration = 1;
        int[] x = new int[REACTANTS];
        double[] a = new double[REACTIONS];
        double[] k = new double[REACTIONS];

        for (int i = 0; i < REACTANTS; ++i) {
            writer.write("Substrate " + (i + 1) + ",");
            x[i] = 100000;
        }

        writer.write("\n");

        for (int i = 0; i < REACTIONS; ++i) {
            k[i] = 1;
        }

        long startTime = System.currentTimeMillis();

        while (t < 0.3) {
            for (int j = 0; j < REACTIONS; ++j) {
                a[j] = k[j];
            }

            int[] itr = new int[REACTIONS];

            for (int i = 0; i < REACTIONS; ++i) {
                itr[i] = 0;

                for (int j = 0; j < REACTANTS; ++j) {

                    if (SMatrix[i][j] == -1) {
                        a[i] *= x[j];
                        itr[i] = itr[i] + 1;
                    }

                    if (SMatrix[i][j] == -2) {
                        // System.out.println("hello!");
                        a[i] *= 0.5 * x[j] * (x[j] - 1);
                        itr[i] = itr[i] + 1;
                    }

                }
            }

            double a0 = summation(a);

            double r1 = Math.random();
            double r2 = Math.random();

            int reaction_index = reaction_selector(a, a0, r2);

            tao = 1 / a0 * Math.log(1 / r1);
            t += tao;

            for (int i = 0; i < REACTANTS; ++i) {
                x[i] += SMatrix[reaction_index][i];
            }

            if (reaction_index == -1)
                break;

            if((iteration % 100) == 0) {
                for (int i = 0; i < REACTANTS; ++i) {
                    writer.write(x[i] + ",");
                }
                writer.write("\n");
            }

            ++iteration;


        }



        long endTime   = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        System.out.println(totalTime);


        writer.flush();
        writer.close();
    }
}
