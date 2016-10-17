import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

/**
 * Created by user on 6/29/2016.
 */
public class Binomial_Leap {


    private static int REACTIONS = 100;
    private static int REACTANTS = 101;
    private static int nc = 10;
    private static double EPSILON = 0.03;
    private static int INFINITY = 99999999;

    public static int minint(int [] a) {
        int min = a[0];

        for (int i = 1; i < a.length; ++i)
            if (a[i] < min)
                min = a[i];

        return min;
    }

    public static int maxindex(double [] a ) {
        int max = 0;

        for (int i = 1; i < a.length; ++i) {
            if (a[i] > a[max]) {
                max = i;
            }
        }
        return max;
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

    public static int Bernoulli_Randon_Generator(double p) {
        double r = Math.random();

        if (r <= p)
            return 1;
        else
            return 0;

    }

    public static int Binomial_Random_Generator(int n, double p) {
        int[] Y = new int[n];
        int X = 0;

        for (int i = 0; i < n; ++i) {
            Y[i] = Bernoulli_Randon_Generator(p);
            X += Y[i];
        }

        return X;
    }

    public static void main(String[] args) throws IOException {

        Scanner sc = new Scanner(System.in);

        System.out.println("Enter the number of reactants: ");
        int REACTANTS = sc.nextInt();

        System.out.println("Enter the number of reactions: ");
        int REACTIONS = sc.nextInt();

        int [][] SMatrix = new int[REACTIONS][REACTANTS];

//        System.out.println("Enter the directory of Stiochiometry Matrix : ");
//        String xyz = sc.nextLine();
        FileReader fr = new FileReader("PC_AM_MSM/Matrices/"+5+".csv");
        BufferedReader br = new BufferedReader(fr);
        String str;

        int row_index  = 0;
        while ((str = br.readLine()) != null) {
            int col_index = 0;
            char[] arr = str.toCharArray();

            for (int i = 0; i < arr.length; ++i) {

                if (arr[i] == '-') {

                    ++i;
                    SMatrix[col_index++][row_index] = -arr[i] + 48;
                }

                else if (arr[i] >= 48 && arr[i] <= 57) {
                    SMatrix[col_index++][row_index] = arr[i] - 48;
                }


                // System.out.println("Hello!");

                // c = br.read();

            }

            //System.out.println(row_index);
            ++row_index;
        }

        double t = 0.0, tao = 0.0;
        int iteration = 0;
        int [] x = new int [REACTANTS];

        double [] a = new double [REACTIONS];
        double [] k = new double [REACTIONS];

        //FileWriter writer = new FileWriter("Output_Sheet_2.csv",true);
        //FileWriter writer = new FileWriter("Binomial_Leap_Ramkrishna.csv");
        FileWriter writer = new FileWriter("Binomial_Leap_CA.csv");

        String buffer= "";

        FileReader xreader = new FileReader("conc.txt");
        BufferedReader concreader = new BufferedReader(xreader);
        for (int i = 0; i < REACTANTS; ++i) {
            buffer = concreader.readLine();
            writer.write("Substrate " + (i + 1) + ",");
            x[i] = Integer.parseInt(buffer);
        }
        concreader.close();
        writer.write("\n");


        FileReader kreader = new FileReader("k.txt");
        BufferedReader ratereader = new BufferedReader(kreader);


        System.out.println("Enter the leap size:");
        int  coarse_factor = sc.nextInt();
        for (int i = 0; i < REACTIONS; ++i)
        {
            //System.out.println("Enter the value of k" + (i + 1) + " :\n");
            //k[i] = sc.nextDouble();
            //k[i] = Double.parseDouble(ratereader.readLine());
            k[i] = 1.0;
        }


        int B_Max = 0;
        long startTime = System.currentTimeMillis();
        System.out.println("Enter the simulation time : ");
        double SIMULATION_TIME = sc.nextDouble();
        while (t < SIMULATION_TIME) {

            long startTime_Bin = System.currentTimeMillis();

            for (int i = 0; i < REACTIONS; ++i) {
                a[i] = k[i];
            }

            int[] itr = new int[REACTIONS];

            for (int i = 0; i < REACTIONS; ++i) {
                itr[i] = 0;

                for (int j = 0; j < REACTANTS; ++j) {
                    if (SMatrix[i][j] < 0) {
                        a[i] *= x[j];
                        itr[i] = itr[i] + 1;
                    }
                }
            }


            double a0 = summation(a);

            int[][] reaction_array = new int[REACTIONS][];
            int[] L = new int[REACTIONS];

            // Calculating tao_dash

            for (int j = 0; j < REACTIONS; ++j) {
                reaction_array[j] = new int[itr[j]];
                int itr2 = 0;
                for (int i = 0; i < REACTANTS && itr2 <= itr[j]; ++i) {
                    if (SMatrix[j][i] < 0) {
                        reaction_array[j][itr2] = (int) Math.floor(-x[i] / SMatrix[j][i]);
                        //System.out.println("reaction_arrray = " + reaction_array[j][itr2] );
                        ++itr2;
                    }

                }
            }

            int L0 = 0;

            for (int j = 0; j < REACTIONS; ++j) {

                L[j] = minint(reaction_array[j]);
                L0 += L[j];
                //System.out.println("L[j] = " + L[j]);
            }

            //System.out.println("L0 = " + L0);

  /*          double [][] f = new double [REACTIONS][REACTIONS];

            for (int j1 = 0 ; j1 < REACTIONS ; ++j1)
            {
                for (int j2 = 0; j2 < REACTIONS; ++j2)
                {
                    f[j2][j1] = 0.0;


                    for (int i = 0; i < REACTANTS ; ++i)
                    {
                        if (SMatrix[j1][i] < 0)
                        {
                            f[j1][j2] += (Math.abs(SMatrix[j1][i]) * a[j1] / x[i] ) * SMatrix[j2][i];
                            //System.out.println("(Math.abs(SMatrix[j1][i]) * a[j1] / x[i] ) = " +(Math.abs(SMatrix[j1][i]) * a[j1] / x[i] ));
                            //if ( f[j1][j2] != 0)
                           // System.out.println("f[j1][j2] = " + f[j1][j2]);
                           // System.out.println("SMatrix[j1][i] = " + SMatrix[j1][i] + "");
                        }
                    }


                }

            }

            double [] mean = new double [REACTIONS];

            for (int j1 = 0 ; j1 < REACTIONS ; ++j1)
            {
                mean[j1] = 0.0;

                for (int j2 = 0; j2 < REACTIONS; ++j2)
                {
                    mean[j1] += f[j1][j2] * a[j2];
                    //System.out.println("mean[j1] = " + mean[j1]);

                }
                mean[j1] = Math.abs(mean[j1]);
                //System.out.println("mean[j1] = " + mean[j1]);
            }

            int max_mean_index = maxindex(mean);
            //double max_mean = maximum ()
            System.out.println("mean[max_mean_index] = " + mean[max_mean_index]);

            tao = EPSILON * a0 / Math.abs(mean[max_mean_index]);
            //tao = 1 / a0 * Math.log(1 / Math.random());
            System.out.println("tao = " + tao);
*/
            int [] B = new int [REACTIONS];

            if (iteration == 0) {
                tao =  (double)coarse_factor / a0;/* Math.log(1 / Math.random())*/
            }
            //else {
//                tao = B_Max / a0; // Math.log(1 / Math.random());
//            }

            for (int j = 0; j < REACTIONS; ++j)
            {
                //System.out.println("tao = " + tao);
                //tao = L[j] / a[j] * Math.log(1 / Math.random());
                //System.out.println("tao = " + tao);
                double p = a[j] * tao / L[j];
                //System.out.println("p = " + p);
                B[j] = Binomial_Random_Generator(L[j], p);
                //System.out.println("B[j] = " + B[j]);
            }

            for (int i = 0; i < REACTANTS; ++i)
            {
                for (int j = 0; j < REACTIONS; ++j)
                {
                    x[i] += B[j] * SMatrix[j][i];
                }
            }

            for (int i = 0; i < REACTANTS; ++i)
            {
                writer.write(x[i] + ",");
            }
            writer.write("\n");

            B_Max = 0;

            for (int j = 0; j < REACTIONS; ++j)
            {
                B_Max += B[j];
            }
            System.out.println("B_Max = " + B_Max);

            t += tao;
            ++iteration;
            long endTime_Bin   = System.currentTimeMillis();
            long totalTime_Bin = endTime_Bin - startTime_Bin;
            System.out.println("Total Time = " + totalTime_Bin);

        }

        long endTime   = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        System.out.println("Total Time = " + totalTime);
        System.out.println("Iteration = " + iteration);

        writer.flush();
        writer.close();
    }
}




