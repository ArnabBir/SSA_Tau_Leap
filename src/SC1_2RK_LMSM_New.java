import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

/**
 * Created by user on 6/13/2016.
 */
public class SC1_2RK_LMSM_New {

    //    private static int REACTIONS = 99;
//    private static int REACTANTS = 100;
    private static int nc = 10;
    private static double EPSILON = 0.03;

    public static int minint(int[] a) {
        int min = a[0];

        for (int i = 1; i < a.length; ++i)
            if (a[i] < min)
                min = a[i];

        return min;
    }

    public static int maxindex(double[] a) {
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



        //System.out.println(SMatrix[4][4]);
//        for (int j = 0; j < REACTIONS; ++j)
//            for (int i = 0; i < REACTANTS - 1; ++i) {
//                if (i == j) {
//                    SMatrix[j][i] = -1;
//                    SMatrix[j][i + 1] = 1;
//                }
//            }
        //  System.out.println(SMatrix[4][5]);

        Scanner sc = new Scanner(System.in);

        System.out.println("Enter the number of reactants: ");
        int REACTANTS = sc.nextInt();

        System.out.println("Enter the number of reactions: ");
        int REACTIONS = sc.nextInt();

        int [][] SMatrix = new int[REACTIONS][REACTANTS];

//        System.out.println("Enter the directory of Stiochiometry Matrix : ");
//        String xyz = sc.nextLine();
        FileReader fr = new FileReader("Strongly_Coupled_1_Matrix.csv");
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
        int[] x = new int[REACTANTS];
        int[][] xp = new int[4][REACTANTS];

        double[][] ap = new double[4][REACTIONS];
        double[] a = new double[REACTIONS];
        double[] k = new double[REACTIONS];

        //FileWriter writer = new FileWriter("Output_Sheet_2.csv",true);
        //FileWriter writer = new FileWriter("Binomial_Leap_Ramkrishna.csv");
        FileWriter writer = new FileWriter("SC1_2RK_4AMM_PC_Binomial_Leap.csv");


        for (int i = 0; i < REACTANTS; ++i) {
            writer.write("Substrate " + (i + 1) + ",");
            x[i] = 100000;
        }

        writer.write("\n");

        for (int i = 0; i < REACTIONS; ++i) {
            //System.out.println("Enter the value of k" + (i + 1) + " :\n");
            //k[i] = sc.nextDouble();
            k[i] = 1;
        }

        int steps = 0;
        int B_Max = 0;
        long startTime = System.currentTimeMillis();

        while (steps < 10000000) {

            for (int i = 0; i < REACTIONS; ++i) {
                a[i] = k[i];
            }

            //System.out.println("hello!!");

            int[] itr = new int[REACTIONS];

            for (int i = 0; i < REACTIONS; ++i) {
                itr[i] = 0;

                for (int j = 0; j < REACTANTS; ++j) {
                    if (SMatrix[i][j] == -1) {
                        a[i] *= x[j];
                        itr[i] = itr[i] + 1;
                    }
                    if( SMatrix[i][j] == -2){
                        a[i] *= 0.5 * x[j] * (x[j] - 1);
                        itr[i] = itr[i] + 1;
                    }
                }
            }

            double a0 = summation(a);

            int[][] reaction_array = new int[REACTIONS][];
            int[] L = new int[REACTIONS];

            for (int j = 0; j < REACTIONS; ++j) {

                reaction_array[j] = new int[itr[j]];
                //System.out.println("itr[j] = " + itr[j] );

                int itr2 = 0;

                for (int i = 0; i < REACTANTS && itr2 <= itr[j] ; ++i) {
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
                //System.out.println("L[j] = " + L[j]);
                L0 += L[j];
            }
            int L_min = minint(L);

            //System.out.println("L0 = " + L0);

//            double[][] f = new double[REACTIONS][REACTIONS];
//
//            for (int j1 = 0; j1 < REACTIONS; ++j1) {
//                for (int j2 = 0; j2 < REACTIONS; ++j2) {
//                    f[j2][j1] = 0.0;
//
//
//                    for (int i = 0; i < REACTANTS; ++i) {
//                        if (SMatrix[j1][i] < 0) {
//                            f[j1][j2] += (Math.abs(SMatrix[j1][i]) * a[j1] / x[i]) * SMatrix[j2][i];
//                            //System.out.println("(Math.abs(SMatrix[j1][i]) * a[j1] / x[i] ) = " +(Math.abs(SMatrix[j1][i]) * a[j1] / x[i] ));
//                            //if ( f[j1][j2] != 0)
//                            // System.out.println("f[j1][j2] = " + f[j1][j2]);
//                            // System.out.println("SMatrix[j1][i] = " + SMatrix[j1][i] + "");
//                        }
//                    }
//                }
//            }
//
//            double[] mean = new double[REACTIONS];
//            double[] variance = new double[REACTIONS];
//
//            for (int j1 = 0; j1 < REACTIONS; ++j1) {
//                mean[j1] = 0.0;
//                variance[j1] = 0.0;
//
//                for (int j2 = 0; j2 < REACTIONS; ++j2) {
//                    mean[j1] += f[j1][j2] * a[j2];
//                    variance[j1] += f[j1][j2] * f[j1][j2] * a[j2];
//                    //System.out.println("mean[j1] = " + mean[j1]);
//
//                }
//                mean[j1] = Math.abs(mean[j1]);
//                //System.out.println("mean[j1] = " + mean[j1]);
//            }
//
//            int max_mean_index = maxindex(mean);
//            int max_variance_index = maxindex(variance);
//            //double max_mean = maximum ()
//            //System.out.println("mean[max_mean_index] = " + mean[max_mean_index]);
//
//            tao = EPSILON * a0 / Math.abs(mean[max_mean_index]);
//            //tao = 1 / a0 * Math.log(1 / Math.random());
//            System.out.println("tao = " + tao);
//            System.out.println("a0 * tao = " + a0 * tao);
//            //System.out.println("3 * Math.sqrt(variance[max_variance_index] * tao = " + 3 * Math.sqrt(variance[max_variance_index] * tao));


////for(int j = 0; j < REACTIONS; ++j)
//            if (a0 * tao + 3 * Math.sqrt(a0 * tao ) >=  L0) {
////                tao = (Math.sqrt(9 * variance[max_variance_index] + 4 * mean[max_mean_index] * L_min) - 3 * Math.sqrt(variance[max_variance_index])) / (2 * mean[max_mean_index]);
////                tao *= tao;
//                tao = (Math.sqrt( 9 + 4 * L0) - 3) / (2 * a0) ;
//            }

            int[] B = new int[REACTIONS];
            int[] Bp = new int[REACTIONS];

//
            if (iteration == 0) {
                tao = 1.8E6 / a0;
            }
//            else {
//                tao = B_Max / a0; // Math.log(1 / Math.random());
//            }

            //System.out.println("tao_final = " + tao);

            if (iteration < 4) {

                long startTrap = System.currentTimeMillis();

                double[] l = new double[REACTIONS];

                for (int j = 0; j < REACTIONS; ++j) {

                    double p = 2.00 * a[j] * tao / (L[j] * 3.00);
                    //System.out.println("p = " + p);
                    B[j] = Binomial_Random_Generator(L[j], p);
                    //System.out.println("B[j] = " + B[j]);
                }

                for (int i = 0; i < REACTANTS; ++i) {
                    xp[iteration][i] = x[i];
                    for (int j = 0; j < REACTIONS; ++j) {
                        xp[iteration][i] += B[j] * SMatrix[j][i];
                    }
                }


                for (int j = 0; j < REACTIONS; ++j) {
                    ap[iteration][j] = k[j];
                    for (int i = 0; i < REACTANTS; ++i) {
                        if (SMatrix[j][i] == -1) {
                            ap[iteration][j] *= xp[iteration][i];
                        }
                        if (SMatrix[j][i] == -2){
                            ap[iteration][j] *= 0.5 * xp[iteration][i] * ( xp[iteration][i] - 1);
                        }
                    }

//                    l[j] = 2 * ap[iteration][j] - a[j];
//
//                    if (l[j] < 0)
//                        l[j] = 0;
                }

                for (int j = 0; j < REACTIONS; ++j) {
                    //System.out.println("tao = " + tao);
                    //tao = L[j] / a[j] * Math.log(1 / Math.random());
                    //System.out.println("tao = " + tao);
                    double p_dash = ( 1.00 / 4.00 * a[j] + 3.00 / 4.00 * ap[iteration][j]) * tao / L[j] ;
                    //System.out.println("p = " + p);
                    Bp[j] = Binomial_Random_Generator(L[j], p_dash);
                    //System.out.println("B[j] = " + B[j]);
                }

                for (int i = 0; i < REACTANTS; ++i) {
//                    x[i] = xp[iteration][i];
                    for (int j = 0; j < REACTIONS; ++j) {
                        x[i] += Bp[j] * SMatrix[j][i];
                    }
                }

                for (int j = 0; j < REACTIONS; ++j) {
                    a[j] = k[j];
                    for (int i = 0; i < REACTANTS; ++i) {
                        if (SMatrix[j][i] == -1) {
                            a[j] *= x[i];
                        }
                        if (SMatrix[j][i] == -2){
                            a[j] *= 0.5 * x[i] * (x[i] - 1);
                        }
                    }
                }

                for (int j = 0; j < REACTIONS; ++j) {
                    ap[iteration][j] = a[j];
                    //xp[iteration][i] = x[i];
                }

                for (int i = 0; i < REACTANTS; ++i) {
                 //   ap[iteration][j] = a[j];
                    xp[iteration][i] = x[i];
                }

                long endTrap = System.currentTimeMillis();
                long totalTrap = endTrap - startTrap;
                System.out.println("totalTrap = " + totalTrap);


            }

            else {
                // tao = 0.5;
                long startMulti = System.currentTimeMillis();

                for (int i = 0; i < REACTIONS; i++) {
                    ap[3][i] = a[i];
                    //xp[2][i] = x[i];
                }

                int[] x_predicted = new int[REACTANTS];

                for (int j = 0; j < REACTIONS; ++j) {


                    //System.out.println("product[j] = " + x_predicted[j] * tao);
                    double mean_equivalent = (55.00 / 24.00 * ap[3][j] * tao - 59.00 / 24.00 * ap[2][j] * tao + 37.00 / 24.00 * ap[1][j] * tao - 9.00 / 24.00 * ap[0][j]);
                    //System.out.println("mean_equivalent = " + mean_equivalent);
                    double p = mean_equivalent / L[j];
                    //System.out.println("p = " + p);
                    Bp[j] = Binomial_Random_Generator(L[j], p);
                    //System.out.println("B[j] = " + Bp[j]);
                }

                for (int i = 0; i < REACTANTS; ++i) {
                    x_predicted[i] = x[i];

                    //System.out.println("x_predicted[i] = " + x_predicted[i]);

                    for (int j = 0; j < REACTIONS; ++j) {
                        x_predicted[i] += SMatrix[j][i] * Bp[j];
                    }
                    //System.out.println("delta = " + (x_predicted[i] - x[i]));
                }
                //System.out.println("bool = " + (xp[2][j] - xp[1][j]));

                //System.out.println("a = " + ap[0][j]);
                //System.out.println("delta = " + SMatrix[j][i] *  ((23 / 12 * ap[2][j])  - (16 / 12 * ap[1][j])  + (5 / 12 * ap[0][j])));
                //System.out.println("delta = " + SMatrix[j][i] *  ((23 / 12 * ap[2][j])  - (16 / 12 * ap[1][j])  + (5 / 12 * ap[0][j])));

                //System.out.println("tao = " + tao);
                //System.out.println("SMatrix[j][i] = " + SMatrix[j][i]);
                //System.out.println("delta1 = " + SMatrix[j][i] *  ((23 / 12 * ap[2][j]) /* - (16 / 12 * ap[1][j]) */ + (5 / 12 * ap[0][j])));
                //System.out.println("delta2 = " + SMatrix[j][i] *  (16 / 12 * ap[1][j]));


                double[] a_predicted = new double[REACTIONS];

                for (int j = 0; j < REACTIONS; ++j) {
                    a_predicted[j] = k[j];
                    for (int i = 0; i < REACTANTS; ++i) {
                        if (SMatrix[j][i] == -1) {
                            a_predicted[j] *= x_predicted[i];
                        }
                        if (SMatrix[j][i] == -2){
                            a_predicted[j] *= 0.5 * x_predicted[i] * (x_predicted[i] - 1);
                        }
                    }
                    //System.out.println("a_predicted[j] = " + a_predicted[j]);
                }

                for (int j = 0; j < REACTIONS; ++j) {
                    //System.out.println("product[j] = " + x_predicted[j] * tao);
                    double mean_equivalent = (9.00 / 24.00 * a_predicted[j] * tao + 19.00 / 24.00 * ap[3][j] * tao - 5.00 / 24.00 * ap[2][j] * tao + 1.00 / 24.00 * ap[1][j] * tao);
                    //System.out.println("mean_equivalent = " + mean_equivalent);
                    double p = mean_equivalent / L[j];
                    //System.out.println("p = " + p);
                    B[j] = Binomial_Random_Generator(L[j], p);
                    //System.out.println("B[j] = " + Bp[j]);
                }

                for (int i = 0; i < REACTANTS; ++i) {
                    for (int j = 0; j < REACTIONS; ++j)
                        x[i] += SMatrix[j][i] * B[j];
                    //System.out.println("x[i] = " + x[i]);
                }


                for (int i = 0; i < REACTIONS; i++) {
                    ap[0][i] = ap[1][i];
                    ap[1][i] = ap[2][i];
                    ap[2][i] = ap[3][i];
                }

                long endMulti = System.currentTimeMillis();
                long totalMulti = endMulti - startMulti;
                System.out.println("totalMulti = " + totalMulti);

            }


            for (int i = 0; i < REACTANTS; ++i) {
                writer.write(x[i] + ",");
                //System.out.println(x[i] + ",");

            }
            writer.write("\n");

            B_Max = 0;

            if (iteration < 4)
                for (int j = 0; j < REACTIONS; ++j) {
                    B_Max += Bp[j];
                }

            else
                for (int j = 0 ; j < REACTIONS; ++j)
                {
                    B_Max += B[j];
                }

            //System.out.println("B_Max = " + B_Max);
            steps += B_Max;

            t += tao;
            ++iteration;

        }

        long endTime = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        System.out.println("Total Time = " + totalTime);
        System.out.println("Iteration = " + iteration);
        System.out.println("Steps = " + steps);

        writer.flush();
        writer.close();

    }







}
