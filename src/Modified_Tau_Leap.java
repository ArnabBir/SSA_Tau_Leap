/**
 * Created by user on 6/5/2016.
 */
import java.io.FileWriter;
import java.io.IOException;

public class Modified_Tau_Leap
{


    private static int REACTIONS = 100;
    private static int REACTANTS = 101;
    private static int nc = 10;
    private static double EPSILON = 0.03;
    private static int INFINITY = 99999999;

    public static int minint(int[] a) {
        int min = a[0];

        for (int i = 1; i < a.length; ++i)
            if (a[i] < min)
                min = a[i];

        return min;
    }

    public static int maxindex(double [] a )
    {
        int max = 0;

        for (int i = 1; i < a.length; ++i)
        {
            if (a[i] > a[max])
            {
                max = i;
            }
        }

        return max;
    }

    public static double less_double(double a, double b)
    {
        if (a < b)
            return a;
        else
            return b;
    }

    public static double summation(double [] a )
    {
        int sum = 0;

        for (int i = 0; i < a.length ; ++i)
            sum += a[i];

        return sum;
    }

    public static int Poisson_Generator(double lambda)
    {
        double L = Math.exp(-lambda);
        System.out.println("Poi = " + L);
        double p = 1.0;
        int k = 0;

        do {
            k++;
            p *= Math.random();
        } while (p > L);

        return k - 1;
    }

    public static void main (String [] args) throws IOException
    {

        int[][] SMatrix = new int[REACTIONS][REACTANTS];
        //System.out.println(SMatrix[4][4]);
        for (int j = 0 ; j < REACTIONS; ++j)
            for (int i = 0 ; i < REACTANTS -1 ; ++i)
            {
                if ( i == j)
                {
                    SMatrix[j][i] = -1;
                    SMatrix[j][i + 1] = 1;
                }
            }
        //System.out.println(SMatrix[4][5]);

        double t = 0.0, tao = 0.0, tao_dash = 0.0, tao_2dash = 0.0;
        int iteration = 0;
        int [] x = new int [REACTANTS];
        double [] a = new double [REACTIONS];
        double [] k = new double [REACTIONS];

        //FileWriter writer = new FileWriter("Output_Sheet_2.csv",true);
        FileWriter writer = new FileWriter("Tau.csv");

        for (int i = 0; i < REACTANTS; ++i) {
            writer.write("Substrate " + (i + 1) + ",");
            x[i] = 100000;
        }

        writer.write("\n");

        for (int i = 0; i < REACTIONS; ++i)
        {
            //System.out.println("Enter the value of k" + (i + 1) + " :\n");
            //k[i] = sc.nextDouble();
            k[i] = 1;
        }

        long startTime = System.currentTimeMillis();

        while (x[0] > 5000)
        {
            try{
                for (int i = 0; i < REACTIONS ; ++i)
                {
                    a[i] = k[i];
                }

                int [] itr = new int [REACTIONS];

                for (int i = 0; i < REACTIONS; ++i)
                {
                    itr[i] = 0;

                    for (int j = 0; j < REACTANTS ; ++j)
                    {
                        if (SMatrix[i][j] < 0)
                        {
                            a[i] *= x[j];
                            itr[i] = itr[i] + 1;
                        }
                    }
                }



                double a0 = summation(a);

                int [][] reaction_array = new int [REACTIONS][];
                int [] L = new int [REACTIONS];

                // Calculating tao_dash

                for (int j = 0; j < REACTIONS; ++j)
                {
                    reaction_array[j] = new int [itr[j]];
                    int itr2 = 0;
                    for (int i = 0; i < REACTANTS && itr2 <= itr[j] ; ++i)
                    {
                        if ( SMatrix[j][i] < 0)
                        {
                            reaction_array[j][itr2] = (int) Math.floor(- x[i] / SMatrix[j][i]);
                            //writer.write ("reaction_arrray = " + reaction_array[j][itr2] + "\n");
                            ++itr2;
                        }

                    }
                }

                for (int j = 0; j < REACTIONS; ++j)
                {
                    L[j] = minint(reaction_array[j]);
                    //writer.write ("\nL[j] = " + L[j] + " ");
                }

                boolean [] is_critical = new boolean[REACTIONS];

                for (int j = 0; j < REACTIONS; ++j)
                {
                    if ( L[j] < nc)
                        is_critical[j] = true;
                    else
                        is_critical[j] = false;
                }

                boolean flag = false;

                for (int j = 0; j < REACTIONS; ++j)
                {
                    if (is_critical[j]  == false)
                        flag = true;
                }

                if(!flag)
                    tao_dash = INFINITY;
                else
                {
                    double [][] f = new double [REACTIONS][REACTIONS];

                    for (int j1 = 0 ; j1 < REACTIONS ; ++j1)
                    {
                        for (int j2 = 0; j2 < REACTIONS; ++j2)
                        {
                            f[j2][j1] = 0.0;

                            if (is_critical[j2] == false)
                            {
                                for (int i = 0; i < REACTANTS ; ++i)
                                {
                                    if (SMatrix[j1][i] < 0)
                                    {
                                        f[j1][j2] += (Math.abs(SMatrix[j1][i]) * a[j1] / x[i] ) * SMatrix[j2][i];
                                       // System.out.println("\n(Math.abs(SMatrix[j1][i]) * a[j1] / x[i] ) = " +(Math.abs(SMatrix[j1][i]) * a[j1] / x[i] ));

                                       // System.out.println("\nf[j1][j2] = " + f[j1][j2]);
                                        //System.out.println("\nSMatrix[j1][i] = " + SMatrix[j1][i] + "\n");
                                    }
                                }

                            }
                        }

                    }

                    double [] mean = new double [REACTIONS];
                    double [] variance = new double [REACTIONS];

                    for (int j1 = 0 ; j1 < REACTIONS ; ++j1)
                    {
                        mean[j1] = 0.0;
                        variance[j1] = 0.0;

                        for (int j2 = 0; j2 < REACTIONS; ++j2)
                        {
                            if (is_critical[j2] == false)
                            {
                                mean[j1] += f[j1][j2] * a[j2];
                                variance[j1] += f[j1][j2] * f[j1][j2] * a[j2];
                            }
                        }

                        //System.out.println("mean[j1] = " + mean[j1]);
                        //System.out.println(" variance[j1] = " +  variance[j1]);

                    }

                    int max_mean_index = maxindex(mean);
                    int max_variance_index = maxindex(variance);

                    //System.out.println("max_variance_index = " + max_variance_index);
                    //System.out.println("Math.abs(mean[max_mean_index] = " + Math.pow((EPSILON * a[max_variance_index]), 2) / variance[max_variance_index]  );

                    tao_dash = less_double( EPSILON * a0 / Math.abs(mean[max_mean_index]), Math.pow((EPSILON * a0), 2) / variance[max_variance_index]  );

                }
                //Calculating tao_2dash

                double ac0 = 0.0;

                for (int i = 0; i < REACTIONS; ++i)
                {
                    if (is_critical[i] == true)
                        ac0 += a[i];

                }

                tao_2dash = 1 / ac0 * Math.log(1 / Math.random());

                //Calculating tao

                tao = less_double(tao_dash, tao_2dash);
                System.out.println("tao = " + tao);

                int [] P = new int [REACTIONS];

                if(tao == tao_dash)
                {
                    for (int j = 0; j < REACTIONS; ++j)
                    {
                        if (is_critical[j] == true)
                        {
                            P[j] = 0;
                        }

                        else
                        {
                            //writer.write ("\nLambda = " + tao * a[j] + "\n");
                            P[j] = Poisson_Generator(a[j] * tao);
                            //System.out.println("P = " + P[j]);
                            //System.out.println("Mean = " + (a[j] * tao));
                        }
                    }
                }

                if (tao == tao_2dash)
                {
                    double r = Math.random();
                    double r_ac0 = ac0 * r;
                    double sum_ac = 0.0;
                    int index = 0;;

                    while ( r_ac0 > sum_ac)
                    {
                        if (is_critical[index] == true)
                        {
                            sum_ac += a[index];
                        }
                        ++index;
                    }

                    //writer.write ("\nTao_dash = " + tao_dash + "\n");
                    //writer.write ("\nTao_2dash = " + tao_2dash + "\n");
                    //writer.write ("\nTao = " + tao + "\n");

                    for (int j = 0; j < REACTIONS; ++j)
                    {

                        if (is_critical[j] == true)
                        {
                            P[j] = 0;
                        }

                        else
                        {
                            //System.out.println("\nLambda = " + tao * a[j]);
                            P[j] = Poisson_Generator(a[j] * tao);
                        }
                    }

                    P[index - 1] = 1;
                }

                //writer.write("\n" + P[0] + " " + P[1] + " " + P[2] +  " " + P[3] + " " + P[4] + "\n");

                for (int i = 0; i < REACTANTS; ++i)
                {
                    for (int j = 0; j < REACTIONS; ++j)
                    {
                 //       System.out.println("\nLambda = " + tao * a[j]);
                        x[i] += P[j] * SMatrix[j][i];
                    }
                }

                //writer.write ("\nTao_dash = " + tao_dash + "\n");
                //writer.write ("\nTao_2dash = " + tao_2dash + "\n");
                System.out.println("\nTao = " + tao + "\n");
                t += tao;

               // writer.write("\nAfter iteration : " + (iteration) + "\n");

                //writer.write("\nAfter iteration : " + (iteration) + "\nA(t" + (iteration) + ") = " + (x[0]) + " B(t" + (iteration) + ") = " + (x[1]) + " C(t" + (iteration) +") = " + (x[2]) +  " D(t" + (iteration) + ") = " + (x[3]) + " E(t" + (iteration) + ") = " + (x[4]) + " F(t" + (iteration) + ") = " + (x[5]) + " G(t" + (iteration) + ") = " + (x[6]));
                ++iteration;
            }
            catch(Exception e)
            {
                break;
            }
        }

        for (int i = 0; i < REACTANTS; ++i)
        {
            writer.write(x[i] + ",");
        }
        writer.write("\n");



        long endTime   = System.currentTimeMillis();
        long totalTime = endTime - startTime;
        System.out.println(totalTime);
        System.out.println(iteration);


        writer.flush();
        writer.close();
    }

}
