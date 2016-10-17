/**
 * Created by user on 6/6/2016.
 */
public class Binomial_random_generator {

    public static int Binomial_Random_Generator(int n, double p) {
        int X = 0;

        for (int i = 0; i < n; ++i) {
            if (p <= Math.random())
                ++X;
        }
        return X;
    }
    public static void main (String[] args)
    {
        long start_time = System.currentTimeMillis();
        for(int i = 0 ; i < 100000; ++i)
        {
            Binomial_Random_Generator(500, 0.5);
        //System.out.println(Binomial_Random_Generator(500, 0.5));
        }
        long end_time = System.currentTimeMillis();
        System.out.println("Total time = " + (end_time - start_time));

    }

}
