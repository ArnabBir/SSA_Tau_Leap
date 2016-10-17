import java.math.BigDecimal;

/**
 * Created by user on 6/6/2016.
 */
public class Poisson_Generator {
    public static long Poisson_Generator(double lambda)
    {
        BigDecimal L = new BigDecimal (Math.exp(-lambda));
        //System.out.println("L = " + L);
        BigDecimal p = new BigDecimal(1.0000);
        long k = 0;

        do {
            k++;
            //System.out.println(p);
            BigDecimal temp = new BigDecimal ( Math.random());
            p = p.multiply(temp);
        } while (p.compareTo(L) == 1);

        return k;
    }

    public static void main (String[] args)
    {

        long start_time = System.currentTimeMillis();
//        for(int i = 0 ; i < 100000; ++i)
            Poisson_Generator(2E10);
            //System.out.println(Poisson_Generator(250));
        long end_time = System.currentTimeMillis();
        System.out.println("Total time = " + (end_time - start_time));

        //{
        //    System.out.println(Poisson_Generator(250));
        //}
    }
}
