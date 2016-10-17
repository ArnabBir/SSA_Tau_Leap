/**
 * Created by user on 6/6/2016.
 */
public class Epsilon {

    public static void main (String [] args)
    {
        double e = 1;
        for (int i = 1000; i >0; --i)
        {
            e = 1 + 10000000 * e / i;
            System.out.println(e);

        }
        System.out.println(e);
    }
}
