import java.io.FileWriter;
import java.io.IOException;


/**
 * Created by user on 6/10/2016.
 */
public class Colloidal_Aggregation_Matrix {

    private static int REACTANTS = 100;
    private static int DECOMPOSITION_REACTIONS = (int) Math.floor(Math.pow(REACTANTS, 2) / 4);
    private static int REACTIONS = 2 * DECOMPOSITION_REACTIONS;

    public static void main(String[] args) throws IOException {

        int[][] SMatrix = new int[REACTIONS][REACTANTS];
        int reaction = 0;

        while (reaction < DECOMPOSITION_REACTIONS) {
            for (int n = 1; n <= Math.floor(REACTANTS / 2); ++n) {
                for (int m = n; m <= REACTANTS - n; ++m) {
                    --SMatrix[reaction][n - 1];
                    --SMatrix[reaction][m - 1];
                    ++SMatrix[reaction][n + m - 1];
                    ++reaction;
                }
            }
        }

        System.out.println("Reaction = " + reaction);

        while (reaction < REACTIONS) {
            for (int n = 1; n <= Math.floor(REACTANTS / 2); ++n) {
                for (int m = n; m <= REACTANTS - n; ++m) {
                    ++SMatrix[reaction][n - 1];
                    ++SMatrix[reaction][m - 1];
                    --SMatrix[reaction][n + m - 1];
                    ++reaction;
                }
            }
        }

        FileWriter fw = new FileWriter("Colloidal_Aggregation_Matrix.csv");

        for (int j = 0; j < REACTANTS; ++j) {
            for (int i = 0; i < REACTIONS; ++i) {
                fw.write(SMatrix[i][j] + ",");
            }
            fw.write("\n");
        }

        fw.flush();
        fw.close();
    }
}