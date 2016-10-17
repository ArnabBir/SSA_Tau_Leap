/**
 * Created by user on 6/10/2016.
 */

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.util.Scanner;

public class Matrix_Reader {

    public static void main(String[] args) throws IOException {

        Scanner sc = new Scanner(System.in);

        System.out.println("Enter the number of reactants: ");
        int REACTANTS = sc.nextInt();
        System.out.println("Enter the number of reactions: ");
        int REACTIONS = sc.nextInt();
//        System.out.println("Enter the directory of Stiochiometry Matrix : ");
//        String xy = sc.nextLine();

        FileReader fr = new FileReader("xyz.txt");
        BufferedReader br = new BufferedReader(fr);
        String str;
        FileWriter fw = new FileWriter("Strongly_Coupled_3_Matrix.csv");
        int[][] SMatrix = new int[REACTIONS][REACTANTS];
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

            System.out.println(row_index);
            ++row_index;
        }


        for(int j = 0; j < REACTANTS; ++j){
            for (int i = 0; i < REACTIONS; ++i){
                fw.write(SMatrix[i][j] + ",");
            }
            fw.write("\n");
        }

        fw.flush();
        fw.close();

        // br.flush();
        br.close();


    }

}
