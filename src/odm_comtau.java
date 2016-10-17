import java.lang.*;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.io.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.Scanner;

public class odm_comtau {

        Scanner sc = new Scanner(System.in);
        int M = sc.nextInt(), N = sc.nextInt();
        //static int M=10, N = 9;
        int[][] a = new int [M][N];
        double[] k = new double [N];  int[] conc = new int[M];  int[] no_spec_per_reac = new int[N], no_input_spec_per_reac = new int[N], no_depend_per_reac = new int[N];
        int [][] spec_per_reac = new int [N][], input_spec_per_reac = new int [N][], change_input_conc = new int[N][], change_conc = new int[N][], depend_per_reac = new int[N][];
        double[] prop = new double[N]; double totalprop = 0;double transitionfiretime, totaltransitionfiringtime=0; long search = 0;
        //public static int[] exec(int selrxn, int no_spec_per_reac[], int spec_per_reac[][], int change_conc[][], int conc[])
        double prop_c(int mol, int n_mol)
        {
            double prop_value1=1;
            if(mol < n_mol)
                return 0;

            for(int i=0; i < n_mol; i++)
                prop_value1 = prop_value1*(mol-i);
            prop_value1=prop_value1/n_mol;
            return prop_value1;

        }
        public int firingtran()
        {
            double sam = totalprop*Math.random(), selector=0;
            //transitionfiretime=-Math.log1p(Math.random())/totalprop;
            //System.out.print("sam ="+sam);
            int i;
            for(i=0; i<N; i++)
            {
                selector+= prop[i];
                search++;
                if(selector >= sam)
                    break;
            }

            return i;
        }
        public  void exec(int selrxn)
        {
            for(int i=0; i< no_spec_per_reac[selrxn]; i++)
            {
                int x = spec_per_reac[selrxn][i];
                //System.out.println(x);
                conc[x]+= change_conc[selrxn][i];
                if (conc[x] < 0)
                    conc[x] = 0;
            }
            //System.out.println(conc[0]);
            // return conc;

        }
        public  void calcprop()
        {
            for(int i=0; i< N; i++)
                prop[i]=k[i];
            for (int i=0; i< N; i++)
            {
                for(int j=0;j < no_spec_per_reac[i]; j++)
                {
                    if (change_conc[i][j] < 0)
                        //  prop[i]*=k[i]*conc[spec_per_reac[i][j]];
                        //prop[i]*=k[i]*Math.pow(conc[spec_per_reac[i][j]], -change_conc[i][j]);
                        prop[i]*=prop_c(conc[spec_per_reac[i][j]], -change_conc[i][j]);
                }

            }



            for (int i=0; i< N; i++)
                totalprop+=prop[i];
            // System.out.print("totalprop ="+totalprop);
        }
        public  void dependencygraph()
        {


            int flage = 0;
            for(int i=0;i<N;i++)
            {

                for(int j=0;j<N;j++)
                {
                    if (i == j)
                        continue;
                    flage = 0;
                    for(int i1=0;i1<no_spec_per_reac[i];i1++)
                    {
                        for(int j1=0; j1 < no_input_spec_per_reac[j];j1++)
                        {
                            if(spec_per_reac[i][i1] == input_spec_per_reac[j][j1])
                            {
                                no_depend_per_reac[i]++;
                                flage = 1;
                                break;
                            }
                        }
                        if ( flage == 1)
                            break;
                    }
                }
            }
            for(int i=0; i<N; i++)
            {
                depend_per_reac[i] = new int[no_depend_per_reac[i]];
                no_depend_per_reac[i]=0;
            }
            for(int i=0;i<N;i++)
            {

                for(int j=0;j<N;j++)
                {
                    if (i == j)
                        continue;
                    flage = 0;
                    for(int i1=0;i1<no_spec_per_reac[i];i1++)
                    {
                        for(int j1=0; j1 < no_input_spec_per_reac[j];j1++)
                        {
                            if(spec_per_reac[i][i1] == input_spec_per_reac[j][j1])
                            {
                                depend_per_reac[i][no_depend_per_reac[i]++] = j;
                                flage = 1;
                                break;
                            }
                        }
                        if ( flage == 1)
                            break;
                    }
                }
            }
/*  for(int i=0;i<N;i++)
  {

   for(int j=0;j<no_depend_per_reac[i];j++)
   {
    System.out.print(" "+depend_per_reac[i][j]);
   }
   System.out.println();
  }
    */



        }

        public  void updatepropesity(int selrxn)
        {
            int x;
            totalprop = totalprop - prop[selrxn];
            for(int i=0;i<no_depend_per_reac[selrxn];i++)
                totalprop = totalprop - prop[depend_per_reac[selrxn][i]];

            //prop[selrxn]= k[selrxn]*conc[input_spec_per_reac[selrxn][0]];
            //prop[selrxn]= k[selrxn]*Math.pow(conc[input_spec_per_reac[selrxn][0]], - change_input_conc[selrxn][0]);
            prop[selrxn]=k[selrxn]*prop_c(conc[input_spec_per_reac[selrxn][0]], - change_input_conc[selrxn][0]);
            for(int j=1;j< no_input_spec_per_reac[selrxn]; j++)
            {
                //prop[selrxn] = prop[selrxn]*conc[input_spec_per_reac[selrxn][j]];
                //prop[selrxn] = prop[selrxn]*Math.pow(conc[input_spec_per_reac[selrxn][j]], - change_input_conc[selrxn][j]);
                prop[selrxn]=prop[selrxn]*prop_c(conc[input_spec_per_reac[selrxn][j]], - change_input_conc[selrxn][j]);
            }
            totalprop = totalprop + prop[selrxn];
            for(int i = 0; i< no_depend_per_reac[selrxn];i++)
            { x = depend_per_reac[selrxn][i];
                //prop[x] = k[x]*conc[input_spec_per_reac[x][0]];
                //prop[x] = k[x]*Math.pow(conc[input_spec_per_reac[x][0]], - change_input_conc[x][0]);
                prop[x] = k[x]*prop_c(conc[input_spec_per_reac[x][0]], - change_input_conc[x][0]);
                for(int j=1;j< no_input_spec_per_reac[x]; j++)
                {
                    //prop[x] = prop[x]*conc[input_spec_per_reac[x][j]];
                    //prop[x] = prop[x]*Math.pow(conc[input_spec_per_reac[x][j]], - change_input_conc[x][j]);
                    prop[x] = prop[x]*prop_c(conc[input_spec_per_reac[x][j]], - change_input_conc[x][j]);
                }
                // outf << " trans " << i << " prop " << prop[i]<< "total prop "<<totalprop;
                totalprop = totalprop + prop[x];
            }
        }

        public static void main(String args[])
        {
            System.out.println(" Enter the values of the number of species and reactions ");
            odm_comtau obj = new odm_comtau();
            System.out.println("Reading File from Java code");
            //Name of the file
            String fileName="xyz.txt", fileName1="k.txt", fileName2="conc.txt";
            try
            {
                //Create object of FileReader
                FileReader inputFile = new FileReader(fileName), inputFile1 = new FileReader(fileName1), inputFile2 = new FileReader(fileName2);
                //Instantiate the BufferedReader Class
                BufferedReader bufferReader = new BufferedReader(inputFile), bufferReader1 = new BufferedReader(inputFile1), bufferReader2 = new BufferedReader(inputFile2);
                //Variable to hold the one line data
                String line;
                String[] result;

                File file = new File("filenameoutput.csv");
                // if file doesnt exists, then create it
                if (!file.exists())
                {
                    file.createNewFile();
                }
                FileWriter fw = new FileWriter(file.getAbsoluteFile());
                BufferedWriter bw = new BufferedWriter(fw);
                // Read file line by line and print on the console
                int m11=0, m12=0;
                while ((line = bufferReader.readLine()) != null)
                {

                    //System.out.println(line);
                    //result=line.split("\\s");
                    result=line.split(",");

                    for (int i=0; i<obj.N;i++)
                    //double value = Double.parseDouble(result[i]);
                    {
                        int value = Integer.parseInt(result[i]);
                        obj.a[m11][m12++]=value;
                        //System.out.print("  "+value);
                    }
                    m11++;
                    m12=0;
                    //System.out.println();
                    //	bw.write(result[0]);
                    //	bw.write("\n");


                }
                int count_k=0;
                while((line = bufferReader1.readLine()) != null)
                {
                    //System.out.println(line);
                    obj.k[count_k++] = Double.parseDouble(line);
                }
                //for(int i=0; i<obj.N; i++)
                //System.out.println(obj.k[i]);
                int[] conc1 = new int[obj.M];
                int count_c=0;
                while((line = bufferReader2.readLine()) != null)
                {
                    //System.out.println(line);
                    conc1[count_c++] = Integer.parseInt(line.trim());
                }

                for(int j=0;j<obj.N;j++)
                {
                    for(int i=0;i<obj.M;i++)
                    {
                        if (obj.a[i][j]!= 0)
                            obj.no_spec_per_reac[j]++;  // number of species connected per reaction
                        if (obj.a[i][j]<0)
                            obj.no_input_spec_per_reac[j]++; // number of input species connected per reaction
                    }

                }
                for(int i=0; i<obj.N;i++)
                {
                    obj.spec_per_reac[i] = new int[obj.no_spec_per_reac[i]];  // species per reaction
                    obj.change_conc[i] = new int[obj.no_spec_per_reac[i]];    // stoichiometry of the species per reaction
                    obj.input_spec_per_reac[i] = new int[obj.no_input_spec_per_reac[i]]; // input species per reaction
                    obj.change_input_conc[i] = new int[obj.no_input_spec_per_reac[i]]; // stoichiometry of the input species per reaction

                }

                for (int j=0; j< obj.N; j++)
                {	int x=0, x1=0;
                    for ( int i=0; i< obj.M; i++)
                    {
                        if (obj.a[i][j]!= 0)
                        {
                            obj.spec_per_reac[j][x] = i;
                            obj.change_conc[j][x++] = obj.a[i][j];
                        }
                        if (obj.a[i][j] < 0)
                        {
                            obj.input_spec_per_reac[j][x1] = i;
                            obj.change_input_conc[j][x1++] = obj.a[i][j];
                        }
                    }
                }


                long avg_simu_time=0; long searchtime=0, updatetime=0;
                System.out.println(" How many simulations you want ?");
                Scanner sc1 = new Scanner(System.in);
                int simu = sc1.nextInt();
                for(int iter=0; iter <simu;iter++)
                {
                    for (int i=0; i< obj.M; i++)
                    {
                        obj.conc[i]=conc1[i];
                    }

                    obj.dependencygraph();
                    obj.calcprop();





                    //conc = exec(selrxn, no_spec_per_reac, spec_per_reac, change_conc, conc);  // execute the reaction

                    int currenttime=0, endtime=10000000;
                    double curtao=0, endtao=0.000001;
                    long startTime = System.currentTimeMillis();
                    long iteration = 0;

                    //while(currenttime < endtime)
                    while(curtao < endtao)
                    {
                        //calcpropensity(t, p);
                        long startsearchTime = System.currentTimeMillis();
                        int  selrxn = obj.firingtran();
                        if (iteration < 10)
                            System.out.println(obj.totalprop);
                        ++iteration;
                        long stopsearchTime = System.currentTimeMillis();
                        searchtime = searchtime + stopsearchTime - startsearchTime;
                        // if (selrxn == 0 )
                        //outf << selrxn << " ";
                        //sparetran++;
                        obj.transitionfiretime=-Math.log(Math.random())/obj.totalprop;
                        //System.out.println("tau ="+ obj.transitionfiretime);
                        curtao = curtao + obj.transitionfiretime;
                        obj.exec(selrxn);
                        long startupdateTime = System.currentTimeMillis();
                        obj.updatepropesity(selrxn);
                        long stopupdateTime = System.currentTimeMillis();
                        updatetime = updatetime + stopupdateTime - startupdateTime;
                        //currenttime = currenttime + 1;
                        obj.totaltransitionfiringtime+=obj.transitionfiretime;
                    }

                    //System.out.println();
                    //System.out.println(obj.conc[obj.M-1]);
//	for (int i=0; i< M; i++)
//   System.out.print(" "+conc[i]);

                    //for (int i=0; i< obj.M; i++)
                        //System.out.print(" "+obj.conc[i]);
                    long stopTime = System.currentTimeMillis();
                    long elapsedTime = stopTime - startTime;
                    System.out.println(iter+" th = " +elapsedTime+ " ms");
                    avg_simu_time=avg_simu_time + elapsedTime;
                    obj.totalprop=0;
                }
                System.out.println(" Search time "+searchtime+"ms");
                System.out.println(" Update time "+updatetime+"ms");
                System.out.println(" average simulation time "+avg_simu_time/simu+" ms ");
                System.out.println("Total number of searching ="+obj.search);
                bw.close();
                //Close the buffer reader
                bufferReader.close();
                bufferReader1.close();
                bufferReader2.close();
            }catch(Exception e)
            {
                System.out.println("Error while reading file line by line:"+ e.getMessage());
            }
        }



        public  void abc(int xyx)
        {
            System.out.println("Hello World");
        }
    }

