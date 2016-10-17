/**
 * Created by user on 6/24/2016.
 */
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

class odm_newhigherorder {
    // this is optimized direct method it will support higher order reactions. this is relevent to make jar file taking input from keyboard
// input files xyz.txt for incidence matrix
// k.txt for the k values
// conc.txt for the populations

    static class odm_higherorder1
    {
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
                conc[x]+= change_conc[selrxn][i];
                if (conc[x] < 0)
                    conc[x] = 0;

            }
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
            odm_higherorder1 obj = new odm_higherorder1();
            int xyx = 5;
            obj.abc(xyx); // double[] k = new double [N]; int[] conc = new int[M]; int[] no_spec_per_reac = new int[N];
            //int [][] spec_per_reac = new int [N][], change_conc = new int[N][]; double[] prop = new double[N]; double totalprop = 0;
/*		a[0][0]=-1;
  for(int i=1;i<M-1;i++)
  {
    a[i][i-1]=1;
    a[i][i]=-1;   //linear chain model
  }
  a[M-1][N-1]=1;   */
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

		/*	System.out.println(" a[][] ");
			for(int i=0; i<M; i++)
			{
			for(int j=0; j<N; j++)
			  System.out.print(" "+a[i][j]);
			  System.out.println();
			}		*/

/*
  int n,m, t=0,p,q;
	for(n=1; n<M/2; n++)
	{
	  for(m=n; m< M-n; m++)
	  {
	    if(n==m)
	    {
	      a[n][t]=-2;
	      a[n+m][t++]=1;		// colloidal aggregation model
	    }
	    else
	    {
	      a[n][t]=-1;
	      a[m][t]=-1;
	      a[n+m][t++]=1;
	    }
	  }
	}
	for(p=1;p<M;p++)
	{
	  for(q=1;q<=p/2;q++)
	  {
	    a[p][t]=-1;
	    a[q][t]=1;
	    a[p-q][t++]=1;
	  }
	}
  */

  /*
  a[0][0]=-1;
  for(int i=1;i<M;i++)
  {
    a[i][i-1]=1;	// cyclic chain model
    a[i][i]=-1;
  }

  a[0][N-1] = 1; */

                //  a[2][5]=-1; a[7][5]=-1;


   /*
    for(int i=0;i<M;i++)
    {
	  for(int j=0;j<N;j++)
	    System.out.print(" "+ a[i][j]);
	    System.out.println();

    }*/

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

  /* for(int i=0; i<N;i++)
    {
    for (int j=0; j< no_spec_per_reac[i]; j++)
      System.out.print(" "+ spec_per_reac[i][j]);
      System.out.println();
    }
    for(int i=0; i<N;i++)
    {
    for (int j=0; j< no_spec_per_reac[i]; j++)
      System.out.print(" "+ change_conc[i][j]);
      System.out.println();
    }

    for(int i=0; i<N;i++)
    {
    for (int j=0; j< no_input_spec_per_reac[i]; j++)
      System.out.print(" "+ input_spec_per_reac[i][j]);
      System.out.println();
    }
    for(int i=0; i<N;i++)
    {
    for (int j=0; j< no_input_spec_per_reac[i]; j++)
      System.out.print(" "+ change_input_conc[i][j]);
      System.out.println();
    } */



//    for (int i=0; i< obj.N; i++)
                //     obj.k[i]=1.0;					// propensity value
                //   k[N-2]=1000.0;
/*
  int timedur = 100;
	conc[0] =105663*timedur;
	conc[1] =14000*timedur;
	conc[2] =394*timedur;
	conc[3] =00153*timedur;
	conc[4] =97*timedur;
	conc[5] =1492*timedur;
	conc[6] =00061*timedur;
	conc[7] =9992*timedur;
	conc[8] =00653*timedur;
	conc[9] =00005*timedur;
	conc[10]=00001*timedur;
	conc[11]=20601*timedur;
	conc[12]=658*timedur;
	conc[13]=84*timedur;
	conc[14]=109*timedur;
	conc[15]=84*timedur;
	conc[16]=2*timedur;
	conc[17]=4*timedur;
	conc[18]=0*timedur;
	conc[19] =0*timedur;
	conc[20] =25*timedur;
	conc[21] =4*timedur;
	conc[22] =31136*timedur;
	conc[23] =47*timedur;
	conc[24] =14*timedur;
	conc[25] =127*timedur;
	conc[26] =154*timedur;
	conc[27] =63*timedur;
	conc[28] =0*timedur;
	conc[29]=5*timedur;
	conc[30]=0*timedur;
	conc[31]=24*timedur;
	conc[32]=0*timedur;
	conc[33]=24*timedur;
	conc[34]=0*timedur;
	conc[35]=25*timedur;
	conc[36]=25*timedur;
	conc[37]=8*timedur;
	conc[38]=0;
	conc[39]=0;

  double mul=1.0;


	k[0]=mul* 1.0*15.8*conc[0]*(- conc[2]*conc[10]/3900.0 + conc[1] + 33.2*conc[1]*conc[37]/(1.14*15.8))/(1.44*(0.1 + conc[0])*(1.0 + conc[37]/1.03 + conc[1]*(1.0 + conc[37]/1.14)/1.44 + (1.55 + conc[2]/0.0045)*(1.0 + conc[37]/1.03) + (conc[11] + conc[29])/2.7 + conc[37]*(conc[11] + conc[29])/(3.44*1.03)));
	k[1]=mul* 935.0*(- conc[3]/0.3925 + conc[2])/(0.182*(1.0 + conc[3]/0.071) + conc[2]);
	k[2]=mul* 239.0*(- conc[4]*conc[10]/100000.0 + conc[3]*conc[1])/((0.1 + conc[3])*(0.068 + conc[1])*(1.0 + 0.001072*Math.pow((1.0 + conc[36]/0.01),4)*Math.pow((1.0 + conc[37]/0.44),4)/(Math.pow((1.0 + conc[3]/0.1),4)*Math.pow((1.0 + (conc[18] + conc[28])/0.033),4))));
	k[3]=mul* 98.91*(conc[4] - conc[5]*conc[6]/0.114)/(0.0071*(1.0 + conc[4]/0.0071 + conc[6]/0.0572 + conc[4]*conc[6]/(0.0071*0.176) + conc[5]*(0.1906 + conc[6])/(0.0364*0.0572)));
	k[4]=mul* 5456.6*(conc[5] - conc[6]/0.0407)/(conc[5] + 0.838*(1.0 + conc[6]/0.428)) ;
	k[5]=mul* 4300.0*(-conc[9]*conc[16]/0.000192 + conc[6]*conc[8]*conc[7])/(0.005*0.05*3.9*(-1.0 + (1.0 + conc[9]/0.0035)*(1.0 + conc[16]/0.0083) + (1.0 + conc[6]/0.005)*(1.0 + conc[8]/0.05)*(1.0 + conc[7]/3.9)));
	k[6]=mul* 5000.0*(conc[9]*conc[10] - conc[12]*conc[1]/1455.0)/(0.002*0.35*(-1.0 + (1.0 + conc[9]/0.002)*(1.0 + conc[10]/0.35) + (1.0 + conc[12]/1.2)*(1.0 + conc[1]/0.48)));
	k[7]=mul* 76000.0*(conc[9] - (conc[11] + conc[29])/100000.0)/(1.0 + (conc[11] + conc[29])/0.04);
	k[8]=mul* 0.53*(conc[11] - conc[12]/100000.0 + conc[29])/(0.2 + conc[11] + conc[29]);
	k[9]=mul*2000.0*(conc[12] - conc[13]/0.145)/(5.0*(1.0 + conc[13]/1.0) + conc[12]);
	k[10]=mul*1500.0*(conc[13] - conc[14]/1.7)/(conc[13] + 1.0*(1.0 + conc[14]/1.0));
	k[11]=mul*570.0*(conc[10]*conc[14] - conc[1]*conc[15]/13790)/((0.474 + conc[10])*(0.225 + conc[14])*(1.0 + 19.0*Math.pow((1.0 + (conc[36] + conc[1])/3.39),4)/(Math.pow((1.0 + conc[4]/0.005),4)*Math.pow((1.0 + conc[14]/0.225),4))));
	k[12]=mul*2.8e6*(conc[16]*conc[15] - conc[34]*conc[8]/9090.0);
	k[13]=mul*243.4*(conc[17]*conc[15] - conc[34]*conc[19]/14181.8);
	k[14]=mul*1.68*conc[1];
	k[15]=mul*1380.0*(- conc[35]*conc[10]/0.25 + conc[18]*conc[1])/(0.08*0.09*(conc[35]*conc[10]/Math.pow(0.11,2) + (conc[35] + conc[10])/0.11 + (1.0 + conc[18]/0.08)*(1.0 + conc[1]/0.09)));
	k[16]=mul*162.0*(conc[2]*conc[19] - conc[20]*conc[17]/2000.0)/(0.0667*0.00367*(1.0 + (conc[36] + conc[1])/0.749 + (conc[11] + conc[29])/2.289 + (1.0 + conc[2]/0.0667)*conc[19]/0.00367 + conc[17]/0.00312));
	k[17]=mul*1575.0*(conc[20]*conc[19] - conc[17]*conc[23]/141.7)/(0.01*0.018*((conc[36] + conc[1])/0.154 + (1.0 + conc[20]/0.01 + (conc[11] + conc[29])/0.12)*(1.0 + conc[19]/0.018) + ((1.0 + conc[20]/0.058)*conc[17]/0.0045)));
	k[18]=mul*90.0*(-Math.pow(conc[22],2)*conc[19]/(1.04*Math.pow(20.0,2)*0.07) + conc[21]*conc[17]/(0.0652*0.00852))/(1.0 + ((1.0 + conc[22]*(1.0 + conc[22]/20.0)/20.0)*conc[19]/0.07) + ((1.0 + conc[21]/0.0652)*conc[17]/0.00852));
	k[19]=mul*0.03*conc[22];
	k[20]=mul*4634.0*(conc[23] - conc[25]/2.7)/(conc[23] + 0.19*(1.0 + conc[25]/0.5));
	k[21]=mul*730.0*(- conc[24]/3.0 + conc[23])/(0.78*(1.0 + conc[24]/2.2) + conc[23]);
	k[22]=mul*23.5*(- conc[6]*conc[26]/1.05 + conc[24]*conc[25])/(0.00496*conc[26] + conc[6]*(12.432 + 0.41139*conc[26]) + conc[24]*(0.3055 + 0.00774*conc[26]) + 48.8*conc[6]*conc[25] + (0.4177 + conc[24])*conc[25]);
	k[23]=mul*27.2*((conc[6]*conc[26] - conc[27]*conc[3]/1.05))/(0.006095*conc[3] + conc[27]*(0.1733 + 0.8683*conc[3]) + (0.04765 + 0.4653*conc[3])*conc[6] + 2.524*conc[27]*conc[26] + (0.00823 + conc[6])*conc[26]);
	k[24]=mul*1.1*(conc[1]*conc[24] - 1.0*conc[28]/100000.0)/((0.03 + conc[1]) + (0.57 + conc[24]));
	k[25]=mul*23.5*((- conc[3]*conc[6]/1.2 + conc[27]*conc[25]))/(0.0003*conc[3] + conc[27]*(0.3055 + 0.122*conc[3]) + (0.0548 + 0.0287*conc[3])*conc[6] + (0.00184 + conc[27])*conc[25] + 0.215*conc[6]*conc[25]);
	k[26]=mul*100.0*(1.0 - conc[7]/1.0);
	k[27]=mul*10000.0*(1.68 - conc[34]/1.0);
	k[28]=mul*10000.0*(0.084 - conc[15]/1.0);
	k[29]=mul*1.0e7*(conc[1] - conc[36]*conc[37]/0.072);
	k[30]=mul*1.0e7*(conc[10] - conc[35]*conc[37]/0.76);
	k[31]=mul*1.0e7*(conc[28] - conc[18]*conc[37]/16.64);
	k[32]=mul*1.0e7*(conc[29] - conc[11]*conc[37]/1.667);
	k[33]=mul*1.0e7*(conc[30] - conc[19]*conc[38]/0.0002);
	k[34]=mul*1.0e7*(conc[31] - conc[17]*conc[38]/0.00001);
	k[35]=mul*1.0e7*(conc[32] - conc[19]*conc[39]/0.00001);
	k[36]=mul*1.0e7*(conc[33] - conc[17]*conc[39]/0.0002);
	double[] k1 = new double [37];

	for(int i=0; i<37; i++)
	{
	  k1[i]=k[i];

	}
	for(int i=0; i<37; i++)
	{
	  k[i]=k1[i];
	  k[37+i]=k1[i];
	}
	  int reac2=0;
	for(int j=0; j<20; j++)
	{
	  for(int i=0; i<37; i++)
	  {
	    k[reac2+i]=k1[i];
	  //  System.out.print(" "+k[reac2+i]);
	  }
	  reac2=reac2+37;
	//  System.out.println();
	}

	*/



                // species concentration
                long avg_simu_time=0; long searchtime=0, updatetime=0;
                System.out.println(" How many simulations you want ?");
                Scanner sc1 = new Scanner(System.in);
                int simu = sc1.nextInt();
                for(int iter=0; iter <simu;iter++)
                {
//for (int i=0; i< obj.M; i++)
                    // obj.conc[i]=10000000;
                    for (int i=0; i< obj.M; i++)
                    {
                        obj.conc[i]=conc1[i];
                        //System.out.println(obj.conc[i]);
                    }

                    // for (int i=0; i< M; i++)
                    //  System.out.print(" "+conc[i]);
                    //  System.out.println();
                    obj.dependencygraph();
                    obj.calcprop();





                    //conc = exec(selrxn, no_spec_per_reac, spec_per_reac, change_conc, conc);  // execute the reaction

                    int currenttime=0, endtime=10000000;

                    long startTime = System.currentTimeMillis();

                    while(currenttime < endtime)
                    {
                        //calcpropensity(t, p);
                        long startsearchTime = System.currentTimeMillis();
                        int  selrxn = obj.firingtran();
                        long stopsearchTime = System.currentTimeMillis();
                        searchtime = searchtime + stopsearchTime - startsearchTime;
                        // if (selrxn == 0 )
                        //outf << selrxn << " ";
                        //sparetran++;
                        obj.transitionfiretime=-Math.log1p(Math.random())/obj.totalprop;
                        obj.exec(selrxn);
                        long startupdateTime = System.currentTimeMillis();
                        obj.updatepropesity(selrxn);
                        long stopupdateTime = System.currentTimeMillis();
                        updatetime = updatetime + stopupdateTime - startupdateTime;
                        currenttime = currenttime + 1;
                        obj.totaltransitionfiringtime+=obj.transitionfiretime;
                    }

                    System.out.println();
                    System.out.println(obj.conc[obj.M-1]);
//	for (int i=0; i< M; i++)
//   System.out.print(" "+conc[i]);

                    // for (int i=0; i< M; i++)
                    // System.out.print(" "+conc[i]);
                    //   Calendar cal = Calendar.getInstance();
                    //   SimpleDateFormat sdf = new SimpleDateFormat("HH:mm:ss");
                    //     System.out.println( sdf.format(cal.getTime()) );
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
}
