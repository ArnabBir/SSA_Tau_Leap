import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class Colloidial_Aggregation_Direct_Method {
	
	
	
	private static int REACTANTS = 100;
	private static int DECOMPOSITION_REACTIONS = (int) Math.floor(Math.pow(REACTANTS, 2) / 4);
	private static int REACTIONS = 2 * DECOMPOSITION_REACTIONS;
	private static double SIMULATION_TIME = 0.01;
	private static int nc = 10;
	private static double EPSILON = 0.3;
	private static int INFINITY = 99999999;

	public static int minint(int[] a) {
		int min = a[0];

		for (int i = 1; i < a.length; ++i)
			if (a[i] < min)
				min = a[i];

		return min;
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

	public static int[] Multinomial_Random_Generator(int K_sum, double[] theta, int m) {
		int[] K = new int[m];
		double theta_sum = summation(theta);

		K[0] = Binomial_Random_Generator(K_sum, theta[0]);

		for (int i = 1; i < m; ++i) {
			K_sum -= K[i - 1];
			theta_sum -= theta[i - 1];
			K[i] = Binomial_Random_Generator(K_sum, theta[i] / theta_sum);
		}

		return K;

	}

	public static int K_Generator(int[] x, int[][] SMatrix, int[] itr, boolean[] is_critical, boolean bool) {

		double[] max_x_array = new double[REACTANTS];
		int[] max_SMatrix = new int[REACTANTS];
		double[] K_guess = new double[REACTANTS];
		double epsilon_x;

		for (int i = 0; i < REACTANTS; ++i) {

			max_SMatrix[i] = Math.abs(SMatrix[0][i]);

			for (int j = 0; j < REACTIONS; ++j) {
				if (is_critical[j] == bool) {

					/*epsilon_x = EPSILON / itr[j] * x[i];
					System.out.println("epsilon_x = " + String.valueOf(epsilon_x));
					if (epsilon_x > 1) {
						//System.out.println("Greater than 1");
						max_x_array[i] = epsilon_x;
					} else {
						max_x_array[i] = 1;
					}*/
					if (Math.abs(SMatrix[j][i]) != 0)
						if (Math.abs(SMatrix[j][i]) > max_SMatrix[i])
							max_SMatrix[i] = Math.abs(SMatrix[j][i]);

				}
			}
			
			epsilon_x = EPSILON * x[i];
			
			//System.out.println("epsilon_x = " + String.valueOf(epsilon_x));
			if (epsilon_x > 1) {
				//System.out.println("Greater than 1");
				max_x_array[i] = epsilon_x;
			} else {
				max_x_array[i] = 1;
			}
			
				//System.out.println("max_x_array = " + String.valueOf(max_x_array[i]));
				//System.out.println("max_SMatrix = " + String.valueOf(max_SMatrix[i]));
				K_guess[i] = max_x_array[i] / max_SMatrix[i];
				//System.out.println("K_guess = " + String.valueOf(K_guess[i]));

		}

		double K = minimum(K_guess);
		//System.out.println("K_guess = " + String.valueOf(K));

		if (K > 1)
			return (int) K;
		else
			return 1;
	}

	public static long factorial(int n) {
		long fact_value = 1;
		if (n == 0 || n == 1) {
			return 1;
		} else
			for (int i = 2; i <= n; ++i) {
				fact_value *= i;
			}
		return fact_value;

	}

	public static double truncated_poisson(double a, double tao, int Km, int K) {
		double den_sum = 0.0;

		for (int j = 0; j < K; ++j) {
			// System.out.println("Power = " + Math.pow(a * tao, j) );
			// System.out.println("Fact = " + factorial(j));
			// System.out.println("Sum = " + Math.pow(a * tao, j) /
			// factorial(j));
			den_sum += Math.pow(a * tao, j) / factorial(j);
			// System.out.println(den_sum);

		}

		// System.out.println("Den = " + den_sum);
		// System.out.println("Power = " + Math.pow(a * tao, Km));
		double tp = Math.pow(a * tao, Km) / (factorial(Km) * den_sum);
		// System.out.println(tp);

		return tp;
	}

	public static int truncated_poisson_generator(double a, double tao, int K) {

		double r = Math.random();
		double cumulative_tp_sum = 0.0;
		int Km = 0;

		while (r > cumulative_tp_sum) {
			cumulative_tp_sum += truncated_poisson(a, tao, Km, K);
			++Km;
		}

		return Km;

	}

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
				
		File file = new File("SMatrix_Collodial_aggregation_K_Leap.csv");
        System.out.println(file.createNewFile());
        FileWriter writer = new FileWriter(file);
        
        double t = 0.0, tao = 0.0, tao_dash = 0.0, tao_2dash = 0.0;
		int iteration = 1;
		int[] x = new int[REACTANTS];
		double[] a = new double[REACTIONS];
		double[] k = new double[REACTIONS];
		
		for (int i = 0; i < REACTANTS; ++i) {
			writer.write("Substrate " + (i + 1) + ",");
			x[i] = 100000;
		}

		writer.write("\n");

		for (int i = 0; i < REACTIONS; ++i) {
			k[i] = 1;
		}
		
		while (iteration <= 500) {
			try

			{
				for (int j = 0; j < REACTIONS; ++j) {
					a[j] = k[j];
				}

				int[] itr = new int[REACTIONS];

				for (int i = 0; i < REACTIONS; ++i) {
					itr[i] = 0;

					for (int j = 0; j < REACTANTS; ++j) {
						
						if (SMatrix[i][j] == -1) {
							a[i] *= x[j];
							itr[i] = itr[i] + 1;
						}
						
						if (SMatrix[i][j] == -2) {
							//System.out.println("hello!");
							a[i] *= 0.5 * x[j] * ( x[j] - 1);
							itr[i] = itr[i] + 1;
						}
							
					}
				}

				int[][] reaction_array = new int[REACTIONS][];
				int[] L = new int[REACTIONS];

				// Calculating tao_dash

				for (int j = 0; j < REACTIONS; ++j) {
					
					reaction_array[j] = new int[itr[j]];
					int itr2 = 0;
					
					for (int i = 0; i < REACTANTS && itr2 <= itr[j]; ++i) {
						if (SMatrix[j][i] < 0) {
							reaction_array[j][itr2] = (int) Math.floor(-x[i] / SMatrix[j][i]);
							// System.out.println(reaction_array[j][itr2]);
							++itr2;
						}

					}
					// System.out.println(reaction_array[j].length);

				}

				for (int j = 0; j < REACTIONS; ++j) {
					L[j] = minint(reaction_array[j]);
					// System.out.println("L[j] = " + L[j] + " ");
				}

				double a0 = summation(a);

				boolean[] is_critical = new boolean[REACTIONS];
				int critical = 0, not_critical = 0;
				double a0_nc = 0.0, a0_c = 0.0;

				for (int j = 0; j < REACTIONS; ++j) {
					if (L[j] < nc) {
						is_critical[j] = true;
						a0_c += a[j];
						++critical;
					}

					else {
						is_critical[j] = false;
						a0_nc += a[j];
						++not_critical;
					}
				}

				int Kc = K_Generator(x, SMatrix, itr, is_critical, true);
				int Knc = K_Generator(x, SMatrix, itr, is_critical, false);
				//System.out.println(Kc);
				//System.out.println(Knc);

				tao_dash = (double) (Knc) / a0_nc * Math.log(1 / Math.random());
				tao_2dash = (double) (Kc) / a0_c * Math.log(1 / Math.random());

				int[] Km_nc = new int[not_critical];
				int[] Km_c = new int[critical];

				// System.out.println(tao_dash);
				// System.out.println(tao_2dash);
				// System.out.println(not_critical);

				if (tao_dash < tao_2dash) {

					tao = tao_dash;
					double[] theta = new double[not_critical];
					int iteration_nc = 0;
					for (int j = 0; j < REACTIONS/*
													 * && iteration_nc <
													 * not_critical
													 */; ++j) {
						if (is_critical[j] == false) {
							theta[iteration_nc++] = a[j] / a0_nc;
						}
					}

					Km_nc = Multinomial_Random_Generator(Knc, theta, not_critical);

					if (not_critical < REACTIONS) {

						//System.out.println("Bye bye world!!");
						//System.out.println(a0_c + "  " + tao + " " + Kc + " " + Knc);
						int K_critical = truncated_poisson_generator(a0_c, tao, Kc);
						// System.out.println(K_critical);
						double[] theta2 = new double[critical];
						int iteration_c = 0;
						for (int j = 0; j < REACTIONS /*
														 * && iteration_c <
														 * critical
														 */; ++j) {
							if (is_critical[j] == true) {
								theta2[iteration_c++] = a[j] / a0_c;
							}
						}

						Km_c = Multinomial_Random_Generator(K_critical, theta2, critical);
					}
				}

				if (tao_dash > tao_2dash) {

					tao = tao_2dash;
					double[] theta = new double[critical];
					int iteration_c = 0;
					for (int j = 0; j < REACTIONS /*
													 * && iteration_c < critical
													 */; ++j) {
						if (is_critical[j] == true) {
							theta[iteration_c++] = a[j] / a0_c;
						}
					}

					Km_c = Multinomial_Random_Generator(Kc, theta, critical);

					if (critical != REACTIONS) {
						//System.out.println("Hello world!!");
						int K_not_critical = truncated_poisson_generator(a0_nc, tao, Knc);

						double[] theta2 = new double[not_critical];
						int iteration_nc = 0;
						for (int j = 0; j < REACTIONS /*
														 * && iteration_nc <
														 * not_critical
														 */; ++j) {
							if (is_critical[j] == false) {
								theta2[iteration_nc++] = a[j] / a0_nc;
							}
						}

						Km_nc = Multinomial_Random_Generator(K_not_critical, theta2, not_critical);
					}
				}

				if (tao_dash == tao_2dash) {
					tao = tao_dash;

					double[] theta = new double[not_critical];
					int iteration_nc = 0;
					for (int j = 0; j < REACTIONS && iteration_nc < not_critical; ++j) {
						if (is_critical[j] == false) {
							theta[iteration_nc++] = a[j] / a0_nc;
						}
					}

					Km_nc = Multinomial_Random_Generator(Knc, theta, not_critical);

					double[] theta2 = new double[critical];
					int iteration_c = 0;
					for (int j = 0; j < REACTIONS && iteration_c < critical; ++j) {
						if (is_critical[j] == true) {
							theta2[iteration_c++] = a[j] / a0_c;
						}
					}

					Km_c = Multinomial_Random_Generator(Kc, theta2, critical);
				}

				int[] K = new int[REACTIONS];
				int c_itr = 0;
				int nc_itr = 0;

				for (int j = 0; j < REACTIONS && c_itr <= critical && nc_itr <= not_critical; ++j) {
					// System.out.println("Hello");
					if (is_critical[j] == false) {
						K[j] = Km_nc[nc_itr];
						++nc_itr;
					}

					else {
						K[j] = Km_c[c_itr];
						++c_itr;
					}
				}

				t += tao;

				for (int i = 0; i < REACTANTS; ++i) {
					for (int j = 0; j < REACTIONS; ++j) {
						// System.out.println("K = " + K[j]);
						// System.out.println("K = " + Km_nc[j]);
						x[i] += K[j] * SMatrix[j][i];
					}
				}

				for (int i = 0; i < REACTANTS; ++i) {
					writer.write(x[i] + ",");
				}
				writer.write("\n");

				++iteration;

			} catch (Exception e) {
				break;
			}
		}
		
		writer.flush();
		writer.close();
        
        
	}
}