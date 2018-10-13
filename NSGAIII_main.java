package jmetal.metaheuristics.nsgaIII;

import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import java.io.*;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.thetadea.ThetaDEA;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.Fonseca;
import jmetal.problems.Kursawe;
import jmetal.problems.ProblemFactory;
import jmetal.problems.Schaffer;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
import jmetal.problems.M2M.IMB1;
import jmetal.problems.M2M.IMB10;
import jmetal.problems.M2M.IMB2;
import jmetal.problems.M2M.IMB3;
import jmetal.problems.M2M.IMB4;
import jmetal.problems.M2M.IMB5;
import jmetal.problems.M2M.IMB6;
import jmetal.problems.M2M.IMB7;
import jmetal.problems.M2M.IMB8;
import jmetal.problems.M2M.IMB9;
import jmetal.problems.M2M.MOP1;
import jmetal.problems.M2M.MOP2;
import jmetal.problems.M2M.MOP3;
import jmetal.problems.M2M.MOP4;
import jmetal.problems.M2M.MOP5;
import jmetal.problems.M2M.MOP6;
import jmetal.problems.M2M.MOP7;
import jmetal.problems.MaF.MaF1;
import jmetal.problems.MaF.MaF2;
import jmetal.problems.MaF.MaF3;
import jmetal.problems.MaF.MaF4;
import jmetal.problems.MaF.MaF5;
import jmetal.problems.MaF.MaF6;
import jmetal.problems.MaF.MaF7;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG5;
import jmetal.problems.WFG.WFG6;
import jmetal.problems.WFG.WFG7;
import jmetal.problems.WFG.WFG8;
import jmetal.problems.WFG.WFG9;
import jmetal.problems.ZDT.ZDT1;
import jmetal.problems.ZDT.ZDT2;
import jmetal.problems.ZDT.ZDT3;
import jmetal.problems.ZDT.ZDT4;
import jmetal.problems.ZDT.ZDT6;
import jmetal.problems.cec2009Competition.UF1;
import jmetal.problems.cec2009Competition.UF10;
import jmetal.problems.cec2009Competition.UF2;
import jmetal.problems.cec2009Competition.UF3;
import jmetal.problems.cec2009Competition.UF4;
import jmetal.problems.cec2009Competition.UF5;
import jmetal.problems.cec2009Competition.UF6;
import jmetal.problems.cec2009Competition.UF7;
import jmetal.problems.cec2009Competition.UF8;
import jmetal.problems.cec2009Competition.UF9;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculator2;
import jmetal.util.Configuration;
import jmetal.util.JMException;

public class NSGAIII_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	/**
	 * @param args
	 *            Command line arguments.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException
	 *             Usage: three options -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main problemName -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main problemName
	 *             paretoFrontFile
	 */
	public static void printGD(String path, double[] GD) {
		try {
			/* Open the file */
			FileOutputStream fos = new FileOutputStream(path);// java文件输出流，创建文件流
			OutputStreamWriter osw = new OutputStreamWriter(fos);// OutputStreamWriter是字符流通向字节流的桥梁
			BufferedWriter bw = new BufferedWriter(osw);// 缓冲区
			for (int i = 0; i < GD.length; i++) {
				bw.write(GD[i] + " ");// 写到缓冲区
				bw.newLine(); // 换行
			}

			/* Close the file */
			bw.close();
		} catch (IOException e) {
			Configuration.logger_.severe("Error acceding to the file");
			e.printStackTrace();
		}
	} // printGD

	public static void main(String args[]) throws JMException,
			ClassNotFoundException, SecurityException, IOException {
		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("NSGAIII_main.log");
		logger_.addHandler(fileHandler_);

		int m = 4;

		int wfgk = 2 * (m - 1);
		int l = 20;

		int mafk = 20;
		int mafn = m + mafk - 1;

		int dtlzk = 5;
		int dtlzn = m + dtlzk - 1;

		// 120 4D 7 0
		// 165 4D 8 0
		// 100 5D 4 3
		// 252 6D 5 0
		// 210 7D 4 0
		// 120 8D 3 0
		// 156 8D 3 2
		// 330 8D 4 0
		// 220 10D 3 0
		// 275 10D 3 2

		int div1 = 8;
		int div2 = 0;

		int popSize = 165;
		int maxGEN = 500;

		int initi = 9;
		int termi = initi + 1;

		for (int fun = initi; fun < termi; fun++) {
			int runtimes = 30;
			// double[] IGDarray=new double[runtimes];
			// double[] HVarray = new double[runtimes];
			double[] GDarray = new double[runtimes];
			double[] IGDarray = new double[runtimes];
			double[] spreadarray = new double[runtimes];
			double[] Hypervolume = new double[runtimes];
			long Execution_time = 0;

			Problem problem = null; // The problem to solve
			Algorithm algorithm; // The algorithm to use
			Operator crossover; // Crossover operator
			Operator mutation; // Mutation operator
			Operator selection; // Selection operator

			HashMap parameters; // Operator parameters

			Solution referencePoint;

			QualityIndicator indicators;// Object to get quality indicators
			indicators = null;

			if (args.length == 1) {
				Object[] params = { "Real" };
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} // if
			else if (args.length == 2) {
				Object[] params = { "Real" };
				problem = (new ProblemFactory()).getProblem(args[0], params);
				indicators = new QualityIndicator(problem, args[1]);
			} // if
			else { // Default problem

				if (fun == 1) {
					problem = new ZDT1("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT1_501.txt");
				} // problem = new WFG1("Real");
				if (fun == 2) {
					problem = new ZDT2("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT2_501.txt");
				} // problem = new WFG1("Real");
				if (fun == 3) {
					problem = new ZDT3("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT3_269.txt");
				}
				if (fun == 4) {
					problem = new ZDT4("Real", 10);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT4_501.txt");
				} // problem = new WFG1("Real");
				if (fun == 5) {
					problem = new ZDT6("Real", 10);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT6_774.txt");
				} // problem = new WFG1("Real");
				if (fun == 6) {
					problem = new DTLZ1("Real", dtlzn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ1mm_"
									+ m + ".txt");
				}
				if (fun == 7) {
					problem = new DTLZ2("Real", dtlzn, m);
					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ2mm_"
									+ m + ".txt");
				} // problem = new WFG1("Real");
				if (fun == 8) {
					problem = new DTLZ3("Real", dtlzn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ3mm_"
									+ m + ".txt");
				} // problem = new WFG1("Real");
				if (fun == 9) {
					problem = new DTLZ4("Real", dtlzn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ4mm_"
									+ m + ".txt");
				}
				if (fun == 10) {
					problem = new DTLZ5("Real", dtlzn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ5mm_"
									+ m + ".txt");
				} // problem = new WFG1("Real");
				if (fun == 11) {
					problem = new DTLZ6("Real", dtlzn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ6mm_"
									+ m + ".txt");
				} // problem = new WFG1("Real");
				if (fun == 12) {
					problem = new DTLZ7("Real", dtlzn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ7mm_"
									+ m + ".txt");
				}
				if (fun == 13) {
					problem = new WFG1("Real", wfgk, l, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG1mm_"
									+ m + ".txt");
				} // problem = new WFG1("Real");
				if (fun == 14) {
					problem = new WFG2("Real", wfgk, l, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG2mm_"
									+ m + ".txt");
				} // problem = new WFG1("Real");
				if (fun == 15) {
					problem = new WFG3("Real", wfgk, l, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG3mm_"
									+ m + ".txt");
				}
				if (fun == 16) {
					problem = new WFG4("Real", wfgk, l, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG4mm_"
									+ m + ".txt");
				} // problem = new WFG1("Real");
				if (fun == 17) {
					problem = new WFG5("Real", wfgk, l, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG5mm_"
									+ m + ".txt");
				}
				if (fun == 18) {
					problem = new WFG6("Real", wfgk, l, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG6mm_"
									+ m + ".txt");
				} // problem = new WFG1("Real");
				if (fun == 19) {
					problem = new WFG7("Real", wfgk, l, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG7mm_"
									+ m + ".txt");
				} // problem = new WFG1("Real");
				if (fun == 20) {
					problem = new WFG8("Real", wfgk, l, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG8mm_"
									+ m + ".txt");
				}
				if (fun == 21) {
					problem = new WFG9("Real", wfgk, l, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG9mm_"
									+ m + ".txt");
				} // problem = new WFG1("Real");
				if (fun == 22) {
					problem = new Fonseca("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\Fonseca.pf");
				} // problem = new WFG1("Real");
				if (fun == 23) {
					problem = new Kursawe("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\Kursawe.pf");
				}
				if (fun == 24) {
					problem = new Schaffer("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\Schaffer.pf");
				} // problem = new WFG1("Real");
				if (fun == 25) {
					problem = new UF1("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\CEC2009_UF1.dat");
				}
				if (fun == 26) {
					problem = new UF2("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\CEC2009_UF2.dat");
				}
				if (fun == 27) {
					problem = new UF3("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\CEC2009_UF3.dat");
				}
				if (fun == 28) {
					problem = new UF4("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\CEC2009_UF4.dat");
				}
				if (fun == 29) {
					problem = new UF5("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\CEC2009_UF5.dat");
				}
				if (fun == 30) {
					problem = new UF6("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\CEC2009_UF6.dat");
				}
				if (fun == 31) {
					problem = new UF7("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\CEC2009_UF7.dat");
				}
				if (fun == 32) {
					problem = new MOP1("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MOP1.dat");
				}
				if (fun == 33) {
					problem = new MOP2("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MOP2.dat");
				}
				if (fun == 34) {
					problem = new MOP3("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MOP3.dat");
				}
				if (fun == 35) {
					problem = new MOP4("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MOP4.dat");
				}
				if (fun == 36) {
					problem = new MOP5("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MOP5.dat");
				}
				if (fun == 37) {
					problem = new UF8("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\CEC2009_UF8.DAT");
				}
				if (fun == 38) {
					problem = new UF9("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\CEC2009_UF9.DAT");
				}
				if (fun == 39) {
					problem = new UF10("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\CEC2009_UF10.DAT");
				}

				if (fun == 40) {
					problem = new MOP6("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MOP6.dat");
				}
				if (fun == 41) {
					problem = new MOP7("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MOP7.dat");
				}

				if (fun == 42) {
					problem = new IMB1("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\IMB1.dat");
				}

				if (fun == 43) {
					problem = new IMB2("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\IMB2.dat");
				}

				if (fun == 44) {
					problem = new IMB3("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\IMB3.dat");
				}

				if (fun == 45) {
					problem = new IMB7("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\IMB7.dat");
				}

				if (fun == 46) {
					problem = new IMB8("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\IMB8.dat");
				}

				if (fun == 47) {
					problem = new IMB9("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\IMB9.dat");
				}

				if (fun == 48) {
					problem = new IMB4("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\IMB4.dat");
				}

				if (fun == 49) {
					problem = new IMB5("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\IMB5.dat");
				}

				if (fun == 50) {
					problem = new IMB6("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\IMB6.dat");
				}

				if (fun == 51) {
					problem = new IMB10("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\IMB10.dat");
				}

				if (fun == 52) {
					problem = new MaF1("Real", mafn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MaF1mm_"
									+ m + ".txt");
				}

				if (fun == 53) {
					problem = new MaF2("Real", mafn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MaF2mm_"
									+ m + ".txt");
				}

				if (fun == 54) {
					problem = new MaF3("Real", mafn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MaF3mm_"
									+ m + ".txt");
				}

				if (fun == 55) {
					problem = new MaF4("Real", mafn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MaF4mm_"
									+ m + ".txt");
				}

				if (fun == 56) {
					problem = new MaF5("Real", mafn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MaF5mm_"
									+ m + ".txt");
				}

				if (fun == 57) {
					problem = new MaF6("Real", mafn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MaF6mm_"
									+ m + ".txt");
				}

				if (fun == 58) {
					problem = new MaF7("Real", mafn, m);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\MaF7mm_"
									+ m + ".txt");
				}

				// System.out.print(problem.getName());
			} // else

			/*
			 * referencePoint = new Solution(problem.getNumberOfObjectives());
			 * for (int i = 0; i < referencePoint.numberOfObjectives(); i++) {
			 * referencePoint.setObjective(i, 2*i+10); }
			 */

			algorithm = new NSGAIII_SBX(problem);
			// SBX version is an original version

			algorithm.setInputParameter("populationSize", popSize);
			// algorithm.setInputParameter("maxEvaluations", 300000);
			algorithm.setInputParameter("normalize", true);

			algorithm.setInputParameter("div1", div1);
			algorithm.setInputParameter("div2", div2);

			algorithm.setInputParameter("FEsRecord", 5000);
			algorithm.setInputParameter("maxGenerations", maxGEN);

			// Mutation and Crossover for Real codification
			parameters = new HashMap();
			parameters.put("probability", 1.0);
			parameters.put("distributionIndex", 30.0);
			crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",
					parameters);// original
			// crossover =
			// CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover",
			// parameters);

			parameters = new HashMap();
			parameters.put("probability", 1.0 / problem.getNumberOfVariables());
			parameters.put("distributionIndex", 20.0);
			mutation = MutationFactory.getMutationOperator(
					"PolynomialMutation", parameters);

			parameters = null;
			selection = SelectionFactory.getSelectionOperator(
					"RandomSelection", parameters);

			// Add the operators to the algorithm
			algorithm.addOperator("crossover", crossover);
			algorithm.addOperator("mutation", mutation);
			algorithm.addOperator("selection", selection);

			for (int i = 0; i < runtimes; i++) {
				long initTime = System.currentTimeMillis();

				System.out.print(" " + (i + 1) + "_");

				SolutionSet population = algorithm.execute();

				Execution_time += (System.currentTimeMillis() - initTime);

				// IGDarray[i]=indicators.getIGD(population);
				// GDarray[i]=indicators.getGD(population);
				IGDarray[i] = indicators.getCECIGD(population);
				// spreadarray[i]=indicators.getSpread(population);
				// Hypervolume[i] = indicators.getHypervolume(population);
				Hypervolume[i] = indicators.getSuperHV(population);
				// System.out.println("IGD: "+IGDarray[i]);
				// System.out.println("HV: "+Hypervolume[i]);

				if (i == 0) {
					// printGD("IGD" + "-NSGAIII-" + problem.getName() + ".txt",
					// IGDarray[i],false);
					printGD("HV" + "-NSGAIII-" + problem.getName() + ".txt",
							Hypervolume[i], false);
				} else {

					// printGD("IGD" + "-NSGAIII-" + problem.getName() + ".txt",
					// IGDarray[i],true);
					printGD("HV" + "-NSGAIII-" + problem.getName() + ".txt",
							Hypervolume[i], true);
				}
				// wfghvCalculator2 Indicators2 = new
				// wfghvCalculator2(population,
				// fun);
				// HVarray[i] = Indicators2.calculatewfghv();
				// population.printObjectivesToFile("Run" + i +
				// problem.getName()
				// + m + "-NSGAIII");
				// population.printVariablesToFile("Variables" + i
				// + problem.getName() + "NSGAⅢ");

				// System.out.println(indicators.getHypervolume(population));
				// System.out.println(indicators.getSuperHV(population));

			}

			// printGD("IGD"+ "-NSGAIII-"+problem.getName(),IGDarray);
			// printGD("HV"+ "-NSGAIII-"+problem.getName(),Hypervolume);
			double sumIGD = 0;
			double sumHV = 0.0;
			double stdIGD = 0;
			double stdHV = 0;
			for (int i = 0; i < runtimes; i++) {
				sumIGD += IGDarray[i];
				// sumHV+=HVarray[i];
				sumHV += Hypervolume[i];
			}

			// logger_.info("Total execution time: " + Execution_time + "ms");
			System.out.println();
			// System.out.println("avrIGD-fun" + fun + "= " + sumIGD /
			// runtimes);
			// System.out.println("avrHV-fun" + fun + "= " + sumHV / runtimes);
			System.out.println(sumHV / runtimes);
			for (int i = 0; i < runtimes; i++) {
				// stdIGD += Math.pow(IGDarray[i] - sumIGD / runtimes, 2);
				stdHV += Math.pow(Hypervolume[i] - sumHV / runtimes, 2);
			}
			// System.out.println("stdIGD-fun" + fun + "= "
			// + Math.sqrt(stdIGD / runtimes));
			System.out.println(Math.sqrt(stdHV / runtimes));
			// System.out.println("stdHV-fun" + fun + "= "
			// + Math.sqrt(stdHV / runtimes));
			System.out.println("NSGA3_" + m + "D");

		} // for-fun

	}// main

	public static void printGD(String path, double GD, boolean appendWrite) {
		try {
			/* Open the file */
			// FileOutputStream fos = new FileOutputStream(path, true);//
			// java文件输出流，创建文件流
			FileOutputStream fos = new FileOutputStream(path, appendWrite);// java文件输出流，创建文件流
			OutputStreamWriter osw = new OutputStreamWriter(fos);// OutputStreamWriter是字符流通向字节流的桥梁
			BufferedWriter bw = new BufferedWriter(osw);// 缓冲区
			bw.write(GD + " ");// 写到缓冲区
			/* Close the file */
			bw.close();
		} catch (IOException e) {
			Configuration.logger_.severe("Error acceding to the file");
			e.printStackTrace();
		}
	} // printGDs
}
