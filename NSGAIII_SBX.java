package jmetal.metaheuristics.nsgaIII;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

import jmetal.core.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.*;

import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.vector.OneLevelWeightVectorGenerator;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;

public class NSGAIII_SBX extends Algorithm {

	private int populationSize_;

	private int div1_;
	private int div2_;

	private SolutionSet population_;
	SolutionSet offspringPopulation_;
	SolutionSet union_;

	int generations_;

	Operator crossover_;
	Operator mutation_;
	Operator selection_;

	double[][] lambda_; // reference points

	boolean normalize_; // do normalization or not

	public NSGAIII_SBX(Problem problem) {

		super(problem);

	} // NSGAII

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		int maxGenerations_;
		double IGDarray = 10;
		double Hypervolume = 0;
		generations_ = 0;
		int fes = 0;
		maxGenerations_ = ((Integer) this.getInputParameter("maxGenerations"))
				.intValue();

		int FEsRecord_ = ((Integer) this.getInputParameter("FEsRecord"))
				.intValue();

		div1_ = ((Integer) this.getInputParameter("div1")).intValue();

		div2_ = ((Integer) this.getInputParameter("div2")).intValue();

		normalize_ = ((Boolean) this.getInputParameter("normalize"))
				.booleanValue();

		VectorGenerator vg = new TwoLevelWeightVectorGenerator(div1_, div2_,
				problem_.getNumberOfObjectives());
		lambda_ = vg.getVectors();

		// System.out.println(vg.getVectors().length);

		// populationSize_ = vg.getVectors().length;
		populationSize_ = ((Integer) this.getInputParameter("populationSize"))
				.intValue();
		if (populationSize_ % 2 != 0)
			populationSize_ += 1;

		String testName = getProblem().getName();

		// QualityIndicator indicators = new QualityIndicator(getProblem(),
		// "E:\\new_multiobjective\\jMetal\\Pareto_front\\" + testName +
		// ".dat");

		System.out.print(getProblem().getName());

		mutation_ = operators_.get("mutation");
		crossover_ = operators_.get("crossover");
		selection_ = operators_.get("selection");

		initPopulation();
		fes += populationSize_;

		while (generations_ < maxGenerations_) {

			// System.out.println(fes);

			offspringPopulation_ = new SolutionSet(populationSize_);
			Solution[] parents1 = new Solution[2];
			Solution[] parents2 = new Solution[2];
			Solution[] parents = new Solution[3];
			for (int i = 0; i < (populationSize_ / 2); i++) {
				if (generations_ < maxGenerations_) {
					// obtain parents

					parents = (Solution[]) selection_.execute(population_);
					Solution[] offSpring = (Solution[]) crossover_
							.execute(parents);

					// parents1 = (Solution[]) selection_.execute(population_);
					// parents2 = (Solution[]) selection_.execute(population_);
					// parents[0] = parents1[0];
					// parents[1] = parents1[1];
					// parents[2] = parents2[0];
					// Solution offSpring1 = (Solution) crossover_.execute(new
					// Object[] { parents2[0], parents });
					// Solution offSpring2 = (Solution) crossover_.execute(new
					// Object[] { parents2[0], parents });
					// Solution[] offSpring=new Solution[2];
					// offSpring[0]=offSpring1;
					// offSpring[1]=offSpring2;

					mutation_.execute(offSpring[0]);
					mutation_.execute(offSpring[1]);

					problem_.evaluate(offSpring[0]);
					problem_.evaluateConstraints(offSpring[0]);
					fes++;

					// if (fes % FEsRecord_ == 0 && fes != maxGenerations_ *
					// populationSize_) {
					// // System.out.println(fes);
					// IGDarray = indicators.getCECIGD(population_);
					// Hypervolume = indicators.getHypervolume(population_);
					// printGD("IGD" + "-NSGAIII-" + testName + "-process.txt",
					// IGDarray);
					// printGD("HV" + "-NSGAIII-" + testName + "-process.txt",
					// Hypervolume);
					// }

					problem_.evaluate(offSpring[1]);
					problem_.evaluateConstraints(offSpring[1]);
					fes++;

					// if (fes % FEsRecord_ == 0 && fes != maxGenerations_ *
					// populationSize_) {
					// // System.out.println(fes);
					// IGDarray = indicators.getCECIGD(population_);
					// Hypervolume = indicators.getHypervolume(population_);
					// printGD("IGD" + "-NSGAIII-" + testName + "-process.txt",
					// IGDarray);
					// printGD("HV" + "-NSGAIII-" + testName + "-process.txt",
					// Hypervolume);
					// }

					offspringPopulation_.add(offSpring[0]);
					offspringPopulation_.add(offSpring[1]);

				} // if
			} // for

			union_ = ((SolutionSet) population_).union(offspringPopulation_);

			// Ranking the union
			Ranking ranking = new NondominatedRanking(union_);

			int remain = populationSize_;
			int index = 0;
			SolutionSet front = null;
			population_.clear();

			// Obtain the next front
			front = ranking.getSubfront(index);

			while ((remain > 0) && (remain >= front.size())) {

				for (int k = 0; k < front.size(); k++) {
					population_.add(front.get(k));
				} // for

				// Decrement remain
				remain = remain - front.size();

				// Obtain the next front
				index++;
				if (remain > 0) {
					front = ranking.getSubfront(index);
				} // if
			}

			if (remain > 0) { // front contains individuals to insert

				new Niching(population_, front, lambda_, remain, normalize_)
						.execute();
				remain = 0;
			}

			// if(getProblem().getNumberOfObjectives()==2){
			// if(fes%FEsRecord_==0)
			// {
			// IGDarray=indicators.getCECIGD(population_);
			// Hypervolume=indicators.getHypervolume(population_);
			// printGD("IGD"+ "-NSGAIII-"+testName+".txt",IGDarray);
			// printGD("HV"+ "-NSGAIII-"+testName+".txt",Hypervolume);
			// }
			// }else{
			// if((generations_+1)%500==0)
			// {
			// IGDarray=indicators.getCECIGD(population_);
			// Hypervolume=indicators.getHypervolume(population_);
			// printGD("IGD"+
			// "-NSGAIII-"+testName+"_"+(generations_+1)+".txt",IGDarray);
			// printGD("HV"+
			// "-NSGAIII-"+testName+"_"+(generations_+1)+".txt",Hypervolume);
			// }
			// }

			generations_++;

		}

		// IGDarray = indicators.getCECIGD(population_);
		// Hypervolume = indicators.getHypervolume(population_);
		// printGD("IGD" + "-NSGAIII-" + testName + "-process.txt", IGDarray);
		// printGD("HV" + "-NSGAIII-" + testName + "-process.txt", Hypervolume);
		// printGD("IGD" + "-NSGAIII-" + testName + "-process.txt", "\n");
		// printGD("HV" + "-NSGAIII-" + testName + "-process.txt", "\n");

		Ranking ranking = new NondominatedRanking(population_);
		return ranking.getSubfront(0);

	}

	public void initPopulation() throws JMException, ClassNotFoundException {

		population_ = new SolutionSet(populationSize_);

		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problem_);

			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);

			population_.add(newSolution);
		} // for
	} // initPopulation

	public static void printGD(String path, double GD) {
		try {
			/* Open the file */
			FileOutputStream fos = new FileOutputStream(path, true);// java文件输出流，创建文件流
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

	public static void printGD(String path, String GD) {
		try {
			/* Open the file */
			FileOutputStream fos = new FileOutputStream(path, true);// java文件输出流，创建文件流
			OutputStreamWriter osw = new OutputStreamWriter(fos);// OutputStreamWriter是字符流通向字节流的桥梁
			BufferedWriter bw = new BufferedWriter(osw);// 缓冲区
			bw.write(GD);// 写到缓冲区
			/* Close the file */
			bw.close();
		} catch (IOException e) {
			Configuration.logger_.severe("Error acceding to the file");
			e.printStackTrace();
		}
	} // printGDs
} // NSGA-III
