/**
 * MOEADACD.java
 * 
 * MOEADACD.java
 * 
 * This is main implementation of MOEA/D-ACD.
 * 
 * Author:
 * 		Luping Wang <zjlxwlp8412@126.com>
 * 
 * Reference:
 * 		Luping Wang, Qingfu Zhang, Aimin Zhou, Maoguo Gong and Licheng Jiao 
 * 		"Constrained Subproblems in Decomposition based Multiobjective Evolutionary Algorithm"
 * 		IEEE Transactions on Evolutionary Computation.
 * 
 * 
 * Copyright (c) 2015 Luping Wang
* 
 * Note: This is a free software developed based on the open source project 
 * jMetal<http://jmetal.sourceforge.net>. The copy right of jMetal belongs to 
 * its original authors, Antonio J. Nebro and Juan J. Durillo. Nevertheless, 
 * this current version can be redistributed and/or modified under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */
//
package jmetal.metaheuristics.moead;

import jmetal.core.*;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.metaheuristics.moead.Utils;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.util.Vector;

public class MOEADACD_many extends Algorithm {

	private int populationSize_;
	
	private SolutionSet population_;	// Stores the population
	
	int nr_;	// maximal number of solutions replaced by each child solution
	int T_;		// neighborhood size
	
	double delta_;	// probability that parent solutions are selected from neighbourhood
	
	double[] z_;	// Z vector (ideal point)
	
	double[][] lambda_;		// weight vectors
	double[][] directionVector_;
	int[][] neighborhood_;	// neighborhood structure
	
	Solution[] indArray_;
	String functionType_;
	int evaluations_;

	Operator crossover_;
	Operator mutation_;

	String dataDirectory_;

	double[] thetas_;
	double[] thetas_min_;
	double[] thetas_max_;
	int kindex_;
	String   algName_;
	String   savePF_;

	/**
	 * Constructor
	 * @param problem: Problem to solve
	 */
	public MOEADACD_many(Problem problem, String functionType) {
		super(problem);

		this.functionType_ = functionType;

	} // MOEA/D

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		int maxEvaluations;
		evaluations_    = 0;
		maxEvaluations  = ((Integer) this.getInputParameter("maxEvaluations")).intValue();
		populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
		dataDirectory_  = this.getInputParameter("dataDirectory").toString();
		population_ = new SolutionSet(populationSize_);
		T_ 	   = ((Integer) this.getInputParameter("T")).intValue();
		nr_    = ((Integer) this.getInputParameter("nr")).intValue();
		delta_ = ((Double) this.getInputParameter("delta")).doubleValue();
		neighborhood_ = new int[populationSize_][T_];
		z_ 		= new double[problem_.getNumberOfObjectives()];
		lambda_ = new double[populationSize_][problem_.getNumberOfObjectives()];
		directionVector_ = new double[populationSize_][problem_.getNumberOfObjectives()];
		crossover_ = operators_.get("crossover");	// default: DE crossover
		mutation_  = operators_.get("mutation"); 	// default: polynomial mutation
		thetas_ = new double[populationSize_];
		thetas_min_ = new double[populationSize_];
		thetas_max_ = new double[populationSize_];
		algName_ = this.getInputParameter("algName").toString();
		savePF_ = this.getInputParameter("savePF").toString();
		if (algName_.equals("MOEADCD"))
		{
			kindex_ = ((Integer) this.getInputParameter("kindex")).intValue();
		}

		// STEP 1. Initialization
		initUniformWeight();
		initNeighborhood();
		initPopulation();
		initIdealPoint();
		
		// STEP 2. Update
		do {
			int[] permutation = new int[populationSize_];
			Utils.randomPermutation(permutation, populationSize_);

			for (int i = 0; i < populationSize_; i++) {
				int n = permutation[i]; // or int n = i;
				int type;
				double rnd = PseudoRandom.randDouble();
				// STEP 2.1. Mating selection based on probability
				if (rnd < delta_)	// if (rnd < realb)
					type = 1; // neighborhood
				else
					type = 2; // whole population
				Vector<Integer> p = new Vector<Integer>();
				matingSelection(p, n, 1, type);
				Solution[] parents   = new Solution[2];
				parents[0] = population_.get(p.get(0));
				parents[1] = population_.get(n);
				// STEP 2.2. Reproduction
				//SBX
				Solution[] offSpring = (Solution[]) crossover_.execute(parents);
				// Apply mutation
				mutation_.execute(offSpring[0]);
				mutation_.execute(offSpring[1]);

				// Evaluation
				problem_.evaluate(offSpring[0]);
				problem_.evaluate(offSpring[1]);
				evaluations_ +=2;

				// STEP 2.3. Repair. If necessary

				// STEP 2.4. Update z_
				updateReference(offSpring[0]);
				updateReference(offSpring[1]);

				// STEP 2.5. Update of solutions
				updateProblem(offSpring[0], n, type);
				updateProblem(offSpring[1], n, type);
						
			} // for
			if (algName_.equals("MOEADACD"))
				updateTheta();
		} while (evaluations_ < maxEvaluations);
		return population_;
	}

	/**
	 * initUniformWeight
	 */
	public void initUniformWeight() {
		if ((problem_.getNumberOfObjectives() == 2)) {
			for (int n = 0; n < populationSize_; n++) {
				double a = 1.0 * n / (populationSize_ - 1);
				directionVector_[n][0] = a;
				directionVector_[n][1] = 1 - a;
			} // for
		} // if
		else {
			String dataFileName;
			dataFileName = "DirectionVector" + problem_.getNumberOfObjectives() + "D_"
					+ populationSize_ + ".dat";
			try {
				// Open the file
				FileInputStream fis = new FileInputStream(dataDirectory_ + "/"
						+ dataFileName);
				InputStreamReader isr = new InputStreamReader(fis);
				BufferedReader br = new BufferedReader(isr);

				int numberOfObjectives = 0;
				int i = 0;
				int j = 0;
				String aux = br.readLine();
				while (aux != null) {
					StringTokenizer st = new StringTokenizer(aux);
					j = 0;
					numberOfObjectives = st.countTokens();
					while (st.hasMoreTokens()) {
						double value = (new Double(st.nextToken()))
								.doubleValue();
						directionVector_[i][j] = value;
						// System.out.println("lambda["+i+","+j+"] = " + value)
						// ;
						j++;
					}
					aux = br.readLine();
					i++;
				}
				br.close();
			} catch (Exception e) {
				System.out
						.println("initUniformWeight: failed when reading for file: "
								+ dataDirectory_ + "/" + dataFileName);
				e.printStackTrace();
			}
		} // else

		if (functionType_.equals("_PBI") || functionType_.equals("_WS"))
		{
			for (int i=0;i<populationSize_ ;i++ )
			{
				for (int j=0;j<problem_.getNumberOfObjectives() ;j++ )
				{
					lambda_[i][j] = directionVector_[i][j];
				}
			}
		}else if (functionType_.equals("_TE"))
		{
			for (int i=0;i<populationSize_ ;i++ )
			{
				for (int j=0;j<problem_.getNumberOfObjectives() ;j++ )
				{
					lambda_[i][j] = 1.0/directionVector_[i][j];
				}
			}
		}
	} // initUniformWeight

	/**
	 * Initialize the neighborhood structure
	 */
	public void initNeighborhood() {
		double[] x = new double[populationSize_];
		int[] idx  = new int[populationSize_];

		for (int i = 0; i < populationSize_; i++) {
			// calculate the distances based on weight vectors
			for (int j = 0; j < populationSize_; j++) {
				x[j]   = norm(minus(lambda_[i],lambda_[j]));
				idx[j] = j;
			} // for

			// find 'niche' nearest neighboring subproblems
			Utils.minFastSort(x, idx, populationSize_, T_);

			System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
		} // for
	} // initNeighborhood

	/**
	 * Initialize the population
	 * @throws JMException
	 * @throws ClassNotFoundException
	 */
	public void initPopulation() throws JMException, ClassNotFoundException {
		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problem_);
			if (algName_.equals("MOEADACD") || algName_.equals("MOEADCD"))
			{
				thetas_min_[i] = angle( directionVector_[i], directionVector_[neighborhood_[i][1]]);
				thetas_max_[i] = 2.0*maxAngle2Axis( directionVector_[i] );
				if (algName_.equals("MOEADACD"))
				{
					thetas_[i] = thetas_min_[i];
				}else if (algName_.equals("MOEADCD"))
				{
					thetas_[i] = thetas_min_[i] + ((double)kindex_/9.0)* (thetas_max_[i] - thetas_min_[i]) ;
				}			
				
			}
			problem_.evaluate(newSolution);
			evaluations_++;
			population_.add(newSolution);
		} // for
	} // initPopulation

	/**
	 * Initialize the ideal objective vector
	 * @throws JMException
	 * @throws ClassNotFoundException
	 */
	void initIdealPoint() throws JMException, ClassNotFoundException {
		for (int i = 0; i < populationSize_; i++) {
			updateReference(population_.get(i));
		} // for
	} // initIdealPoint

	/**
	 * Select mating parents
	 * @param list: the set of the indexes of selected mating parents
	 * @param cid : the id of current subproblem
	 * @param size: the number of selected mating parents
	 * @param type: 1 - neighborhood; otherwise - whole population
	 */
	public void matingSelection(Vector<Integer> list, int cid, int size,
			int type) {
		int ss, r, p;

		ss = neighborhood_[cid].length;
		while (list.size() < size) {
			if (type == 1) {
				r = PseudoRandom.randInt(0, ss - 1);
				p = neighborhood_[cid][r];
			} else {
				p = PseudoRandom.randInt(0, populationSize_ - 1);
			}
			boolean flag = true;
			for (int i = 0; i < list.size(); i++) {
				if (list.get(i) == p) {	// p is in the list
					flag = false;
					break;
				}
			}

			if (flag) {
				list.addElement(p);
			}
		}
	} // matingSelection


	/**
	 * update ideal objective vector
	 * @param individual
	 */

	void updateReference(Solution indiv) {
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			if (indiv.getObjective(i) < z_[i]) {
				z_[i] = indiv.getObjective(i);
			}
		}
	} // updateReference

	void updateNadirPoint(Solution indiv, double[] nz_) {
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			if (indiv.getObjective(i) > nz_[i])
				nz_[i] = indiv.getObjective(i);
		}
	} // updateNadirPoint

	void updateTheta(){
		double[] temp = new double[populationSize_];
		double sum_temp = 0.0;
		for (int i=0; i < populationSize_ ; i++ )
		{
			temp[i] = angle(directionVector_[i], minus(population_.get(i).getObjectives(), z_));
			sum_temp += temp[i];
		}
		for (int i=0; i < populationSize_ ; i++ )
		{
			if ((populationSize_*temp[i]/sum_temp) > 1.0)
			{
				thetas_[i] = thetas_[i] - thetas_min_[i];
				thetas_[i] = thetas_[i] < thetas_min_[i]?thetas_min_[i]:thetas_[i];
			}else if ((populationSize_*temp[i]/sum_temp) < 1.0)
			{
				thetas_[i] = thetas_[i] + thetas_min_[i];
				thetas_[i] = thetas_[i] > thetas_max_[i]?thetas_max_[i]:thetas_[i];
			}
		}
	} // updateTheta

	/**
	 * update the subproblem
	 * @param indiv: child solution
	 * @param id   : the id of current subproblem
	 * @param type : update solutions in - neighborhood (1) or whole population (otherwise)
	 */
	void updateProblem(Solution indiv, int id, int type) {

		int size, time;

		time = 0;

		if (type == 1) {
			size = neighborhood_[id].length;
		} else {
			size = population_.size();
		}
		int[] perm = new int[size];

		randomPermutation(perm, size);

		for (int i = 0; i < size; i++) {
			int k;
			if (type == 1) {
				k = neighborhood_[id][perm[i]];
			} else {
				k = perm[i]; // calculate the values of objective function
			}
			double f1, f2;
			f1 = fitnessFunction(population_.get(k), lambda_[k],directionVector_[k], thetas_[k]);
			f2 = fitnessFunction(indiv, lambda_[k],directionVector_[k], thetas_[k]);

			if (f2 < f1) {
				population_.replace(k, new Solution(indiv));
				time++;
			}
			// the maximal number of solutions updated is not allowed to exceed 'limit'
			if (time >= nr_) {
				return;
			}
		}
	} // updateProblem

	double fitnessFunction(Solution indiv, double[] lambda, double[] directionVector, double theta) {
		double fitness;
		fitness = 0.0;

		if (functionType_.equals("_TE")) {
			double maxFun = -1.0e+30;

			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				double diff = Math.abs(indiv.getObjective(n) - z_[n]);

				double feval;
				if (lambda[n] == 0) {
					feval = 0.0001 * diff;
				} else {
					feval = diff * lambda[n];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for

			fitness = maxFun;
		} else if (functionType_.equals("_PBI")) {
			double beta = 5.0;
			double d1 = dot(minus(indiv.getObjectives(),z_),lambda)/norm(lambda);
			double d2 = norm(minus(minus(indiv.getObjectives(),z_), multiple((d1/norm(lambda)),lambda)));
			fitness = d1 + beta * d2;

		} else {
			System.out.println("MOEAD.fitnessFunction: unknown type "
					+ functionType_);
			System.exit(-1);
		}


		if (algName_.equals("MOEADACD") || algName_.equals("MOEADCD"))
		{
			if (angle(  directionVector, minus(indiv.getObjectives(),z_)  ) > 0.5*theta)
			{
				fitness = Double.POSITIVE_INFINITY;
			}
		}
		return fitness;
	} // fitnessEvaluation

/**
===============================PRIVATE METHODS========================================	
*/

	private double[] minus(double[] vec1,double[] vec2){
		int n = vec1.length;
		double[] value = new double[n];
		for (int i=0;i<n ;i++ )
		{
			value[i] = vec1[i] - vec2[i];
		}
		return value;

	}

	private double[] multiple(double s, double[] vec1){
		int n = vec1.length;
		double[] value = new double[n];
		for (int i=0;i<n ;i++ )
		{
			value[i] = s*vec1[i];
		}
		return value;
	}

	private double dot(double[] vec1, double[] vec2) {
		double sum = 0;

		for (int i = 0; i < vec1.length; i++)
			sum += vec1[i] * vec2[i];

		return sum;
	}

	private double norm(double[] z) {
		double sum = 0;

		for (int i = 0; i < z.length; i++)
			sum += z[i] * z[i];

		return Math.sqrt(sum);
	}

	private double angle(double[] vec1, double[] vec2){
		return Math.acos(dot(vec1,vec2)/(norm(vec1)*norm(vec2)));	
	}

  private double maxAngle2Axis(double[] vector1){
				
	  double max_angle = -1000.0;
	  double angle;
	  double cos_angle;

	  for (int i=0;i<vector1.length ;i++ )
	  {
		  cos_angle = vector1[i]/norm(vector1);
		  cos_angle = cos_angle<-1.0?-1.0:cos_angle;
		  cos_angle = cos_angle>1.0?1.0:cos_angle;
		  angle = Math.acos(cos_angle);
		  if (max_angle < angle)
		  {
			  max_angle = angle;
		  }
	  }
	  return max_angle;
  }



  private void randomPermutation(int[] perm, int size) {
    int[] index = new int[size];
    boolean[] flag = new boolean[size];
    for (int n = 0; n < size; n++) {
      index[n] = n;
      flag[n] = true;
    }

    int num = 0;
    while (num < size) {
      int start = PseudoRandom.randInt(0, size - 1);
      //int start = int(size*nd_uni(&rnd_uni_init));
      while (true) {
        if (flag[start]) {
          perm[num] = index[start];
          flag[start] = false;
          num++;
          break;
        }
        if (start == (size - 1)) {
          start = 0;
        } else {
          start++;
        }
      }
    } // while
  } // randomPermutation
	
} // MOEADACD

