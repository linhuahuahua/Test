package jmetal.metaheuristics.MOEAD;

import jmetal.core.*;
import jmetal.util.PseudoRandom;

public class MOEAD_CA extends Algorithm{
	
	public MOEAD_CA(Problem problem){
		super(problem);
	}
	
	public SolutonSet execute() throws JMException ,ClassNotFoundException {
		
		int maxEvaluations;
		int type;
		int delta_;
		int populationSize;
		int evaluations_=0;
		
		initUniformWeight();
		initNeighborhood();
		initPopulation();
		initIdealPoint();
		do{
			
			
			int[] permutation=new int[populationSize];
			Utils.randPermutation(permutation, populationSize);
			
			for(int i=0;i<populationSize;i++){
				
				double rnd=PseudoRandom.randDouble();
				if(rnd<delta_){
					type=1;
				}
					
				else{
					type=2;
				}
					
				
				
			}
			
			
			
			
		} while(evaluations_ < maxEvaluations);
	}

}
