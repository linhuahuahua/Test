

package jmetal.problems.F;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.Scanner;

import jmetal.core.*;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;
import jmetal.util.wrapper.XReal;

public class F1 extends Problem {

	/**
	 * Constructor. Creates a default instance of problem CEC2009_UF1 (31
	 * decision variables)
	 * 
	 * @param solutionType
	 *            The solution type must "Real" or "BinaryReal".
	 */
	public F1(String solutionType) throws ClassNotFoundException {
		this(solutionType, 31); // 31 variables by default
	} 
	  
	public F1(String solutionType, Integer numberOfVariables)
			throws ClassNotFoundException {
		numberOfVariables_ = numberOfVariables.intValue();
		numberOfObjectives_ = 2;
		numberOfConstraints_ = 0;
		problemName_ = "F1";
		
		upperLimit_ = new double[numberOfVariables_];
		lowerLimit_ = new double[numberOfVariables_];

		for (int var = 0; var < numberOfVariables_; var++) {
			if (var<numberOfObjectives_-1){
				lowerLimit_[var] = 0.0;
				upperLimit_[var] = 1.0;
			}
			else{
				lowerLimit_[var] = -1.0;
				upperLimit_[var] = 1.0;
			}

		} // for
		
		if (solutionType.compareTo("BinaryReal") == 0)
	    	solutionType_ = new BinaryRealSolutionType(this) ;
	    else if (solutionType.compareTo("Real") == 0)
	    	solutionType_ = new RealSolutionType(this) ;
	    else {
	    	System.out.println("Error: solution type " + solutionType + " invalid") ;
	    	System.exit(-1) ;
	    }
	} // F1

	/**
	 * Evaluates a solution.
	 * @param solution The solution to evaluate.
	 * @throws JMException
	 */
	public void evaluate(Solution solution) throws JMException {
		double A=1.0;
		double B=4.0;
		double C=2.0;
		double D=4.0;
		double E=3.0;
		double F=1.0;
		
		XReal x = new XReal(solution) ;
		//File file =new File("I:\\data.txt");
		//FileReader fr=new FileReader(file);
		try {
            FileInputStream fos = new FileInputStream("I:\\matlab\\x.txt");
            InputStreamReader osw = new InputStreamReader(fos, "UTF-8");
            BufferedReader bw = new BufferedReader(osw);
            
            for (int i=0;i<numberOfVariables_;i++){
            	String ss=bw.readLine();
            	x.setValue(i, Double.parseDouble(ss));
            }
            bw.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
		
	    double[] y = new double[numberOfObjectives_];
	    for(int k=0;k<numberOfObjectives_;k++){
	    	if (k==0){
	    		y[k]=1-x.getValue(k);
	    	}
	    	else if(k==numberOfObjectives_-1){
	    		y[k]=1.0;
	    		for(int i=0;i<k;i++){
	    			y[k]=y[k]*x.getValue(i);
	    		}
	    	}
	    	else{
	    		double temp=1.0;
	    		for(int i=0;i<k;i++){
	    			temp=temp*x.getValue(i);
	    		}
	    		y[k]=(1-x.getValue(k))*temp;
	    	}
	    }
	    
	    double[] h = new double[numberOfObjectives_];
	    for(int i=0;i<numberOfObjectives_;i++){
	    	h[i]=Math.pow(y[i], F);
	    }
	    
	    double n=(double)numberOfVariables_;
	    double m=(double)numberOfObjectives_;
	    double[][] N=new double[numberOfObjectives_][numberOfObjectives_];
	    for (int i=0;i<numberOfObjectives_;i++){
		    for (int j=0;j<numberOfObjectives_;j++){
		    	if(i==j){
		    		N[i][j]=-1/m;
		    	}
		    	else{
		    		N[i][j]=1/(m*(m-1));
		    	}
		    }
	    }
	    
	    double[] nn=new double[numberOfObjectives_];
	    double n_max=-1.0e+6;
	    for (int i=0;i<numberOfObjectives_;i++){
	    	for (int j=0;j<numberOfObjectives_;j++){
	    		nn[i]=nn[i]+N[i][j]*y[j];
	    	}
	    	if(nn[i]>n_max){
	    		n_max=nn[i];
	    	}
	    }
	    double R=1/Math.sqrt(m*(m-1));
	    double L=n_max/R/R;
	    
	    double b_D=Math.pow(Math.sin(Math.PI/2*Math.pow(L, m-1)),D);
	    double b_B=Math.pow(Math.sin(Math.PI/2*Math.pow(L, m-1)),B);
	    
	    double[] g = new double[numberOfObjectives_];
	    int size_J=(numberOfVariables_-numberOfObjectives_+1)/numberOfObjectives_;
	    int[][] J=new int[numberOfObjectives_][size_J];
	    for (int k=0;k<numberOfObjectives_;k++){
	    	double[] t=new double[size_J];
	    	for(int j=0;j<size_J;j++){
	    		J[k][j]=(j+1)*numberOfObjectives_+k-1;
	    		t[j]=x.getValue(J[k][j])-0.9*b_B*Math.cos(E*Math.PI*L+(n+2)*(double)(J[k][j]+1)*Math.PI/(2*n));
	    		//System.out.println(x.getValue(J[k][j]));
	    		//System.out.println(Math.cos(E*Math.PI*L+(n+2)*(double)(J[k][j]+1)*Math.PI/(2*n)));
	    	}
	    	double temp=0.0;
	    	for (int i=0;i<size_J;i++){
	    		temp=temp+Math.pow(Math.abs(t[i]),C);
	    	}
	    	g[k]=(A*b_D+1)/(double)size_J*temp;
	    }
	    
	    double[] f=new double[numberOfObjectives_];
	    for(int k=0;k<numberOfObjectives_;k++){
	    	f[k]=g[k]+h[k];
	    	solution.setObjective(k, f[k]);
	    	System.out.println(f[k]);
	    }
	    System.out.println("=========");

	} // evaluate
}
