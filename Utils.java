//  Utils.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.metaheuristics.moead;

/**
 * Utilities methods to used by MOEA/D
 */
public class Utils {

  public static double[] vectorAdd(double[] vector1, double[] vector2) {
	  
    int dim = vector1.length;
	double[] value = new double[dim];
    for (int n = 0; n < dim; n++) {
      value[n] = vector1[n] + vector2[n];
    }
    return value;
  } 

  public static double[] vectorMinus(double[] vector1, double[] vector2) {
    int dim = vector1.length;
	double[] value = new double[dim];
    for (int n = 0; n < dim; n++) {
      value[n] = vector1[n] - vector2[n];
    }
    return value;
  } 
  public static double vectorDot(double[] vector1, double[] vector2) {
    int dim = vector1.length;
	double value = 0.0;
    for (int n = 0; n < dim; n++) {
      value += vector1[n]*vector2[n];
    }
    return value;
  } 
  public static double distVector(double[] vector1, double[] vector2) {
    int dim = vector1.length;
    double sum = 0;
    for (int n = 0; n < dim; n++) {
      sum += (vector1[n] - vector2[n]) * (vector1[n] - vector2[n]);
    }
    return Math.sqrt(sum);
  } // distVector

  public static double vectorNorm(double[] vector1){
	double norm_value = 0.0;
	for (int i=0;i<vector1.length ;i++ )
	{
		norm_value += Math.pow(vector1[i],2);
	}
	return Math.sqrt(norm_value);
  }

  public static double vectorAngle(double[] vector1, double[] vector2){
	  double cos_angle = 0.0;
	  cos_angle = vectorDot(vector1,vector2)/ ( vectorNorm(vector1) * vectorNorm(vector2) );
	  cos_angle = cos_angle<-1.0?-1.0:cos_angle;
	  cos_angle = cos_angle>1.0?1.0:cos_angle;
	  return Math.acos(cos_angle);
  }

  public static double vector_maxAngle2Axis(double[] vector1){
				
	  double max_angle = -1000.0;
	  double angle;
	  double cos_angle;

	  for (int i=0;i<vector1.length ;i++ )
	  {
		  cos_angle = vector1[i]/vectorNorm(vector1);
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

  public static void minFastSort(double x[], int idx[], int n, int m) {
    for (int i = 0; i < m; i++) {
      for (int j = i + 1; j < n; j++) {
        if (x[i] > x[j]) {
          double temp = x[i];
          x[i] = x[j];
          x[j] = temp;
          int id = idx[i];
          idx[i] = idx[j];
          idx[j] = id;
        } // if
      }
    } // for

  } // minFastSort

  public static void randomPermutation(int[] perm, int size) {
    int[] index = new int[size];
    boolean[] flag = new boolean[size];

    for (int n = 0; n < size; n++) {
      index[n] = n;
      flag[n] = true;
    }

    int num = 0;
    while (num < size) {
      int start = jmetal.util.PseudoRandom.randInt(0, size - 1);
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



}
