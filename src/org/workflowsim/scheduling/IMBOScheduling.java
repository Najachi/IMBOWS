package org.workflowsim.scheduling;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;

import org.workflowsim.WorkflowEngine;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.RandomGenerator;
import org.workflowsim.WorkflowEngine;

public class IMBOScheduling{
	
	// Paramètres classique de la metaheuristique IMBO
    public static int populationSize; // Taille de la population de papillons = nombre de tâches
    public static int maxIterations; // Nombre maximal d'itérations
    public static double migrationRatio; // Ratio de migration entre Land1 et Land2
    public static double adjustmentRate; // Taux d'ajustement des papillons
    public static double maxStep; //Smax ou taux de pas
    public static double migrationPeriod; //Taux de migration
    public static double lowerBound = 0.1; // Borne inférieure de l'espace de recherche
    public static double upperBound = 0.9; // Borne supérieure de l'espace de recherche
    public static int dimension; // Dimension du problème
    public static Random random; // Générateur de nombres aléatoires
    
    //Paramètre d'ordonnancmeent
    public static List<int[]> population=new ArrayList<int[]>(); //Liste des papillons à optimiser (premier mappage T-VM)
    public static List<int[]> newPopulation=new ArrayList<int[]>();
    public static int initFlag=0;
    
    public static List<int[]> pbest_schedule=new ArrayList<int[]>();
	public static int[] gbest_schedule;
	public static double[] pbest_fitness;
	public static double gbest_fitness=Double.MAX_VALUE;
	public static int[] pbestSchedule;
	public static List<int[]> newSchedules=new ArrayList<int[]>();
    
    //Paramètres de l'environnement
    public static int taskNum;
    public static int vmNum;
   
    
    /**Génération de la population initiale**/
    public static void initializePopulation(int jobNum, int maxVmNum) {
    	taskNum = jobNum;
    	vmNum = maxVmNum;
    	int[] schedule=new int[taskNum];
    	
    	//Itérer sur la taille de la population initiale
        for (int i = 0; i < populationSize; i++) {
        	//Générer des positions aléatoires pour les taches
        	for(int j = 0; j < taskNum; j++) {
        		schedule[j]=new Random().nextInt(vmNum);
        	}
        	population.add(schedule);
        	// Calculer la fitness du papillon
        	//double fitness = calculateFitness(schedule);
        	
        	// Mettre à jour la fitness du papillon
        	//solution.get(i).setFitness(fitness);
        }
        initFlag = 1;
    }
    
    /**Calcul de la dominance de Pareto**/
    /**public static void epsilonFuzzyDominance(List<int[]> population) {
    	for(int l = 0; l<population.size(); l++) {
    		List<int[]> pareto = new ArrayList<int[]>();
    		for(int w = 0; w<population.size(); w++) {
    			if(l==w) {
    				w++;
    			}
    			
    			int u = 1;
    			
    			//if(population[l]>population[w]) {
    				
    			//}
    		}
    	}
    }**/
    
    
    public static List<int[]> epsilonFuzzyDominate(List<int[]> population, double epsilon) {
        List<int[]> dominatedSolutions = new ArrayList<int[]>();
        List<int[]> solutiontemp = population;

        // Construit une pile de solutions non dominées
        PriorityQueue<int[]> nonDominatedSolutions = new PriorityQueue<int[]>((a, b) -> a.length);
        for (int i = 0; i< population.size(); i++){
            nonDominatedSolutions.add(solutiontemp.get(i));
        }

        // Parcours la pile de solutions non dominées
        while (!nonDominatedSolutions.isEmpty()) {
             solutiontemp = nonDominatedSolutions.poll();

            // Vérifie si la solution est dominée par une autre solution
            for (Solution otherSolution : population) {
                if (epsilonFuzzyDominates(solution, otherSolution, epsilon)) {
                    dominatedSolutions.add(solutiontemp);
                    break;
                }
            }
        }

        return dominatedSolutions;
    }
    
   
    
    /**Mise à jour de la population**/
	public static void updateParticles() {
		int SP1, SP2; //les sous-populations
		int stepSize;
		double[] deltaX;
		int lambda = 20;
		int numvar = 20;
		List<int[]> population1 = new ArrayList<int[]>(); //SPopulation1 à T = i
		List<int[]> population2 = new ArrayList<int[]>(); //SPopulation2 à T = i
		List<int[]> Newpopulation1 = new ArrayList<int[]>(); //SPopulation1 à T= i+1
		List<int[]> Newpopulation2 = new ArrayList<int[]>(); //SPopulation2 à T = i+1
		double r1, r2, r3, r4, alpha;
		double randomNumber = random.nextDouble();
		double SelfAdaptiveStrategy, a, b;
		
		//Calcul de SP1 et SP2
		SP1 = (int) Math.ceil(migrationRatio * populationSize);
		SP2 = populationSize - SP1;
		
		//RealMatrix Land1 = new Array2DRowRealMatrix(SP1, numvar).scalarMultiply(0);
		//RealMatrix Land2 = new Array2DRowRealMatrix(SP2, numvar).scalarMultiply(0);
		int [][] Land1 = new int [SP1][numvar];
		int [][] Land2 = new int [SP2][numvar];
		
		//RealMatrix Population1 = Land1;
		//RealMatrix Population2 = Land2;
		//double FitnessPopulation1; //= WorkflowEngine.caculatefitness();
		//double FitnessPopulation2;
		
		a = (lowerBound * maxIterations - upperBound) / maxIterations-1;
		b = (upperBound - lowerBound) / maxIterations - 1;
		
		//calcul du coefficient de pondération alpha et du walk step de la solution j
		alpha = maxStep / (maxIterations^2);
		stepSize = (int) Math.ceil(exprnd(2 * maxIterations));
        deltaX = levyFlight(stepSize, lambda);
		
		//Creation de la SP1 et de la SP2
		for(int i=0; i<populationSize; i++) {
			if(i<=SP1) {
				population1.add(population.get(i));
			}
			else
				population2.add(population.get(i));
		}
		
		
		/**Mise à jour de le SP1 par le BMO et du Self-Adaptive Strategy**/
			for(int i = 0; i<=SP1; i++) {
				for(int k = 0; k<populationSize; k++) {
					r1 = randomNumber * migrationPeriod;
				
					//Calcul du Self-Adaptive Strategy
					SelfAdaptiveStrategy = a + b*i;
				
					if(r1<=SelfAdaptiveStrategy) {
						r2 = (int) Math.round(SP1 * randomNumber + 0.5);
						//Land1.setEntry(i,k, solution1.get((int)r2)[i]);
						Land1[i][k] = population1.get((int)r2)[i];
					
					}else {
						r3 = (int) Math.round(SP2* randomNumber + 0.5);
						//Land1.setEntry(i,k, solution2.get((int)r3)[i]);
						Land1[i][k] = population2.get((int)r3)[i];
					}
				}
				//Calcul du greedy strategy
				//if(getCrowdingDistance(Land1)<getCrowdingDistance(solution1)) {
					
				//}
				//Newsolution1.addAll(Land1);
				//Newpopulation1.addAll(Land1);
			}
		
			/**Mise à jour de la SP2 par le Butterfly Adjusting Operator**/
			for(int j = 0; j<SP2; j++) {
                SelfAdaptiveStrategy = a + b*j;
                
                for(int k = 0; k<populationSize; k++) {
                	if(randomNumber<=SelfAdaptiveStrategy) {
                		//equation15
                		Land2[j][k] = population.get(k)[j];
                	}
                	else if(randomNumber>SelfAdaptiveStrategy && randomNumber<adjustmentRate) {
                		r4 = (int) Math.round(SP2 * randomNumber + 0.5);
                		//equation16;
                		Land2[j][k] = population2.get((int)r4)[j];
                		
                	}else if(randomNumber>SelfAdaptiveStrategy && randomNumber>adjustmentRate) {
                		//equation 18;
                		Land2[j][k] = (Land2[j][k] + alpha*(deltaX-0.5));
                	}
                }
                
                //Newpopulation2.add(Land2);
			}
			
			/**Fusionner les deux SP*/
			for(int i = 0; i<SP1; i++) {
				newPopulation.addAll(Newpopulation1);
			}
			
			for(int i = SP1; i<SP2; i++) {
				newPopulation.addAll(Newpopulation2);
			}
		}
    
	/**Fonction du vol de Levy**/
	 public static double[] levyFlight(int stepSize, int dim) {
	     //Allocation d'une matrice pour la solution   
		 double[] deltaX = new double[dim];
		 
	        for (int i = 0; i < dim; i++) {
	        	
	        	//distribution de cauchy
	            double[] fx = new double[stepSize];
	            for (int j = 0; j < stepSize; j++) {
	                fx[j] = Math.tan(Math.PI * Math.random());
	            }
	            
	            //somme de la distribution
	            deltaX[i] = sum(fx);
	        }
	        return deltaX;
	    }
	 
	 /**Fonction qui calcule la somme**/
	 public static double sum(double[] array) {
	        double sum = 0.0;
	        for (double value : array) {
	            sum += value;
	        }
	        return sum;
	    }
	 
	 /**Fonction qui genère une valeur aléatoire exponentielle*/
	 public static double exprnd(double lambda) {

		    // Generate a random number from a uniform distribution
		    double u = Math.random();

		    // Calculate the exponential distribution
		    return -lambda * Math.log(u);
		}
	 
	 public static void clear() {
			gbest_fitness = Double.MAX_VALUE;
		    initFlag = 0;
	        population.removeAll(population);
	        pbest_schedule.removeAll(pbest_schedule);
	        //velocity.removeAll(velocity);
	        newPopulation.removeAll(newPopulation);
	        pbest_schedule.removeAll(pbest_schedule);
		}
}