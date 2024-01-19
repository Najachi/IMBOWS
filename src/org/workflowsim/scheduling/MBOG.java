// Importer les bibliothèques nécessaires
package org.workflowsim.scheduling;

import java.util.Random;

// Classe principale MBOG
public class MBOG {
    // Attributs de classe
    private int dims;
    private MinimizationFunction minFunc;
    private double tolerance;
    private boolean terminationConditionMet;
    private int popSize;
    private int t;
    private Butterfly[] butterflyPopulation;
    private Butterfly[] np1;
    private Butterfly[] np2;
    private int maxGen;
    private int sMax;
    private double bar;
    private double peri;
    private double p;
    private int numElite;
    private Butterfly[] elite;

    // Constructeur de classe
    public MBOG(int dims, MinimizationFunction minFunc, int popSize, int maxGen, int sMax, double bar, double peri, double p, int numElite) {
        this.dims = dims;
        this.minFunc = minFunc;
        this.tolerance = 1e-3;
        this.terminationConditionMet = false;
        this.popSize = popSize;
        this.t = 1;
        this.butterflyPopulation = null;
        this.np1 = null;
        this.np2 = null;
        this.maxGen = maxGen;
        this.sMax = sMax;
        this.bar = bar;
        this.peri = peri;
        this.p = p;
        this.numElite = numElite;
        this.elite = null;
    }

    // Méthode pour initialiser la population de papillons
    private Butterfly[] initializePopulation(int popSize, int dim) {
        Butterfly[] population = new Butterfly[popSize];
        for (int i = 0; i < popSize; i++) {
            double[] x = new double[dim];
            for (int j = 0; j < dim; j++) {
                x[j] = new Random().nextGaussian();
            }
            population[i] = new Butterfly(x);
        }
        return population;
    }

    // Méthode pour trier la population de papillons
    private Butterfly[] sortPopulation(Butterfly[] pop) {
        for (int i = 0; i < pop.length - 1; i++) {
            for (int j = 0; j < pop.length - i - 1; j++) {
                if (pop[j].getFitness() > pop[j + 1].getFitness()) {
                    Butterfly temp = pop[j];
                    pop[j] = pop[j + 1];
                    pop[j + 1] = temp;
                }
            }
        }
        return pop;
    }

    // Méthode pour évaluer la fitness de chaque papillon
    private Butterfly[] fitnessEval(Butterfly[] pop) {
        for (int i = 0; i < pop.length; i++) {
            double fitness = minFunc.evaluate(pop[i].getX());
            pop[i].setFitness(fitness);
        }
        return pop;
    }

    // Méthode pour diviser la population en deux groupes
    private void splitPopulation() {
        int splitIndex = (int) Math.ceil(p * popSize);
        np1 = new Butterfly[splitIndex];
        np2 = new Butterfly[popSize - splitIndex];
        System.arraycopy(butterflyPopulation, 0, np1, 0, splitIndex);
        System.arraycopy(butterflyPopulation, splitIndex, np2, 0, popSize - splitIndex);
    }

    // Méthode pour effectuer l'opérateur de migration
    private void migrationOperator() {
        for (int i = 0; i < np1.length; i++) {
            Butterfly butterfly = np1[i];
            for (int k = 0; k < butterfly.getX().length; k++) {
                if (k == 0) {
                    continue;
                } else {
                    double rand = new Random().nextDouble();
                    if (rand * peri <= p) {
                        Butterfly randomButterfly = np1[new Random().nextInt(np1.length)];
                        butterfly.getX()[k] = randomButterfly.getX()[k];
                    } else {
                        Butterfly randomButterfly = np2[new Random().nextInt(np2.length)];
                        butterfly.getX()[k] = randomButterfly.getX()[k];
                    }
                }
            }
            double fitness = minFunc.evaluate(butterfly.getX());
            if (fitness < butterfly.getFitness()) {
                butterfly.setFitness(fitness);
            }
        }
        np1 = fitnessEval(np1);
        np1 = sortPopulation(np1);
    }

    // Méthode pour effectuer l'opérateur d'ajustement des papillons
    private void butterflyAdjustingOperator() {
        for (int i = 0; i < np2.length; i++) {
            Butterfly butterfly = np2[i];
            double dx = levy();
            double alpha = sMax / Math.pow(t, 2);
            for (int k = 0; k < butterfly.getX().length; k++) {
                if (k == 0) {
                    continue;
                } else {
                    double rand = new Random().nextDouble();
                    if (rand <= p) {
                        Butterfly bestButterfly = np1[0].getFitness() < np2[0].getFitness() ? np1[0] : np2[0];
                        butterfly.getX()[k] = bestButterfly.getX()[k];
                    } else {
                        Butterfly randomButterfly = np2[new Random().nextInt(np2.length)];
                        butterfly.getX()[k] = randomButterfly.getX()[k];
                        if (rand > bar) {
                            butterfly.getX()[k] += alpha * (dx - 0.5);
                        }
                    }
                }
            }
            double fitness = minFunc.evaluate(butterfly.getX());
            if (fitness < butterfly.getFitness()) {
                butterfly.setFitness(fitness);
            }
        }
        np2 = fitnessEval(np2);
        np2 = sortPopulation(np2);
    }

    // Méthode pour effectuer un vol de Lévy
    private double levy() {
        double sum = 0;
        for (int i = 0; i < sMax; i++) {
            sum += Math.tan(Math.PI * new Random().nextDouble());
        }
        return sum;
    }

    // Méthode pour combiner les populations np1 et np2
    private void combinePopulation() {
        butterflyPopulation = new Butterfly[np1.length + np2.length];
        System.arraycopy(np1, 0, butterflyPopulation, 0, np1.length);
        System.arraycopy(np2, 0, butterflyPopulation, np1.length, np2.length);
    }

    // Méthode pour vérifier la condition d'arrêt
    private void checkTerminationCondition() {
        double delta = Math.abs(butterflyPopulation[0].getFitness() - elite[0].getFitness());
        if (delta <= tolerance && t != 1) {
            terminationConditionMet = true;
        }
    }

    // Méthode principale pour exécuter l'algorithme MBO
    private void mainLoop() {
        while (!terminationConditionMet && t <= maxGen) {
            butterflyPopulation = sortPopulation(butterflyPopulation);
            elitism("apply");
            elitism("save");
            checkTerminationCondition();
            splitPopulation();
            migrationOperator();
            butterflyAdjustingOperator();
            combinePopulation();
            t++;
        }
    }

    // Méthode pour gérer l'élitisme
    private void elitism(String action) {
        if (action.equals("save")) {
            elite = new Butterfly[numElite];
            for (int i = 0; i < numElite; i++) {
                elite[i] = butterflyPopulation[i];
            }
        } else if (elite != null && action.equals("apply")) {
            for (int i = 0; i < elite.length; i++) {
                elite[elite.length - 1 - i] = elite[i];
            }
        }
    }

    // Méthode pour exécuter l'algorithme MBO et renvoyer la meilleure solution trouvée
    public double[] run() {
        butterflyPopulation = initializePopulation(popSize, dims);
        butterflyPopulation = fitnessEval(butterflyPopulation);
        mainLoop();
        return sortPopulation(butterflyPopulation)[0].getX();
    }

    // Méthode principale pour exécuter l'algorithme MBO avec les paramètres spécifiés
    public static void main(String[] args) {
        // Spécifier les paramètres de l'algorithme MBO
        String funcType = "f6";
        int dims = 2;
        int n = 10;
        int ncpu = 4;

        // Créer une instance de la fonction de minimisation
        MinimizationFunction minFunc = new MinimizationFunction(funcType);

        // Créer une instance de l'algorithme MBO
        MBOG mbog = new MBOG(dims, minFunc, n, ncpu);

        // Afficher un message de démarrage
        System.out.println("\n");
        System.out.println("L'algorithme MBO démarre sur la fonction " + funcType + "...");

        // Exécuter l'algorithme MBO plusieurs fois et stocker les résultats
        double[][] results = new double[n][];
        for (int i = 0; i < n; i++) {
            results[i] = mbog.run();
        }

        // Calculer la moyenne des résultats
        double[] avgR = averageResults(results, dims);

        // Afficher la solution trouvée
        System.out.println("Solution trouvée!\n");
        System.out.println(avgR[0] + " @ x=" + avgR[1] + "\n");
        System.out.println("TERMINÉ!\n");
    }

    // Méthode pour calculer la moyenne des résultats
    private static double[] averageResults(double[][] results, int dim) {
        double[] avg = new double[dim + 1];
        for (double[] r : results) {
            for (int i = 0; i < dim + 1; i++) {
                avg[i] += r[i];
            }
        }
        for (int i = 0; i < dim + 1; i++) {
            avg[i] /= results.length;
        }
        return avg;
    }
}

// Classe Butterfly pour représenter un individu dans la population de papillons
class Butterfly {
    private double[] x;
    private double fitness;

    public Butterfly(double[] x) {
        this.x = x;
        this.fitness = Double.MAX_VALUE;
    }

    public double[] getX() {
        return x;
    }

    public double getFitness() {
        return fitness;
    }

    public void setFitness(double fitness) {
        this.fitness = fitness;
    }
}

// Classe MinimizationFunction pour évaluer les fonctions de minimisation spécifiques
class MinimizationFunction {
    private String funcType;
    private double bounds;

    public MinimizationFunction(String funcType) {
        this.funcType = funcType;
        this.bounds = 0;
        setBounds();
    }

    public double evaluate(double[] x) {
        switch (funcType) {
            case "f6":
                return f6(x);
            case "f7":
                return f7(x);
            case "ackley":
                return ackley(x);
            case "schwefel":
                return schwefel(x);
            case "rastrigin":
                return rastrigin(x);
            default:
                return 0;
        }
    }

    private void setBounds() {
        switch (funcType) {
            case "f6":
                bounds = 100;
                break;
            case "f7":
                bounds = 100;
                break;
            case "ackley":
                bounds = 32.768;
                break;
            case "schwefel":
                bounds = 500;
                break;
            case "rastrigin":
                bounds = 5.12;
                break;
        }
    }

    public double getBounds() {
        return bounds;
    }

    private double f6(double[] x) {
        double sumXSquared = 0;
        for (double xi : x) {
            sumXSquared += Math.pow(xi, 2);
        }
        return 0.5 + (Math.sin(Math.pow(Math.sqrt(sumXSquared), 2)) - 0.5) / Math.pow(1 + 0.001 * sumXSquared, 2);
    }

    private double f7(double[] x) {
        double normalizer = 1.0 / (x.length - 1);
        double fitness = 0;
        for (int i = 0; i < x.length - 1; i++) {
            double si = Math.pow(x[i], 2) + Math.pow(x[i + 1], 2);
            fitness += normalizer * Math.sqrt(si) * Math.pow(Math.sin(50.0 * Math.pow(si, 0.2)) + 1, 2);
        }
        return fitness;
    }

    private double ackley(double[] x) {
        double a = 20;
        double b = 0.2;
        double c = 2 * Math.PI;
        double sumXSquared = 0;
        double sumXCos = 0;
        for (double xi : x) {
            sumXSquared += Math.pow(xi, 2);
            sumXCos += Math.cos(c * xi);
        }
        return -a * Math.exp(-b * Math.sqrt(1.0 / x.length * sumXSquared)) - Math.exp(1.0 / x.length * sumXCos) + a + Math.exp(1);
    }

    private double schwefel(double[] x) {
        double sum = 0;
        for (double xi : x) {
            sum += xi * Math.sin(Math.sqrt(Math.abs(xi)));
        }
        return 418.9829 * x.length - sum;
    }

    private double rastrigin(double[] x) {
        double sum = 0;
        for (double xi : x) {
            sum += Math.pow(xi, 2) - 10 * Math.cos(2 * Math.PI * xi);
        }
        return 10 * x.length + sum;
    }
}