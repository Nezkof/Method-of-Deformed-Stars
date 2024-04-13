import javax.security.auth.kerberos.KerberosTicket;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.function.BinaryOperator;

public class DeformedStarsMethod {
    private ArrayList<double[]> P_population;
    private ArrayList<Double> F_population;

    private ArrayList<double[]> Pz_population;
    private ArrayList<Double> Fz_population;

    private ArrayList<double[]> Ps_population;
    private ArrayList<Double> Fs_population;

    private ArrayList<double[]> Pw_population;
    private ArrayList<Double> Fw_population;

    int populationSize;
    int compressionCoefficient;

    private int itersNumber;
    private double populationDistance;
    private double populationValuesDistance;

    private final Function function;
    private double[] bounds;
    private final Random random;
    private double bestResult;

    public DeformedStarsMethod(int populationSize, int compressionCoefficient, int itersNumber, double populationDistance, double populationValuesDistance, double[] bounds,Function function){
        this.random = new Random();
        this.itersNumber = itersNumber;
        this.populationDistance = populationDistance;
        this.populationValuesDistance = populationValuesDistance;
        this.populationSize = populationSize;
        this.compressionCoefficient = compressionCoefficient;
        this.bounds = bounds;

        this.P_population = createPopulation();
        this.F_population = new ArrayList<>();
        this.Pz_population = new ArrayList<>();
        this.Fz_population = new ArrayList<>();
        this.Ps_population = new ArrayList<>();
        this.Fs_population = new ArrayList<>();
        this.Pw_population = new ArrayList<>();
        this.Fw_population = new ArrayList<>();

        this.function = function;
    }

    public void startOptimization(){
        for (int i = 0; i < itersNumber; ++i) {
            calculateFunctionValues(P_population, F_population);
            this.Pz_population = createParallelTransferPopulation();
            calculateFunctionValues(Pz_population, Fz_population);
            this.Ps_population = createRotationPopulation();
            calculateFunctionValues(Ps_population, Fs_population);
            this.Pw_population = createCompressionPopulation();
            calculateFunctionValues(Pw_population, Fw_population);
            selectBestPopulation();
            System.out.println(bestResult);
        }
    }

    private ArrayList<double[]> createPopulation() {
        ArrayList<double[]> population = new ArrayList<>();

        for (int i = 0; i < populationSize; ++i) {
            population.add(new double[] {
                    random.nextDouble(bounds[1] - bounds[0] + 1) + bounds[0],
                    random.nextDouble(bounds[1] - bounds[0] + 1) + bounds[0]
            });
        }
        return population;
    }

    private void calculateFunctionValues(ArrayList<double[]> population, ArrayList<Double> results) {
        for (double[] doubles : population) {
            results.add(function.calculate(doubles[0], doubles[1]));
        }
    }

    private ArrayList<double[]> createParallelTransferPopulation() {
        ArrayList<double[]> population = new ArrayList<>();

        for (int individ = 0; individ < P_population.size(); ++individ) {
            int i = random.nextInt(P_population.size());
            int j;
            double a = random.nextDouble(bounds[1]/2.0);
            double alpha = Math.toRadians(random.nextDouble(180));

            do {
                j = random.nextInt(P_population.size());
            } while (j == i);

            double[] firstPoint = P_population.get(i);
            double[] secondPoint = P_population.get(j);

            firstPoint[0] += a * Math.cos(alpha);
            firstPoint[1] += a * Math.sin(alpha);

            firstPoint = validatePoint(firstPoint);

            secondPoint[0] += a * Math.cos(alpha);
            secondPoint[1] += a * Math.sin(alpha);

            secondPoint = validatePoint(secondPoint);

            double firstPointValue = function.calculate(firstPoint[0], firstPoint[1]);
            double secondPointValue = function.calculate(secondPoint[0], secondPoint[0]);

            if (firstPointValue > secondPointValue)
                population.add(secondPoint);
            else
                population.add(firstPoint);

        }

        return population;
    }

    private ArrayList<double[]> createRotationPopulation() {
        ArrayList<double[]> population = new ArrayList<>();

        for (int individ = 0; individ < P_population.size(); ++individ) {
            int i = random.nextInt(P_population.size());
            int j;
            double beta = Math.toRadians(random.nextDouble(360));

            do {
                j = random.nextInt(P_population.size());
            } while (j == i);

            double[] firstPoint = P_population.get(i);
            double[] secondPoint = P_population.get(j);

            if (function.calculate(firstPoint[0], firstPoint[1]) > function.calculate(secondPoint[0], secondPoint[0])) {
                secondPoint[0] += (secondPoint[0] - firstPoint[0])*Math.cos(beta) - (secondPoint[1] - firstPoint[1])*Math.sin(beta);
                beta = Math.toRadians(random.nextDouble(360));
                secondPoint[0] += (secondPoint[0] - firstPoint[0])*Math.sin(beta) - (secondPoint[1] - firstPoint[1])*Math.cos(beta);
            }
            else {
                firstPoint[0] += (firstPoint[0] - secondPoint[0])*Math.cos(beta) - (firstPoint[1] - secondPoint[1])*Math.sin(beta);
                beta = Math.toRadians(random.nextDouble(360));
                firstPoint[0] += (firstPoint[0] - secondPoint[0])*Math.sin(beta) - (firstPoint[1] - secondPoint[1])*Math.cos(beta);
            }

            firstPoint = validatePoint(firstPoint);
            secondPoint = validatePoint(secondPoint);

            if (function.calculate(firstPoint[0], firstPoint[1]) > function.calculate(secondPoint[0], secondPoint[0]))
                population.add(secondPoint);
            else
                population.add(firstPoint);
        }

        return population;
    };

    private ArrayList<double[]> createCompressionPopulation() {
        ArrayList<double[]> population = new ArrayList<>();

        for (int individ = 0; individ < P_population.size(); ++individ) {
            int i = random.nextInt(P_population.size());
            int j;

            do {
                j = random.nextInt(P_population.size());
            } while (j == i);

            double[] firstPoint = P_population.get(i);
            double[] secondPoint = P_population.get(j);

            if (function.calculate(firstPoint[0], firstPoint[1]) < function.calculate(secondPoint[0], secondPoint[0])) {
                secondPoint[0] = (secondPoint[0] + firstPoint[0])/compressionCoefficient;
                secondPoint[1] = (secondPoint[1] + firstPoint[1])/compressionCoefficient;
            }
            else {
                firstPoint[0] = (firstPoint[0] + secondPoint[0])/compressionCoefficient;
                firstPoint[1] = (firstPoint[1] + secondPoint[1])/compressionCoefficient;
            }

            firstPoint = validatePoint(firstPoint);
            secondPoint = validatePoint(secondPoint);

            if (function.calculate(firstPoint[0], firstPoint[1]) > function.calculate(secondPoint[0], secondPoint[0]))
                population.add(secondPoint);
            else
                population.add(firstPoint);
        }

        return population;
    }

    private void selectBestPopulation() {
        ArrayList<double[]> combinedPopulation = new ArrayList<>();
        combinedPopulation.addAll(P_population);
        combinedPopulation.addAll(Pz_population);
        combinedPopulation.addAll(Ps_population);
        combinedPopulation.addAll(Pw_population);

        combinedPopulation.sort((o1, o2) -> Double.compare(function.calculate(o1[0], o1[1]), function.calculate(o2[0], o2[1])));

        P_population.clear();
        P_population.addAll(combinedPopulation.subList(0, populationSize));

        this.bestResult = function.calculate(P_population.get(0)[0], P_population.get(0)[1]);
    }
    public double[] validatePoint(double[] point) {
        if (point[0] < bounds[0] && point[1] < bounds[0]) {
            point[0] -= (bounds[0] - bounds[1]);
            point[1] -= (bounds[0] - bounds[1]);
        }

        if (point[0] < bounds[0] && point[1] > bounds[1]) {
            point[0] -= (bounds[0] - bounds[1]);
            point[1] += (bounds[0] - bounds[1]);
        }

        if (point[0] > bounds[1] && point[1] > bounds[1]) {
            point[0] += (bounds[0] - bounds[1]);
            point[1] += (bounds[0] - bounds[1]);
        }

        if (point[0] > bounds[1] && point[1] < bounds[0]) {
            point[0] += (bounds[0] - bounds[1]);
            point[1] -= (bounds[0] - bounds[1]);
        }

        if(point[0] >= bounds[0] && point[0] <= bounds[1]) {
            if (point[1] > bounds[1])
                point[1] += bounds[0] - bounds[1];
            if (point[1] < bounds[0])
                point[1] -= bounds[0] - bounds[1];
        }

        if (point[1] >= bounds[0] && point[1] <= bounds[1]) {
            if (point[0] > bounds[1])
                point[0] += bounds[0] - bounds[1];
            if (point[0] < bounds[0])
                point[0] -= bounds[0] - bounds[1];
        }

        return point;
    }
}
