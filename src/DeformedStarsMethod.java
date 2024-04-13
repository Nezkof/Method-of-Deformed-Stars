import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class DeformedStarsMethod {
    private final ArrayList<double[]> P_population;
    private final ArrayList<Double> F_population;

    private ArrayList<double[]> Pz_population;
    private final ArrayList<Double> Fz_population;

    private ArrayList<double[]> Ps_population;
    private final ArrayList<Double> Fs_population;

    private ArrayList<double[]> Pw_population;
    private final ArrayList<Double> Fw_population;

    int populationSize;
    int compressionCoefficient;

    private final int itersNumber;
    private final double populationDistance;
    private final double populationValuesDistance;

    private final Function function;
    private final double[] bounds;
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
            createParallelTransferPopulation();
            calculateFunctionValues(Pz_population, Fz_population);
            createRotationPopulation();
            calculateFunctionValues(Ps_population, Fs_population);
            createCompressionPopulation();
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

    private void createParallelTransferPopulation() {
        Pz_population.clear();
        Fz_population.clear();

        for (int individ = 0; individ < P_population.size(); ++individ) {
            int i, j;
            double a, alpha;

            do {
                i = random.nextInt(P_population.size());
                j = random.nextInt(P_population.size());
            } while (j == i);

            a = random.nextDouble() * (bounds[1] - bounds[0]) + bounds[0];
            alpha = random.nextDouble() * 2 * Math.PI;

            double[] firstPoint = P_population.get(i);
            double[] secondPoint = P_population.get(j);

            firstPoint[0] += a * Math.cos(alpha);
            firstPoint[1] += a * Math.sin(alpha);

            firstPoint = validatePoint(firstPoint);

            secondPoint[0] += a * Math.cos(alpha);
            secondPoint[1] += a * Math.sin(alpha);

            secondPoint = validatePoint(secondPoint);

            Pz_population.add(function.calculate(firstPoint[0], firstPoint[1]) < function.calculate(secondPoint[0], secondPoint[1]) ? firstPoint : secondPoint);
        }
    }

    private void createRotationPopulation() {
        Ps_population.clear();
        Fs_population.clear();

        for (int individ = 0; individ < P_population.size(); ++individ) {
            int i, j;
            double beta;

            do {
                i = random.nextInt(P_population.size());
                j = random.nextInt(P_population.size());
            } while (j == i);

            beta = Math.toRadians(random.nextDouble(360));

            double[] firstPoint = P_population.get(i);
            double[] secondPoint = P_population.get(j);

            double xDiff = secondPoint[0] - firstPoint[0];
            double yDiff = secondPoint[1] - firstPoint[1];

            double newX = secondPoint[0] + xDiff * Math.cos(beta) - yDiff * Math.sin(beta);
            double newY = secondPoint[1] + xDiff * Math.sin(beta) + yDiff * Math.cos(beta);

            secondPoint[0] = newX;
            secondPoint[1] = newY;

            firstPoint = validatePoint(firstPoint);
            secondPoint = validatePoint(secondPoint);

            Ps_population.add(function.calculate(firstPoint[0], firstPoint[1]) > function.calculate(secondPoint[0], secondPoint[1]) ? secondPoint : firstPoint);
        }
    }


    private void createCompressionPopulation() {
        Pw_population.clear();
        Fw_population.clear();

        for (int individ = 0; individ < P_population.size(); ++individ) {
            int i, j;

            do {
                i = random.nextInt(P_population.size());
                j = random.nextInt(P_population.size());
            } while (j == i);

            double[] firstPoint = P_population.get(i);
            double[] secondPoint = P_population.get(j);

            if (function.calculate(firstPoint[0], firstPoint[1]) < function.calculate(secondPoint[0], secondPoint[1])) {
                secondPoint[0] = (secondPoint[0] + firstPoint[0]) / compressionCoefficient;
                secondPoint[1] = (secondPoint[1] + firstPoint[1]) / compressionCoefficient;
            } else {
                firstPoint[0] = (firstPoint[0] + secondPoint[0]) / compressionCoefficient;
                firstPoint[1] = (firstPoint[1] + secondPoint[1]) / compressionCoefficient;
            }

            firstPoint = validatePoint(firstPoint);
            secondPoint = validatePoint(secondPoint);

            Pw_population.add(function.calculate(firstPoint[0], firstPoint[1]) > function.calculate(secondPoint[0], secondPoint[1]) ? secondPoint : firstPoint);
        }
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
