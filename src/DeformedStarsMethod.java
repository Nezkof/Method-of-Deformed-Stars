import java.util.ArrayList;
import java.util.Random;

public class DeformedStarsMethod {
    private final ArrayList<double[]> P_population;
    private ArrayList<double[]> Pz_population;
    private ArrayList<double[]> Ps_population;
    private ArrayList<double[]> Pw_population;

    int populationSize;
    int compressionCoefficient;

    private final int itersNumber;
    private final double populationDistance;
    private final double populationValuesDistance;

    private final Function function;
    private final double[] bounds;
    private final Random random;
    private long duration;
    private double avgResult = Double.MAX_VALUE;
    private double bestResult;
    private double oldAvgResult = 0;


    public DeformedStarsMethod(int populationSize, int compressionCoefficient, int itersNumber, double populationDistance, double populationValuesDistance, double[] bounds,Function function){
        this.random = new Random();
        this.itersNumber = itersNumber;
        this.populationDistance = populationDistance;
        this.populationValuesDistance = populationValuesDistance;
        this.populationSize = populationSize;
        this.compressionCoefficient = compressionCoefficient;
        this.bounds = bounds;

        this.P_population = createPopulation();
        this.Pz_population = new ArrayList<>();
        this.Ps_population = new ArrayList<>();
        this.Pw_population = new ArrayList<>();

        this.function = function;
    }

    public void startOptimization() {
        long startTime = System.currentTimeMillis();

        if (itersNumber != 0) {
            for (int i = 0; i < itersNumber; ++i) {
                createParallelTransferPopulation();
                createRotationPopulation();
                createCompressionPopulation();
                selectBestPopulation();
            }
        } else
        if (populationDistance != 0) {
            while(Math.abs(getIndividsDistance()) > populationDistance){
                createParallelTransferPopulation();
                createRotationPopulation();
                createCompressionPopulation();
                selectBestPopulation();
            }
        } else
        if (populationValuesDistance != 0) {
            while(Math.abs(oldAvgResult - avgResult) > populationValuesDistance) {
                oldAvgResult = avgResult;
                createParallelTransferPopulation();
                createRotationPopulation();
                createCompressionPopulation();
                selectBestPopulation();
            }
        }

        duration = System.currentTimeMillis() - startTime;
    }

    private double getIndividsDistance() {
        return Math.sqrt(Math.pow(P_population.get(0)[0] - P_population.get(1)[0],2) + Math.pow(P_population.get(0)[1] - P_population.get(1)[1],2));
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

    private void createParallelTransferPopulation() {
        Pz_population.clear();

        for (int individ = 0; individ < P_population.size(); ++individ) {
            int i, j;
            double a, alpha;

            do {
                i = random.nextInt(P_population.size());
                j = random.nextInt(P_population.size());
            } while (j == i);

            a = random.nextDouble() * (bounds[1] - bounds[0]) + bounds[0];
            alpha = random.nextDouble() * 2 * Math.PI;


            double[] firstPoint = new double[2];
            firstPoint[0] = P_population.get(i)[0];
            firstPoint[1] = P_population.get(i)[1];

            double[] secondPoint = new double[2];
            secondPoint[0] = P_population.get(j)[0];
            secondPoint[1] = P_population.get(j)[1];

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

        for (int individ = 0; individ < P_population.size(); ++individ) {
            int i, j;
            double beta;

            do {
                i = random.nextInt(P_population.size());
                j = random.nextInt(P_population.size());
            } while (j == i);

            beta = Math.toRadians(random.nextDouble(360));

            double[] firstPoint = new double[2];
            firstPoint[0] = P_population.get(i)[0];
            firstPoint[1] = P_population.get(i)[1];

            double[] secondPoint = new double[2];
            secondPoint[0] = P_population.get(j)[0];
            secondPoint[1] = P_population.get(j)[1];

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

        for (int individ = 0; individ < P_population.size(); ++individ) {
            int i, j;

            do {
                i = random.nextInt(P_population.size());
                j = random.nextInt(P_population.size());
            } while (j == i);

            double[] firstPoint = new double[2];
            firstPoint[0] = P_population.get(i)[0];
            firstPoint[1] = P_population.get(i)[1];

            double[] secondPoint = new double[2];
            secondPoint[0] = P_population.get(j)[0];
            secondPoint[1] = P_population.get(j)[1];

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

        evaluateAvgPopulationScore();
        bestResult = function.calculate(P_population.get(0)[0], P_population.get(0)[1]);

    }

    public void evaluateAvgPopulationScore() {
        avgResult = 0;
        for (var individ : P_population){
            avgResult += function.calculate(individ[0], individ[1]);
        }
        avgResult /= P_population.size();
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

    public int getPopulationSize() {
        return populationSize;
    }

    public long getDuration(){
        return duration;
    }

    public double getAvgResult() {
        return avgResult;
    }

    public int getCompressionCoefficient() {
        return compressionCoefficient;
    }

    public double getPopulationDistance() {
        return populationDistance;
    }

    public double getPopulationValuesDistance() {
        return populationValuesDistance;
    }

    public double getBestResult() {return bestResult;}
}
