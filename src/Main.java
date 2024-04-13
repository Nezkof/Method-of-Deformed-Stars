public class Main {
    public static void main(String[] args) {
        Function[] functions = new Function[] {
                (x, y) -> 20 + (x * x - 10 * Math.cos(2 * Math.PI * x)) + (y * y - 10 * Math.cos(2 * Math.PI * y)), // Rastrigin Function
                (x, y) -> Math.pow(Math.sin(3 * Math.PI * x), 2) + Math.pow(x - 1, 2) * (1 + Math.pow(Math.sin(3 * Math.PI * y), 2)) +  Math.pow(y - 1, 2) * (1 + Math.pow(Math.sin(2 * Math.PI * y), 2)), // Levy 13 Function
                (x, y) -> 0.5 * (Math.pow(x, 4) - 16 * x * x + 5 * x + Math.pow(y, 4) - 16 * y * y + 5 * y) // Stibinski-Tanga Function
        };
        double[][] functionsBounds ={ {-5.12, 5.12}, {-10, 10}, {-5, 5} };
        double[][] tableArguments = { {0,0}, {1,1}, {-2.903534,-2.903534}};

        int populationSize = 10;
        int compressionCoefficient = 2;
        int iterNumbers = 30;
        double populationDistance = 0;
        double populationValuesDistance = 0;

        DeformedStarsMethod algorithm = new DeformedStarsMethod(populationSize, compressionCoefficient, iterNumbers, populationDistance, populationValuesDistance, functionsBounds[2], functions[2]);
        algorithm.startOptimization();
    }
}