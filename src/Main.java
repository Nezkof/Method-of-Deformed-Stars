public class Main {
    public static void main(String[] args) {
        Function[] functions = new Function[] {
                (x, y) -> 20 + (x * x - 10 * Math.cos(2 * Math.PI * x)) + (y * y - 10 * Math.cos(2 * Math.PI * y)), // Rastrigin Function
                (x, y) -> Math.pow(Math.sin(3 * Math.PI * x), 2) + Math.pow(x - 1, 2) * (1 + Math.pow(Math.sin(3 * Math.PI * y), 2)) +  Math.pow(y - 1, 2) * (1 + Math.pow(Math.sin(2 * Math.PI * y), 2)), // Levy 13 Function
                (x, y) -> 0.5 * (Math.pow(x, 4) - 16 * x * x + 5 * x + Math.pow(y, 4) - 16 * y * y + 5 * y) // Stibinski-Tanga Function
        };
        double[][] functionsBounds ={ {-5.12, 5.12}, {-10, 10}, {-5, 5} };
        double[][] tableArguments = { {0,0}, {1,1}, {-2.903534,-2.903534}};

        int[] populationSize = new int[] {10, 20 ,30};
        int[] compressionCoefficient =  new int[] {2,3,4};
        int iterNumbers = 0;
        double populationDistance = 0.000001;
        double populationValuesDistance = 0;

        for (int i = 0; i < functions.length; ++i) {
            System.out.println("=====================================================================");
            System.out.print("\t\t\t\t\tTesting for <");
            switch (i) {
                case 0:
                    System.out.print("Rastrigin function");
                    break;
                case 1:
                    System.out.print("Levy 13 function");
                    break;
                case 2:
                    System.out.print("Stibinski-Tanga function");
                    break;
            }
            System.out.println(">");
            System.out.println("---------------------------------------------------------------------");
            System.out.printf("%17s%17s%18s%17s\n", "Precision", "Duration", "Result", "Error");
            for (int j = 0; j < 3; ++j) {
                DeformedStarsMethod algorithm = new DeformedStarsMethod(populationSize[2], compressionCoefficient[0], iterNumbers, populationDistance, populationValuesDistance, functionsBounds[i], functions[i]);
                algorithm.startOptimization();
                System.out.printf("%17s%17.5s%18.7f%17.7f\n", algorithm.getPopulationDistance(), algorithm.getDuration() + "ms", algorithm.getBestResult(), getError(algorithm.getBestResult(), functions[i].calculate(tableArguments[i][0], tableArguments[i][1])));

            }
        }
    }

    public static double getError(double bestScore, double tableValue){
        return bestScore - tableValue;
    };
}