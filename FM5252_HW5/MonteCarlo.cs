using System;
using System.Linq;
class MonteCarlo
{
    public static void Main(string[] args)
    {
        // Define variables
        Console.Write("input initial price:");
        double s0 = double.Parse(Console.ReadLine());
        Console.Write("input strike price:");
        double k = double.Parse(Console.ReadLine());
        Console.Write("input risk-free rate:");
        double r = double.Parse(Console.ReadLine());
        Console.Write("input volatility:");
        double sigma = double.Parse(Console.ReadLine());
        Console.Write("input time maturity:");
        double t = double.Parse(Console.ReadLine());
        Console.Write("input steps: ");
        int steps = int.Parse(Console.ReadLine());
        Console.Write("input simulations: ");
        int n = int.Parse(Console.ReadLine());

        double dt = t / steps;
        double[] s = new double[n];

        // Define results output
        double[] callPayoffs = new double[n];
        double[] putPayoffs = new double[n];
        double[] cDeltaUp = new double[n];
        double[] cDeltaDown = new double[n];
        double[] pDeltaUp = new double[n];
        double[] pDeltaDown = new double[n];
        double[] cVegaUp = new double[n];
        double[] cVegaDown = new double[n];
        double[] pVegaUp = new double[n];
        double[] pVegaDown = new double[n];
        double[] cThetaUp = new double[n];
        double[] pThetaUp = new double[n];
        double[] cRhoUp = new double[n];
        double[] cRhoDown = new double[n];
        double[] pRhoUp = new double[n];
        double[] pRhoDown = new double[n];
        

        // Generate random numbers
        Random rand = new Random();

        // Generate pertubation simulations
        double pertubation = rand.NextDouble();

        for (int i = 0; i < n; i++)
        {
            // Generate normal distributed numbers
            double z = rand.NextDouble();
            double x = Math.Sqrt(-2.0 * Math.Log(z)) * Math.Cos(2.0 * Math.PI * rand.NextDouble());

            s[i] = s0 * Math.Exp((r - 0.5 * Math.Pow(sigma, 2)) * t + sigma * Math.Sqrt(t) * x);
            callPayoffs[i] = Math.Max(s[i] - k, 0.0);
            putPayoffs[i] = Math.Max(k - s[i], 0.0);

            // Elements for calculate Delta and Gamma
            cDeltaUp[i] = Math.Max((s0 + pertubation) * Math.Exp((r - 0.5 * Math.Pow(sigma, 2)) * t + sigma * Math.Sqrt(t) * x) - k, 0);
            cDeltaDown[i] = Math.Max((s0 - pertubation) * Math.Exp((r - 0.5 * Math.Pow(sigma, 2)) * t + sigma * Math.Sqrt(t) * x) - k, 0);
            pDeltaUp[i] = Math.Max(k - (s0 + pertubation) * Math.Exp((r - 0.5 * Math.Pow(sigma, 2)) * t + sigma * Math.Sqrt(t) * x), 0);
            pDeltaDown[i] = Math.Max(k - (s0 - pertubation) * Math.Exp((r - 0.5 * Math.Pow(sigma, 2)) * t + sigma * Math.Sqrt(t) * x), 0);

            // Elements for calculate Vega
            cVegaUp[i] = Math.Max(s0 * Math.Exp(r - 0.5 * Math.Pow(sigma + pertubation, 2) * t + (sigma + pertubation) * Math.Sqrt(t) * x) - k, 0);
            cVegaDown[i] = Math.Max(s0 * Math.Exp(r - 0.5 * Math.Pow(sigma - pertubation, 2) * t + (sigma - pertubation) * Math.Sqrt(t) * x) - k, 0);
            pVegaUp[i] = Math.Max(k - s0 * Math.Exp(r - 0.5 * Math.Pow(sigma + pertubation, 2) * t + (sigma + pertubation) * Math.Sqrt(t) * x), 0);
            pVegaDown[i] = Math.Max(k - s0 * Math.Exp(r - 0.5 * Math.Pow(sigma - pertubation, 2) * t + (sigma - pertubation) * Math.Sqrt(t) * x), 0);

            // Elements for calculate Theta
            cThetaUp[i] = Math.Max(s0 * Math.Exp((r - 0.5 * Math.Pow(sigma, 2)) * (t + pertubation) + sigma * Math.Sqrt(t + pertubation) * x) - k, 0);
            pThetaUp[i] = Math.Max(k - s0 * Math.Exp((r - 0.5 * Math.Pow(sigma, 2)) * (t + pertubation) + sigma * Math.Sqrt(t + pertubation) * x), 0);

            // Elements for calculate Rho
            cRhoUp[i] = Math.Max(s0 * Math.Exp((r + pertubation) - 0.5 * Math.Pow(sigma, 2) * t + sigma * Math.Sqrt(t) * x) - k, 0);
            cRhoDown[i] = Math.Max(s0 * Math.Exp((r - pertubation) - 0.5 * Math.Pow(sigma, 2) * t + sigma * Math.Sqrt(t) * x) - k, 0);
            pRhoUp[i] = Math.Max(k - s0 * Math.Exp((r + pertubation) - 0.5 * Math.Pow(sigma, 2) * t + sigma * Math.Sqrt(t) * x), 0);
            pRhoDown[i] = Math.Max(k - s0 * Math.Exp((r - pertubation) - 0.5 * Math.Pow(sigma, 2) * t + sigma * Math.Sqrt(t) * x), 0);

        }

        // Calculate call and put prices
        double callPrice = callPayoffs.Average() * Math.Exp(-r * t);
        double putPrice = putPayoffs.Average() * Math.Exp(-r * t);

        // Calculate call and put standard deviations
        double callStdDev = Math.Sqrt((callPayoffs.Select(p => (p - callPrice) * (p - callPrice)).Sum()) / (n - 1)) / Math.Sqrt(n);
        double putStdDev = Math.Sqrt((putPayoffs.Select(p => (p - putPrice) * (p - putPrice)).Sum()) / (n - 1)) / Math.Sqrt(n);

        // Calculate call and put Greeks
        //Delta
        double callDelta = (cDeltaUp.Average() - cDeltaDown.Average()) / (2 * pertubation);
        double putDelta = (pDeltaUp.Average() - pDeltaDown.Average()) / (2 * pertubation);

        //Gamma
        double callGamma = (cDeltaUp.Average() - 2 * callPayoffs.Average() + cDeltaDown.Average()) / Math.Pow(pertubation, 2);
        double putGamma = (pDeltaUp.Average() - 2 * putPayoffs.Average() + pDeltaDown.Average()) / Math.Pow(pertubation, 2);

        //Vega
        double callVega = (cVegaUp.Average() - cVegaDown.Average()) / (2 * pertubation);
        double putVega = (pVegaUp.Average() - pVegaDown.Average()) / (2 * pertubation);

        //Theta
        double callTheta = (cThetaUp.Average() - callPayoffs.Average()) / pertubation;
        double putTheta = (pThetaUp.Average() - putPayoffs.Average()) / pertubation;

        //Rho
        double callRho = (cRhoUp.Average() - cRhoDown.Average()) / (2 * pertubation);
        double putRho = (pRhoUp.Average() - pRhoDown.Average()) / (2 * pertubation);


        Console.WriteLine("Call price: {0}", callPrice);
        Console.WriteLine("Put price: {0}", putPrice);
        Console.WriteLine("Call standard deviation: {0}", callStdDev);
        Console.WriteLine("Put standard deviation: {0}", putStdDev);
        Console.WriteLine("Call Delta:{0}", callDelta);
        Console.WriteLine("Put Delta:{0}", putDelta);
        Console.WriteLine("Call Gamma: {0}", callGamma);
        Console.WriteLine("Put Gamma: {0}", putGamma);
        Console.WriteLine("Call Vega: {0}", callVega);
        Console.WriteLine("Put Vega: {0}", putVega);
        Console.WriteLine("Call Theta: {0}", callTheta);
        Console.WriteLine("Put Theta: {0}", putTheta);
        Console.WriteLine("Call Rho: {0}", callRho);
        Console.WriteLine("Put Rho: {0}", putRho);
        Console.ReadKey();
    }

}
