using System;
using System.Linq;

class MonteCarlo
{
    public static void Main(string[] args)
    {
        // Define variables
        Console.Write("Input option type; 'call' or 'put'");
        string optionType = Console.ReadLine();
        if (optionType != "call" && optionType != "put") 
        {
            Console.WriteLine("Invalid option type. Please input 'call' or 'put'.");
            return;
        }
        Console.Write("Input initial price:");
        double s0 = double.Parse(Console.ReadLine());
        Console.Write("Input strike price:");
        double k = double.Parse(Console.ReadLine());
        Console.Write("Input risk-free rate:");
        double r = double.Parse(Console.ReadLine());
        Console.Write("Input volatility:");
        double sigma = double.Parse(Console.ReadLine());
        Console.Write("Input time maturity:");
        double t = double.Parse(Console.ReadLine());
        Console.Write("Input steps: ");
        int steps = int.Parse(Console.ReadLine());
        Console.Write("Input simulations: ");
        int n = int.Parse(Console.ReadLine());
        Console.Write("Input base number 1 for Van Der Corput: ");
        double baseNumber1 = double.Parse(Console.ReadLine());
        Console.Write("Input base number 2 for Van Der Corput: ");
        double baseNumber2 = double.Parse(Console.ReadLine());
        Console.Write("Apply antithetic?: 'true' or 'false'");
        bool applyAntithetic = bool.Parse(Console.ReadLine());

        double dt = t / steps;
        double[] s = new double[n];
        double[] sAnti = new double[n];

        // Define Van Der Corput sequences
        double baseNumber = 0;
        double[] rand = new double[2*n];
        double[] rand1 = new double[n];
        double[] rand2 = new double[n];
        double[] vdc1 = new double[n];
        double[] vdc2 = new double[n];

        // Define results output
        double[] payoff = new double[n];
        double sumPayoffs = 0;
        double sumPayoffSquared = 0;
        double stdErrorAnti = 0;
        double stdError = 0;
        double optionPrice = 0;
        double delta = 0;
        double gamma = 0;
        double theta = 0;
        double vega = 0;
        double rho = 0;

        // Run simulations
        if (applyAntithetic = true)
        {
            for (int i = 0; i < n; i++)
            {
                double[] payoffAnti = new double[n];
                double[] payoffAvg = new double[n];

                vdc1[i] = VanDerCorput(i + 1, baseNumber1);
                vdc2[i] = VanDerCorput(i + 1, baseNumber2);
                rand1[i] = Math.Sqrt(-2.0 * Math.Log(vdc1[i])) * Math.Cos(2.0 * Math.PI * vdc2[i]);
                rand2[i] = Math.Sqrt(-2.0 * Math.Log(vdc1[i])) * Math.Sin(2.0 * Math.PI * vdc2[i]);
                Array.Copy(rand1, rand, rand1.Length);
                Array.Copy(rand2, 0, rand, rand1.Length, rand2.Length);

                s = GeometricBrownian(s0,k,r,sigma,dt,n,rand);
                sAnti = GeometricBrownianAnti(s0,k,r,sigma,dt,n,rand);

                payoff[i] = optionType == "call" ? Math.Max(s[i] - k, 0) : Math.Max(k - s[i], 0);
                payoffAnti[i] = optionType == "call" ? Math.Max(sAnti[i] - k, 0) : Math.Max(k - sAnti[i], 0);
                payoffAvg[i] = (payoff[i] + payoffAnti[i]) / 2;

                sumPayoffs += payoffAvg[i];
                sumPayoffSquared += Math.Pow(payoffAvg[i], 2);
            }
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                vdc1[i] = VanDerCorput(i + 1, baseNumber1);
                vdc2[i] = VanDerCorput(i + 1, baseNumber2);
                rand1[i] = Math.Sqrt(-2.0 * Math.Log(vdc1[i])) * Math.Cos(2.0 * Math.PI * vdc2[i]);
                rand2[i] = Math.Sqrt(-2.0 * Math.Log(vdc1[i])) * Math.Sin(2.0 * Math.PI * vdc2[i]);
                Array.Copy(rand1, rand, rand1.Length);
                Array.Copy(rand2, 0, rand, rand1.Length, rand2.Length);

                s = GeometricBrownian(s0,k,r,sigma,dt,n,rand);

                payoff[i] = optionType == "call" ? Math.Max(s[i] - k, 0) : Math.Max(k - s[i], 0);

                sumPayoffs += payoff[i];
                sumPayoffSquared += Math.Pow(payoff[i], 2);
            }
        }

        optionPrice = payoff.Average() * Math.Exp(-r * t);
        stdErrorAnti = Math.Sqrt((sumPayoffSquared / (2 * n)) - Math.Pow(sumPayoffs / (2 * n), 2) / (n - 1)) / Math.Sqrt((2 * n));
        stdError = Math.Sqrt((payoff.Select(p => (p - optionPrice) * (p - optionPrice)).Sum()) / ((2 * n) - 1)) / Math.Sqrt((2 * n));
        delta = GetDelta(s0,k,r,sigma,dt,n,rand,optionType).Average();
        gamma = GetGamma(s0,k,r,sigma,dt,n,rand,optionType).Average();
        vega = GetVega(s0,k,r,sigma,dt,n,rand,optionType).Average();
        theta = GetTheta(s0,k,r,sigma,dt,n,rand,optionType).Average();
        rho = GetRho(s0,k,r,sigma,dt,n,rand,optionType).Average();

        Console.WriteLine("Option price: {0}", optionPrice);
        Console.WriteLine("Standard Error without antithetic: {0}", stdError);
        Console.WriteLine("Standard Error with antithetic: {0}", stdErrorAnti);
        Console.WriteLine("Delta:{0}", delta);
        Console.WriteLine("Gamma: {0}", gamma);
        Console.WriteLine("Vega: {0}", vega);
        Console.WriteLine("Theta: {0}", theta);
        Console.WriteLine("Rho: {0}", rho);
        Console.ReadKey();
    }
    public static double VanDerCorput(int n, double baseNumber)     
    {
        double result = 0;
        double f = 1.0 / baseNumber;

        while (n > 0)
        {
            result += f * (n % baseNumber);
            n /= (int) baseNumber;
            f /= baseNumber;
        }  
        return result;
    }

    public static double[] GeometricBrownian(double s0, double k, double r, double sigma, double dt, int n, double[] rand)
    {
        double[] s = new double[2*n];
        for (int i = 0; i < rand.Length; i++)
        {
            dt = dt < 0 ? Math.Abs(dt) : dt;
            s[i] = s0 * Math.Exp((r - 0.5 * Math.Pow(sigma, 2)) * dt + sigma * Math.Sqrt(dt) * rand[i]);
        }
        return s;
    }
    public static double[] GeometricBrownianAnti(double s0, double k, double r, double sigma, double dt, int n, double[] rand)
    {
        double[] sAnti = new double[2*n];
        for (int i = 0; i < rand.Length; i++)
        {
            sAnti[i] = s0 * Math.Exp((r - 0.5 * Math.Pow(sigma, 2)) * dt - sigma * Math.Sqrt(dt) * rand[i]);
        }
        return sAnti;
    }
    public static double[] GetDelta(double s0, double k, double r, double sigma, double dt, int n, double[] rand, string optionType)
    {
        double[] delta = new double[2*n];
        double[] deltaUp = new double[2*n];
        double[] deltaDown = new double[2*n];
        double[] randDouble = new double[2*n];

        for (int i = 0; i < rand.Length; i++)
        {
            deltaUp = GeometricBrownian(s0 + rand[i], k, r, sigma, dt, n, rand);
            deltaDown = GeometricBrownian(s0 - rand[i], k, r, sigma, dt, n, rand);

            randDouble[i] = rand[i] * 2;

            delta[i] = optionType == "call" ? (Math.Max(deltaUp[i] - k, 0) - Math.Max(deltaDown[i] - k, 0)) / randDouble[i] : (Math.Max(k - deltaUp[i], 0) - Math.Max(k - deltaDown[i], 0)) / randDouble[i];
        }
        return delta;
    }
    public static double[] GetGamma(double s0, double k, double r, double sigma, double dt, int n, double[] rand, string optionType)
    {
        double[] gamma = new double[2*n];
        double[] gammaUp = new double[2*n];
        double[] gammaDown = new double[2*n];
        double[] randSquare = new double[2*n];

        for (int i = 0; i < rand.Length; i++)
        {
            gammaUp = GeometricBrownian(s0 + rand[i], k, r, sigma, dt, n, rand);
            gammaDown = GeometricBrownian(s0 - rand[i], k, r, sigma, dt, n, rand);

            randSquare[i] = Math.Pow(rand[i], 2);

            gamma[i] = optionType == "call" ? (Math.Max(gammaUp[i] - k, 0) + Math.Max(gammaDown[i] - k, 0) - 2 * Math.Max(GeometricBrownian(s0, k, r, sigma, dt, n, rand)[i] - k, 0)) / randSquare[i] : (Math.Max(k - gammaUp[i], 0) + Math.Max(k - gammaDown[i], 0) - 2 * Math.Max(k - GeometricBrownian(s0, k, r, sigma, dt, n, rand)[i], 0)) / randSquare[i];
        }
        return gamma;
    }

    public static double[] GetVega(double s0, double k, double r, double sigma, double dt, int n, double[] rand, string optionType)
    {
        double[] vega = new double[2*n];
        double[] vegaUp = new double[2*n];
        double[] vegaDown = new double[2*n];
        double[] randDouble = new double[2*n];

        for (int i = 0; i < rand.Length; i++)
        {
            vegaUp = GeometricBrownian(s0, k, r, sigma + rand[i], dt, n, rand);
            vegaDown = GeometricBrownian(s0, k, r, sigma - rand[i], dt, n, rand);

            randDouble[i] = rand[i] * 2;

            vega[i] = optionType == "call" ? (Math.Max(vegaUp[i] - k, 0) - Math.Max(vegaDown[i] - k, 0)) / randDouble[i] : (Math.Max(k - vegaUp[i], 0) - Math.Max(k - vegaDown[i], 0)) / randDouble[i];
        }
        return vega;
    }

    public static double[] GetTheta(double s0, double k, double r, double sigma, double dt, int n, double[] rand, string optionType)
    {
        double[] theta = new double[2*n];
        double[] thetaUp = new double[2*n];
        double[] thetaDown = new double[2*n];

        for (int i = 0; i < rand.Length; i++)
        {
            thetaUp = GeometricBrownian(s0, k, r, sigma, dt + rand[i], n, rand);
            thetaDown = GeometricBrownian(s0, k, r, sigma, dt, n, rand);
            
            theta[i] = optionType == "call" ? (Math.Max(thetaUp[i] - k, 0) - Math.Max(thetaDown[i] - k, 0)) / rand[i] : (Math.Max(k - thetaUp[i], 0) - Math.Max(k - thetaDown[i], 0)) / rand[i];          
        }
        return theta;
    }

    public static double[] GetRho(double s0, double k, double r, double sigma, double dt, int n, double[] rand, string optionType)
    {
        double[] rho = new double[2*n];
        double[] rhoUp = new double[2*n];
        double[] rhoDown = new double[2*n];
        double[] randDouble = new double[2*n];

        for (int i = 0; i < rand.Length; i++)
        {
            rhoUp = GeometricBrownian(s0, k, r  + rand[i], sigma, dt, n, rand);
            rhoDown = GeometricBrownian(s0, k, r - rand[i], sigma, dt, n, rand);

            randDouble[i] = rand[i] * 2;

            rho[i] = optionType == "call" ? (Math.Max(rhoUp[i] - k, 0) - Math.Max(rhoDown[i] - k, 0)) / randDouble[i] : (Math.Max(k - rhoUp[i], 0) - Math.Max(k - rhoDown[i], 0)) / randDouble[i];
        }
        return rho;
    }
}
