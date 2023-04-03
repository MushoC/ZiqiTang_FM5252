using System;
using MathNet.Numerics.Random;

namespace NormalRandomGeneratorApp
{
    class NormalRandomGenerator
    {
        public static double SumTwelve(Random rand)
        {
            double sum = 0;

            for (int i = 0; i < 12; i++)
            {
                sum += rand.NextDouble();
            }

            return sum - 6;           
        }

        public static Tuple<double,double> BoxMuller(Random rand)
        {
            double u1 = rand.NextDouble();
            double u2 = rand.NextDouble();

            double z1 = Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);
            double z2 = Math.Sqrt(-2 * Math.Log(u1)) * Math.Sin(2 * Math.PI * u2);

            return Tuple.Create<double,double>(z1,z2);
        }

        public static Tuple<double,double> PolarRejection(Random rand)
        {
            double u1, u2, s;

            do
            {
                u1 = rand.NextDouble() ;
                u2 = rand.NextDouble() ;
                s = u1 * u1 + u2 * u2;
            } while (s >= 1 || s == 0);

            double z1 = u1 * Math.Sqrt(-2 * Math.Log(s) / s);
            double z2 = u2 * Math.Sqrt(-2 * Math.Log(s) / s);
            
            return Tuple.Create<double,double>(z1,z2);

        }

        public static Tuple<double, double> GeneratePairBoxMuller(Random rand, double correlation)
        {
            double u1 = rand.NextDouble();
            double u2 = rand.NextDouble();

            double z1 = Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);
            double z2 = Math.Sqrt(-2 * Math.Log(u1)) * Math.Sin(2 * Math.PI * u2);

            double x = z1;
            double y = z2 * Math.Sqrt(1 - correlation * correlation) + z1 * correlation;

            return Tuple.Create(x, y);
        }

        public static Tuple<double, double> GeneratePairPolarRejection(Random rand, double correlation)
        {
            double u1, u2, z;

            do
            {
                u1 = rand.NextDouble() ;
                u2 = rand.NextDouble() ;
                z = u1 * u1 + u2 * u2;
            } while (z >= 1 || z == 0);

            z = Math.Sqrt(-2 * Math.Log(z) / z);

            double z1 = u1 * z;
            double z2 = z * Math.Sqrt(1 - correlation * correlation) + z1 * correlation;

            return Tuple.Create(z1, z2);
        }

    }

    class Program
    {
        static void Main(string[] args)
        {
            RandomSource source = new SystemRandomSource();
            int seed = Guid.NewGuid().GetHashCode();
            Random rand = new Random(seed);

            double sum12 = NormalRandomGenerator.SumTwelve(rand);
            Console.WriteLine("Value of Sum of 12:{0}",sum12);

            Tuple<double, double> z1z2 = NormalRandomGenerator.BoxMuller(rand);
            double z1 = z1z2.Item1;
            double z2 = z1z2.Item2;
            Console.WriteLine("Box-Muller: z1={0}, z2={1}", z1, z2);

            Tuple<double, double> z1z2PR = NormalRandomGenerator.PolarRejection(rand);
            double z1PR = z1z2PR.Item1;
            double z2PR = z1z2PR.Item2;
            Console.WriteLine("Polar Rejection: z1={0}, z2={1}", z1PR, z2PR);

            Console.WriteLine("Enter correlation (-1 to 1):");
            double correlation = double.Parse(Console.ReadLine());

            Tuple<double, double> pairb = NormalRandomGenerator.GeneratePairBoxMuller(rand, correlation);
            Console.WriteLine("Joint normally distributed random value for Box Muller x: {0}", pairb.Item1);
            Console.WriteLine("Joint normally distributed random value for Box Muller y: {0}", pairb.Item2);

            Tuple<double, double> pairp = NormalRandomGenerator.GeneratePairBoxMuller(rand, correlation);
            Console.WriteLine("Joint normally distributed random value for Polar Rejection x: {0}", pairp.Item1);
            Console.WriteLine("Joint normally distributed random value for Polar Rejection y: {0}", pairp.Item2);

            Console.ReadKey();
        }
    }
}
