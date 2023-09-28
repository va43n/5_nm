using System;
using System.IO;
using System.Diagnostics;

namespace _5_nm_1 
{
    class Program 
    {
        static void Main() 
        {
            //TEST N
            //int[] N = { 2000 };
            int[] N = { 200, 400, 800 };

            int time, numberOfTests = 10;

            //double[] check = { 0 };

            double delta;

            StartParameters sp = new();
            LU_Decomposition lu = new();

            Stopwatch stopwatch = new();

            foreach (int Ni in N) 
            {
                //-----------------------START PARAMETERS-----------------------
                sp.x = StartParameters.GenerateX(Ni);
                //CHECK X
                //Console.WriteLine("x: ");
                //foreach (double xi in sp.x)
                //{
                //    Console.Write("{0} ", xi);
                //}
                //Console.WriteLine();

                sp.A = StartParameters.GenerateA(Ni);
                //sp.A = StartParameters.GenerateAFromFile(Ni);
                //CHECK A
                //Console.WriteLine("A: ");
                //foreach (double a in sp.A)
                //{
                //    Console.Write("{0} ", a);
                //}
                //Console.WriteLine();

                sp.f = StartParameters.GenerateF(Ni, sp.A);
                //CHECK f
                //Console.WriteLine("f: ");
                //foreach (double fi in sp.f)
                //{
                //    Console.Write("{0} ", fi);
                //}
                //Console.WriteLine();
                //-----------------------START PARAMETERS-----------------------

                time = 0;
                for (int i = 0; i < numberOfTests; i++)
                {
                    stopwatch.Restart();

                    lu.U = LU_Decomposition.GenerateU(Ni, sp.A);
                    //CHECK U
                    //Console.WriteLine("U: ");
                    //foreach (double u in lu.U)
                    //{
                    //    Console.Write("{0} ", u);
                    //}
                    //Console.WriteLine();

                    lu.L = LU_Decomposition.GenerateL(Ni, sp.A, lu.U);
                    //CHECK L
                    //Console.WriteLine("L: ");
                    //foreach (double l in lu.L)
                    //{
                    //    Console.Write("{0} ", l);
                    //}
                    //Console.WriteLine();

                    //CHECK = L * U
                    //Array.Resize(ref check, Ni * Ni);
                    //for (int k = 0; k < Ni; k++)
                    //{
                    //    for (int i = 0; i < Ni; i++)
                    //    {
                    //        check[k * Ni + i] = 0;
                    //        for (int j = 0; j < Ni; j++)
                    //        {
                    //            check[k * Ni + i] += lu.L[k * Ni + j] * lu.U[j * Ni + i];
                    //        }
                    //    }
                    //}
                    //L * U = A?
                    //Console.WriteLine("check: ");
                    //foreach (double c in check)
                    //{
                    //    Console.Write("{0} ", c);
                    //}
                    //Console.WriteLine();

                    //for (int i = 0; i < Ni; i++)
                    //{
                    //    for (int j = 0; j < Ni; j++)
                    //    {
                    //        if (Math.Abs(sp.A[i * Ni + j] - check[i * Ni + j]) < 1)
                    //        {
                    //            Console.WriteLine("ok");
                    //        }
                    //        else
                    //        {
                    //            Console.WriteLine("bad: {0}, {1}", sp.A[i * Ni + j], check[i * Ni + j]);
                    //        }
                    //    }
                    //}

                    lu.y = LU_Decomposition.GenerateY(Ni, lu.L, sp.f);
                    //CHECK y
                    //Console.WriteLine("y: ");
                    //foreach (double yi in lu.y)
                    //{
                    //    Console.Write("{0} ", yi);
                    //}
                    //Console.WriteLine();

                    lu.x = LU_Decomposition.GenerateResult(Ni, lu.U, lu.y);
                    //CHECK LU RESULT
                    //Console.WriteLine("x: ");
                    //foreach (double xi in lu.x)
                    //{
                    //    Console.Write("{0} ", xi);
                    //}
                    //Console.WriteLine();

                    //CHECK GAUSS RESULT
                    //lu.x = LU_Decomposition.GenerateGaussResult(Ni, sp.A, sp.f);
                    //Console.WriteLine("x: ");
                    //foreach (double xi in lu.x)
                    //{
                    //    Console.Write("{0} ", xi);
                    //}
                    //Console.WriteLine();

                    stopwatch.Stop();

                    time += Convert.ToInt32(stopwatch.ElapsedMilliseconds);
                }
                Console.WriteLine(SupportingFunctions.SubstractVectors(Ni, sp.x, lu.x));
                Console.WriteLine(SupportingFunctions.GetNorm(Ni, sp.x));
                delta = SupportingFunctions.GetNorm(Ni, SupportingFunctions.SubstractVectors(Ni, sp.x, lu.x)) / SupportingFunctions.GetNorm(Ni, sp.x);

                time /= numberOfTests;
                Console.WriteLine("N = {0}, time of LU-decomposition (sec) = {1}, delta = {2}", Ni, Convert.ToDouble(time) / 1000, delta);
            }
        }

    }
    class StartParameters
    {
        public double[] x;
        public double[] A;
        public double[] f;
        public static double[] GenerateX(int Ni)
        {
            double[] x;
            x = new double[Ni];

            for (int i = 0; i < Ni; i++)
            {
                x[i] = 1;
            }

            return x;
        }
        public static double[] GenerateA(int Ni)
        {
            double[] A;
            A = new double[Ni * Ni];

            for (int i = 0; i < Ni; i++)
            {
                for (int j = 0; j < Ni; j++)
                {
                    if (i != j)
                        A[i * Ni + j] = 1 + 0.5 * i - 0.7 * j;
                    else
                        A[i * Ni + j] = 100;

                    //TEST A
                    //A[i * Ni + j] = (i + 1) * (j + 1) + 1;
                }
            }

            return A;
        }
        public static double[] GenerateAFromFile(int Ni)
        {
            double[] A = new double[Ni * Ni];
            StreamReader sr = new("J:\\С#\\source\\repos\\5_nm_1\\data.txt");
            string line;
            string[] str;
            int i = 0;

            line = sr.ReadLine();
            while (line != null)
            {
                str = line.Split(' ');

                for (int j = 0; j < str.Length; j++)
                {
                    A[i * Ni + j] = Convert.ToDouble(str[j]);
                }
                i++;

                line = sr.ReadLine();
            }
            sr.Close();

            return A;
        }
        public static double[] GenerateF(int Ni, double[] A)
        {
            double[] f;
            f = new double[Ni];

            for (int i = 0; i < Ni; i++)
            {
                for (int j = 0; j < Ni; j++)
                {
                    f[i] += A[i * Ni + j];
                }
            }

            return f;
        }
    }
    class LU_Decomposition
    {
        public double[] U;
        public double[] L;
        public double[] y;
        public double[] x;

        public static double[] GenerateU(int Ni, double[] A)
        {
            double temp;
            double[] U;

            U = new double[Ni * Ni];

            for (int i = 0; i < Ni; i++)
            {
                for (int j = 0; j < Ni; j++)
                {
                    U[i * Ni + j] = A[i * Ni + j];
                }
            }

            for (int i = 1; i < Ni; i++)
            {
                //for (int k = i; k < Ni; k++)
                //{
                //    if (U[(i - 1) * Ni + (i - 1)] < U[k * Ni + (i - 1)])
                //    {
                //        for (int column = i - 1; column < Ni; column++)
                //        {
                //            temp = U[(i - 1) * Ni + column];
                //            U[(i - 1) * Ni + column] = U[k * Ni + column];
                //            U[k * Ni + column] = temp;
                //        }
                //    }
                //}
                for (int row = i; row < Ni; row++)
                {
                    temp = U[row * Ni + (i - 1)] / U[(i - 1) * Ni + (i - 1)];
                    temp *= -1;

                    for (int j = 0; j < Ni; j++)
                    {
                        U[row * Ni + j] += U[(i - 1) * Ni + j] * temp;
                    }
                }
            }

            return U;
        }
        public static double[] GenerateGaussResult(int Ni, double[] A, double[] f)
        {
            double temp;
            double[] U = new double[Ni * Ni], x = new double[Ni];


            for (int i = 0; i < Ni; i++)
            {
                for (int j = 0; j < Ni; j++)
                {
                    U[i * Ni + j] = A[i * Ni + j];
                }
            }

            for (int i = 1; i < Ni; i++)
            {
                for (int k = i; k < Ni; k++)
                {
                    if (U[(i - 1) * Ni + (i - 1)] < U[k * Ni + (i - 1)])
                    {
                        for (int column = i - 1; column < Ni; column++)
                        {
                            temp = U[(i - 1) * Ni + column];
                            U[(i - 1) * Ni + column] = U[k * Ni + column];
                            U[k * Ni + column] = temp;
                        }
                        temp = f[(i - 1)];
                        f[(i - 1)] = f[k];
                        f[k] = temp;
                    }
                }
                for (int row = i; row < Ni; row++)
                {
                    temp = U[row * Ni + (i - 1)] / U[(i - 1) * Ni + (i - 1)];
                    temp *= -1;

                    for (int j = 0; j < Ni; j++)
                    {
                        U[row * Ni + j] += U[(i - 1) * Ni + j] * temp;
                    }
                    f[row] += f[i - 1] * temp;
                }
            }

            for (int i = Ni - 1; i >= 0; i--)
            {
                x[i] = f[i];
                for (int j = Ni - 1; j > i; j--)
                {
                    x[i] -= U[i * Ni + j] * x[j];
                }
                x[i] /= U[i * Ni + i];
            }

            return x;
        }
        public static double[] GenerateL(int Ni, double[] A, double[] U)
        {
            double[] L;
            L = new double[Ni * Ni];

            for(int i = 0; i < Ni; i++)
            {
                for (int j = 0; j < Ni; j++)
                {
                    if (i == j)
                    {
                        L[i * Ni + j] = 1;
                    }
                    else if (j > i)
                    {
                        break;
                    }
                    else
                    {
                        L[i * Ni + j] = A[i * Ni + j];
                        for (int k = 0; k < j; k++)
                        {
                            L[i * Ni + j] -= L[i * Ni + k] * U[k * Ni + j];
                        }
                        L[i * Ni + j] /= U[j * Ni + j];
                    }
                }
            }

            return L;
        }
        public static double[] GenerateY(int Ni, double[] L, double[] f)
        {
            double[] y;

            y = new double[Ni];

            for (int i = 0; i < Ni; i++)
            {
                y[i] = f[i];
                for (int j = 0; j < i; j++)
                {
                    y[i] -= L[i * Ni + j] * y[j];
                }
                y[i] /= L[i * Ni + i];
            }

            return y;
        }

        public static double[] GenerateResult(int Ni, double[] U, double[] y)
        {
            double[] x;

            x = new double[Ni];

            for (int i = Ni - 1; i >= 0; i--)
            {
                x[i] = y[i];
                for (int j = Ni - 1; j > i; j--)
                {
                    x[i] -= U[i * Ni + j] * x[j];
                }
                x[i] /= U[i * Ni + i];
            }

            return x;
        }
    }

    class SupportingFunctions
    {
        public static double[] SubstractVectors(int Ni, double[] x, double[] y)
        {
            double[] z = new double[Ni];

            for (int i = 0; i < Ni; i++)
            {
                z[i] = x[i] - y[i];
            }
            return z;
        }

        public static double GetNorm(int Ni, double[] x)
        {
            double norm = 0;

            for (int i = 0; i < Ni; i++)
            {
                norm += x[i] * x[i];
            }
            norm = Math.Pow(norm, 0.5);

            return norm;
        }
    }
}