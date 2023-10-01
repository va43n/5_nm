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
            int[] N = { 200 };
            //int[] N = { 200, 400, 800 };

            long time;
            int numberOfTests = 1;

            //double[] check = { 0 };

            double delta;

            StartParameters sp = new();
            LU_Decomposition lu = new();
            Householder_Transformation hh = new();

            Stopwatch stopwatch = new();

            foreach (int Ni in N) 
            {
                //-----------------------START PARAMETERS-----------------------

                sp.x = StartParameters.GenerateX(Ni);

                //sp.A = StartParameters.GenerateAFromFile(Ni);
                sp.A = StartParameters.GenerateA(Ni);

                sp.f = StartParameters.GenerateF(Ni, sp.A);

                //-----------------------START PARAMETERS-----------------------

                time = 0;
                for (int i = 0; i < numberOfTests; i++)
                {
                    stopwatch.Restart();

                    lu.U = LU_Decomposition.GenerateU(Ni, sp.A);

                    lu.L = LU_Decomposition.GenerateL(Ni, sp.A, lu.U);

                    lu.y = LU_Decomposition.GenerateY(Ni, lu.L, sp.f);

                    lu.x = LU_Decomposition.GenerateResult(Ni, lu.U, lu.y);

                    stopwatch.Stop();

                    time += stopwatch.ElapsedMilliseconds;
                }
                delta = SupportingFunctions.GetNorm(Ni, SupportingFunctions.SubstractVectors(Ni, sp.x, lu.x)) / SupportingFunctions.GetNorm(Ni, sp.x);

                time /= numberOfTests;
                Console.WriteLine("N = {0}, time of LU-decomposition (sec) = {1}, delta = {2}", Ni, Convert.ToDouble(time) / 1000, delta);


                time = 0;
                for (int i = 0; i < numberOfTests; i++)
                {
                    stopwatch.Restart();

                    hh.GenerateQR(Ni, sp.A);

                    hh.GenerateY(Ni, sp.f);

                    hh.GenerateX(Ni);

                    stopwatch.Stop();

                    time += stopwatch.ElapsedMilliseconds;
                }
                delta = SupportingFunctions.GetNorm(Ni, SupportingFunctions.SubstractVectors(Ni, sp.x, hh.x)) / SupportingFunctions.GetNorm(Ni, sp.x);

                time /= numberOfTests;
                Console.WriteLine("N = {0}, time of Householder transformation (sec) = {1}, delta = {2}", Ni, Convert.ToDouble(time) / 1000, delta);

                //for (int i = 0; i < Ni; i++)
                //{
                //    Console.Write("{0} ", hh.x[i]);
                //}
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

    class Householder_Transformation
    {
        private double[] Q;
        private double[] R;
        private double[] y;
        public double[] x;

        public void GenerateQR(int Ni, double[] A)
        {
            double beta;
            double[] p = new double[1];
            double[] H = new double[Ni * Ni];

            Q = new double[Ni * Ni];
            R = new double[Ni * Ni];

            for (int i = 0; i < Ni * Ni; i++)
            {
                R[i] = A[i];
            }

            for (int i = 0; i < Ni; i++)
            {
                Q[Ni * i + i] = 1;
            }

            for (int column = 0; column < Ni - 1; column++)
            {
                if (-R[Ni * column + column] < 0)
                    beta = -1;
                else
                    beta = 1;

                Array.Resize(ref p, Ni - column);

                for (int i = 0; i < Ni; i++)
                {
                    for (int j = 0; j < Ni; j++)
                    {
                        if (i == j)
                        {
                            H[Ni * i + j] = 1;
                        }
                        else
                        {
                            H[Ni * i + j] = 0;
                        }
                    }
                }

                for (int row = 0; row < Ni - column; row++)
                {
                    p[row] = R[Ni * (row + column) + column];
                }

                beta *= SupportingFunctions.GetNorm(Ni - column, p);
                p[0] -= beta;

                beta = SupportingFunctions.GetNorm(Ni - column, p);
                for (int i = 0; i < Ni - column; i++)
                {
                    p[i] /= beta;
                }

                for (int i = 0; i < Ni - column; i++)
                {
                    for (int j = 0; j < Ni - column; j++)
                    {
                        H[(i + column) * Ni + (j + column)] -= 2 * p[i] * p[j];
                    }
                }

                
                Q = SupportingFunctions.MultiplyMatriсes(Ni, Q, H);
                R = SupportingFunctions.MultiplyMatriсes(Ni, H, R);
            }
        }

        public void GenerateY (int Ni, double[] f)
        {
            y = new double[Ni];

            Q = SupportingFunctions.TransposeMatrix(Ni, Q);

            for (int i = 0; i < Ni; i++)
            {
                for (int j = 0; j < Ni; j++)
                {
                    y[i] += Q[i * Ni + j] * f[j];
                }
            }
        }

        public void GenerateX (int Ni)
        {
            x = new double[Ni];
            for (int i = Ni - 1; i >= 0; i--)
            {
                x[i] = y[i];
                for (int j = Ni - 1; j > i; j--)
                {
                    x[i] -= R[i * Ni + j] * x[j];
                }
                x[i] /= R[i * Ni + i];
            }
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
        
        public static double[] MultiplyMatriсes(int Ni, double[] A, double[] B)
        {
            double[] C = new double[Ni * Ni];

            for (int i = 0; i < Ni; i++)
            {
                for (int j = 0; j < Ni; j++)
                {
                    for (int k = 0; k < Ni; k++)
                    {
                        C[i * Ni + j] += A[i * Ni + k] * B[k * Ni + j];
                    }
                }
            }

            return C;
        }

        public static double[] TransposeMatrix(int Ni, double[] A)
        {
            double temp;

            for (int i = 0; i < Ni; i++)
            {
                for (int j = i; j < Ni; j++)
                {
                    temp = A[i * Ni + j];
                    A[i * Ni + j] = A[j * Ni + i];
                    A[j * Ni + i] = temp;
                }
            }

            return A;
        }
    }
}