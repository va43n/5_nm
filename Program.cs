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
            int[] N = { 100 };
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
                //    Console.Write("{0} ", Math.Round(hh.x[i]));
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
            double beta, dtemp, mu;
            double[] p = new double[Ni];
            double[] H = new double[Ni * Ni];
            //double[] temp;

            Q = new double[Ni * Ni];
            R = new double[Ni * Ni];

            //temp = new double[Ni * Ni];

            //for (int i = 0; i < Ni * Ni; i++)
            //{
            //    R[i] = A[i];
            //}

            for (int i = 0; i < Ni; i++)
            {
                Q[Ni * i + i] = 1;
            }

            //for (int column = 0; column < Ni - 1; column++)
            //{
            //    if (-A[Ni * column + column] < 0)
            //        beta = -1;
            //    else
            //        beta = 1;

            //    for (int row = 0; row < Ni; row++)
            //    {
            //        p[row] = A[Ni * row + column];
            //    }

            //    beta *= SupportingFunctions.GetNorm(Ni - column, p);

            //    if (Math.Abs(beta) - A[column * Ni + column] < Math.Pow(10, -15))
            //    {
            //        continue;
            //    }

            //    p[column] -= beta;

            //    mu = 1.0 / beta / (beta - p[column]);

            //    for (int i = 0; i < Ni; i++)
            //    {
            //        for (int j = 0; j < Ni; j++)
            //        {
            //            if (i == j)
            //            {
            //                H[Ni * i + j] = 1;
            //            }
            //            else
            //            {
            //                H[Ni * i + j] = 0;
            //            }
            //        }
            //    }

            //    for(int i = column; i < Ni; i++)
            //    {
            //        dtemp = 0;
            //        for (int j = column; j < Ni; j++)
            //        {
            //            dtemp += A[j * Ni + column] * p[j];
            //        }
            //        dtemp *= mu;
            //        for (int j = column; j < Ni; j++)
            //        {
            //            A[j * Ni + column] -= dtemp * p[j];
            //        }
            //    }

            //    for (int i = 0; i < Ni; i++)
            //    {
            //        dtemp = 0;
            //        for (int j = column; j < Ni; j++)
            //        {
            //            dtemp += Q[j * Ni + column] * p[j];
            //        }
            //        dtemp *= mu;
            //        for (int j = column; j < Ni; j++)
            //        {
            //            Q[j * Ni + column] -= dtemp * p[j];
            //        }
            //    }
            //}
            //for (int i = 0; i < Ni * Ni; i++)
            //{
            //    R[i] = A[i];
            //}


            for (int i = 0; i < Ni * Ni; i++)
            {
                R[i] = A[i];
            }

                //алгоритм отражений Хаусхолдера
            for (int i = 0; i < Ni - 1; i++)
            {
                //находим квадрат нормы столбца для обнуления
                dtemp = 0;
                for (int I = i; I < Ni; I++) dtemp += Math.Pow(R[I * Ni + i], 2);

                //если есть ненулевые элементы под диагональю, тогда 
                //квадрат нормы вектора для обнуления не совпадает с квадратом диагонального элемента
                if (Math.Sqrt(Math.Abs(dtemp - R[i * Ni + i] * R[i * Ni + i])) > Math.Pow(10, -15))
                {
                    //выбор знака слагаемого beta = sign(-x1)
                    if (R[i * Ni + i] < 0) beta = Math.Sqrt(dtemp);
                    else beta = -Math.Sqrt(dtemp);

                    //вычисляем множитель в м.Хаусхолдера mu = 2 / ||p||^2
                    mu = 1.0 / beta / (beta - R[i * Ni + i]);

                    //формируем вектор p
                    for (int I = 0; I < Ni; I++) { p[I] = 0; if (I >= i) p[I] = R[I * Ni + i]; }

                    //изменяем диагональный элемент
                    p[i] -= beta;

                    //вычисляем новые компоненты матрицы A = Hm * H(m-1) ... * A
                    for (int m = i; m < Ni; m++)
                    {
                        //произведение S = At * p
                        dtemp = 0;
                        for (int n = i; n < Ni; n++) { dtemp += R[n * Ni + m] * p[n]; }
                        dtemp *= mu;
                        //A = A - 2 * p * (At * p)^t / ||p||^2
                        for (int n = i; n < Ni; n++) { R[n * Ni + m] -= dtemp * p[n]; }
                    }

                    //вычисляем новые компоненты матрицы Q = Q * H1 * H2 * ...
                    for (int m = 0; m < Ni; m++)
                    {
                        //произведение Q * p
                        dtemp = 0;
                        for (int n = i; n < Ni; n++) { dtemp += Q[m * Ni + n] * p[n]; }
                        dtemp *= mu;
                        //Q = Q - p * (Q * p)^t
                        for (int n = i; n < Ni; n++) { Q[m * Ni + n] -= dtemp * p[n]; }
                    }
                }
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