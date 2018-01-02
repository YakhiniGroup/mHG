using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace mHGJumper
{
    public enum mHGCorrectionType { Exact, Lipson, Bonferroni, None }
    /// <summary>
    /// This class computes mHG efficiently for traversal of the spatial enrichment space.
    /// It can also offer precomputed jump distances given an mHG score to the nearest possible neighbor
    /// with a better score.
    /// </summary>
    public static class mHGJumper
    {
        /// <summary>
        /// Algorithm details available here http://www.ploscompbiol.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pcbi.0030039&representation=PDF
        /// </summary>
        /// <param name="binVec">relevant vector subset (topK from original)</param>
        /// <param name="tN">universe size (length of original vector)</param>
        /// <param name="tB">|falses|</param>
        /// <returns>mHG p-Value and threshold</returns>
        private const int c_rZone = 0;
        private const double TOLERANCE = double.Epsilon;
        public static readonly double Epsilon = 0.01;
        private static double[,] HGTmat;
        
        private static object concLocker = new object();
        public static double? optHGT = 0.05;

        private static double TotalPaths;
        private static int Ones, Zeros;
        private static Dictionary<double, double> ScoreMap;
        private static Random rnd = new Random();
        
        public static void Initialize(int ones, int zeros)
        {
            optHGT = 0.05;
            if (Ones == ones && Zeros == zeros && HGTmat != null)
            {
                return; //this was pre-initialized sometime
            }

            Ones = ones;
            Zeros = zeros;
            var N = zeros + ones;
            HGTmat = new double[zeros + 1, ones + 1];
            HGTmat[0, 0] = 1.0; //Base condition

            //init zeros
            for (var i = 1; i < zeros + 1; i++)
                HGTmat[i, 0] = HGfromPrevious(HGTmat[i - 1, 0], N, ones, i - 1, 0, false);
            
            //init ones
            for (var j = 1; j < ones + 1; j++)
                HGTmat[0, j] = HGfromPrevious(HGTmat[0, j - 1], N, ones, j - 1, j - 1, true);
            
            //compute HG w/ dynamic program
            for (var i = 1; i < zeros + 1; i++)
            {
                HGTmat[i, 1] = HGfromPrevious(HGTmat[i - 1, 1], N, ones, i, 1, false);
                for (var j = 2; j < ones + 1; j++)
                {
                    HGTmat[i, j] = HGfromPrevious(HGTmat[i, j - 1], N, ones, i + j - 1, j - 1, true);
                }
            }

            var scoreToPval = new HashSet<double>();
            //Sum for the 'upper' HGT
            for (var diagsum = 0; diagsum <= zeros+ones; diagsum++)
            {
                var runsum = 0.0;
                for (var numones = Math.Min(diagsum,ones); numones >= 0; numones--)
                {
                    var numzeros = diagsum - numones;
                    if (numzeros > zeros)
                        continue;
                    runsum += HGTmat[numzeros, numones];
                    HGTmat[numzeros, numones] = runsum;
                    //if(runsum<0.5)
                    scoreToPval.Add(runsum);
                }
            }
            TotalPaths = pathCounting(-0.1, out var pMat); //Todo: should be N choose B
            var tscoremap = new ConcurrentDictionary<double,double>();
            Console.WriteLine("Mapping {0} mHG scores to pvalue.", scoreToPval.Count);
            Parallel.ForEach(scoreToPval, score =>
            {
                var pval = 1.0 - (pathCounting(score, out var pMat1) / TotalPaths);
                tscoremap.AddOrUpdate(score, pval, (a, b) => pval);
            });
            ScoreMap = tscoremap.ToDictionary(t => t.Key, t => t.Value);
            Console.WriteLine("Done initializing HGT matrix of size {0}x{1}",zeros,ones);
        }

        /// <summary>
        /// Computes the recurrence relation for the hypergeometric pdf.
        /// </summary>
        /// <param name="hgt0">previous hg pdf </param>
        /// <param name="N">#items</param>
        /// <param name="K">Total possible successes</param>
        /// <param name="n0">previous num trials</param>
        /// <param name="k0">previous #success</param>
        /// <param name="isSuccess">Is the current trial a succcess</param>
        /// <returns></returns>
        public static double HGfromPrevious(double hgt0, int N, int K, int n0, int k0, bool isSuccess)
        {
            var hgt1 = hgt0 * ((n0 + 1) / (double) (N - n0)) *
                       (isSuccess ? (K - k0) / (double) (k0 + 1) : (N - K - n0 + k0) / (double) (n0 - k0 + 1));
            return hgt1;
        }


        /// <summary>
        /// Counts with dynamic program the number of paths that dont go through a HGT score.
        /// </summary>
        /// <param name="hgtScore"></param>
        /// <returns></returns>
        public static double pathCounting(double hgtScore, out double[,] pMat)
        {
            pMat = new double[Zeros + 1, Ones + 1];
            //There is exactly one path that travel on edges
            for (var i = 0; i < Zeros + 1; i++)
                pMat[i, 0] = HGTmat[i, 0] > hgtScore ? 1 : 0;
            for (var j = 0; j < Ones + 1; j++)
                pMat[0, j] = HGTmat[0, j] > hgtScore ? 1 : 0;

            for (var i = 1; i < Zeros + 1; i++)
                for (var j = 1; j < Ones + 1; j++)
                {
                    var isInR = HGTmat[i, j] <= hgtScore + TOLERANCE;
                    pMat[i, j] = isInR ? 0 : pMat[i - 1, j] + pMat[i, j - 1];
                }
            return pMat[Zeros,Ones];
        }


        public static Tuple<double, int> minimumHypergeometric(bool[] binVec, int tN = -1, int tB = -1, mHGCorrectionType correctMultiHypothesis = mHGCorrectionType.Exact)
        {
            var N = tN > 0 ? tN : binVec.Length;
            var K = tB > 0 ? tB : binVec.Sum(val => val ? 1 : 0);
            var B = tB > 0 ? tB : binVec.Sum(val => !val ? 1 : 0);
            var currHGT = 1.0;
            var mHGT = 1.1;
            var currIndex = 0;
            var k = 0;
            //OptDistVec is a vector that counts for each '1' in the binary vector the minimum number of 1's needed directly after it for a significant p-value
            var OptDistVec = new int[Ones+1];
            for (var i = 0; i < Ones + 1; i++)
                OptDistVec[i] = Ones; // default max step size (int.maxvalue)
            for (var i = 0; i < Ones + 1; i++) if (HGTmat[1, i] <= optHGT.Value) OptDistVec[0] = Math.Min(OptDistVec[0], i);
            for (var n = 0; n < binVec.Length; n++)
            {
                if (binVec[n])
                {
                    k++;
                    currHGT = HGTmat[n - k + 1, k];
                    //currHGT = ScoreMap[currHG];
                    if (currHGT < mHGT)
                    {
                        currIndex = n;
                        mHGT = currHGT;
                    }
                    //check distance to optimum
                    if (optHGT.HasValue && optHGT <= currHGT)
                    {
                        for (var i = k; i < Ones+1; i++)
                            if (HGTmat[n - k + 1, i] <= optHGT.Value)
                                OptDistVec[k] = Math.Min(OptDistVec[k], i - k);
                    }
                    else
                    {
                        optHGT = currHGT;
                    }
                }
            }
            double pval = -1;
            switch (correctMultiHypothesis)
            {
                case mHGCorrectionType.Exact:
                    pval = ScoreMap[mHGT];
                    break;
                case mHGCorrectionType.None:
                    pval = mHGT;
                    break;
                case mHGCorrectionType.Bonferroni:
                    pval = mHGT * N;
                    break;
                case mHGCorrectionType.Lipson:
                    pval = mHGT * B;
                    break;
            }

            return new Tuple<double, int>(pval, currIndex + 1);
        }
    }
}
