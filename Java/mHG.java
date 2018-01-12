package statistics;

import java.io.IOException;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

public class mHG {
	
	private static int threshold = 1000;  //The maximal position in the motif binary vector we calculate the p-value for
	private final static int R_ZONE = 0;
	
	public static HGScore calculate_HGT(int N, int n, int B, int b)
	{
		if (B>N || n>N || b>B || b>n) throw new IllegalStateException("Wrong Parameter Sizes!: N="+N+" n="+n+" B="+B+" b="+b);
		double hg = mHG.HG(N, n, B, b, "upper");
		double hgt = mHG.HGT(hg, n, N, B, b);
		HGScore hgScore = new HGScore(N,B,n,b,hgt);
		return hgScore;
	}
	
	/**
	 * based on Eran's code. Computes the hyper geometric tail.
	 * @param currHG The current HG(N,K,n,k) score.
	 * @param N The number of elements in vector
	 * @param B The number of 1's in vector
	 * @param n The number of elements so far
	 * @param b The number of 1's in the small vector that goes until n.
	 * @return The HG tail score
	 */
	public static double HGT(double currHG, int n,int N,int B, int b)
	{
		if(currHG == 0)
		{
			double enrichment1 = (double)b/(double)n;
			double enrichment2 = (double)B/(double)N;
			if(enrichment1 < enrichment2) return 1;
		}
		int min_nB=(n<B) ? n : B;
		double tail=currHG;
		
		for(int i=b;i<min_nB;i++)
		{
			currHG = currHG * ((n-i)*(B-i)) / ((i+1)*(N-n-B+i+1));
			tail += currHG;
		}
		return tail;
	}

	/**
	 * based on Eran's code. Computes the p-value on a mHG score received
	 * in a B set of N ranked list
	 * @param mHGT The min HG score
	 * @param N The number of elements in vector
	 * @param B The number of 1's in vector
	 * @return The HG tail score
	 */
	public static double calculate_pValue(int N, int B, double mHGT)
	{
		return calculate_pValue(N, B, mHGT, threshold);
	}
	public static double calculate_pValue(int N, int B, double mHGT, int threshold)
	{
		int min_NTHRESHOLD = N < threshold ? N : threshold;
		int min_B_nthresh= B < threshold ? B : threshold; 
		if (B==0) return 1;
		
		double[][] mat = new double[min_B_nthresh+1][min_NTHRESHOLD+1];
		for (int i = 0; i <= min_B_nthresh; i++)
			for (int j = 0; j <= min_NTHRESHOLD; j++)
				mat[i][j] = 0;
		
		mat[0][0] = 1;
		double baseHG = 1;	// holds HG(N,K,n,min(n,K))

		for (int n = 1; n <= min_NTHRESHOLD; n++){
			// n is the number of elemnets in current vector
			int min_nB;
			if (B >= n){
				min_nB = n;
				baseHG = baseHG * (B - n + 1) / (N - n + 1);
			}else{
				min_nB = B;
				baseHG = baseHG * n / (n - B);
			}
			
			if (baseHG <= Double.MIN_VALUE){
				baseHG = Double.MIN_VALUE;
				min_NTHRESHOLD=n;
			}

			double tailHG = baseHG;
			double currHG = baseHG;
			// first loop - sum up the tail, until the sum is bigger than mHGT
			int b;
			for (b = min_nB; tailHG <= mHGT && b > 0; b--){
				// b is the number of ones in current vector
				currHG = currHG * (b*(N-B-n+b)) / ((n-b+1)*(B-b+1));
				//if (currHG == 0) currHG = Double.MIN_VALUE;///
				tailHG += currHG;
				mat[b][n] = R_ZONE;
			}
			// second loop, starts when b is the maximal for which
			// HGT(N,B,n,b)> mHGT
			for (; b > 0; b--){
				// calculate current cell value by two optional cells from
				// which it can be reached
				// 1. last element in vector is 0
				
				mat[b][n]=0; //////////////////////// for priniting reasons
				if (mat[b][n-1]<=1){ //////////////////////// for priniting reasons
					mat[b][n] += mat[b][n-1] * (N-B-n+b+1) / (N-n+1);
				}
				// 2. last element in vector is 1
				if (mat[b-1][n-1]<=1){ //////////////////////// for priniting reasons
					mat[b][n] += mat[b-1][n-1] * (B-b+1) / (N-n+1);
				}
				//if (mat[b][n] == 0) mat[b][n] = Double.MIN_VALUE;///
			}
			mat[0][n] = mat[0][n-1] * (N-B-n+1) / (N-n+1);
			
			if (mat[0][n] == Double.MIN_VALUE){
				min_NTHRESHOLD=n;
				//System.err.println("2: n = "+n);
			}
		}

		double result = 0;
		for (int i = 0; i <= min_B_nthresh; i++)
			result += mat[i][min_NTHRESHOLD];
		return (1 - result);
	}
	
	public static double logFactorial(int n, String mode)
	{
		if(n == 0)
		{
			return 0;
		}
		if(n < 0)
		{
			throw new IllegalArgumentException("n should be positive, but it is equal to " + n);
		}
		if(n > 20)
		{
			return logStirlingApproximation(n, mode);
		}
		else
		{
			return Math.log(n) + logFactorial(n-1, mode);
		}
	}
	
	public static double logFactorialExact(int n)
	{
		if(n == 0)
		{
			return 0;
		}
		return Math.log(n) + logFactorialExact(n-1);
	}
	
	public static double logStirlingApproximation(int n, String mode)
	{
		double addition = mode.equalsIgnoreCase("upper") ? 1.0/(12*n) : 1.0/(12*n+1);
		return 0.5*Math.log(2*Math.PI*n) + n*Math.log(n/Math.E) + addition;
	}
		
	public static double log_nCk(int n, int k, String mode)
	{
		String invMode = mode.equalsIgnoreCase("upper") ? "lower" : "upper";
		return logFactorial(n, mode) - logFactorial(k, invMode) - logFactorial(n-k, invMode);
	}
	
	public static double HG(int N, int n, int B, int b, String mode)
	{
		String invMode = mode.equalsIgnoreCase("upper") ? "lower" : "upper";
		double log = log_nCk(n,b, mode) + log_nCk(N-n, B-b, mode) - log_nCk(N,B, invMode);
		return Math.exp(log);
	}
	
//	public static HGScore calculateMinimalHG(int N, SortedSet<Integer> indices)
//	{
//		int B = indices.size();
//		int n = 0;
//		int b = 1;
//		HGScore minHG = new HGScore(N,B,0,0,1);
//		Iterator<Integer> iter = indices.iterator();
//		int prevn = -1;
//		while(iter.hasNext())
//		{
//			n = iter.next();
//			if(n==prevn) System.err.println("double n here");
//			HGScore s = mHG.calculate_HGT(N, n, B, b);
//			if(s.score < minHG.score)
//			{
//				minHG = s;
//			}
//			b++;
//		}
//		return minHG;
//	}
	
	public static HGScore calculateMinimalHGIfGood(int N, SortedSet<Integer> indices, double threshold)
	{
		int B = indices.size();
		double r = (double)B/(double)N;
		HGScore minHG = new HGScore(N,B,0,0,1);
		int n = 0;
		int b = 1;
		Iterator<Integer> iter = indices.iterator();
		while(iter.hasNext())
		{
			n = iter.next();
			if((double)b/(double)n > r)
			{
				double hg = mHG.HG(N, n, B, b, "upper");
				if(hg <= threshold)
				{
					double hgt = mHG.HGT(hg, n, N, B, b);
					if(hgt <= threshold && hgt < minHG.score)
					{
						minHG = new HGScore(N,B,n,b,hgt);
					}
				}
			}
			b++;
		}
		if(minHG.B > 0)
		{
			minHG.pValue = minHG.score*minHG.B;
		}
		if(minHG.pValue > 1 || minHG.B == 0)
		{
			minHG.pValue = 1;
		}
		return minHG;
	}
	
	public static HGScore fastScore(int N, SortedSet<Integer> indices, int i)
	{
		HGScore min = new HGScore(N,0,0,0,1);
		int B = indices.size();
		for(double j=Math.pow(2,i); j>=2; j/=2)
		{
			int n = (int)((double)N/j);
			int b = indices.headSet(n+1).size();
			HGScore s = calculate_HGT(N, n, B, b);
			if(s.score < min.score)
			{
				min = s;
			}
		}
		if(min.B > 0)
		{
			min.pValue = min.score*min.B;
		}
		return min;
	}
	
	public static HGScore calculateHGT(int N, SortedSet<Integer> indices, int n)
	{
		int B = indices.size();
		int i = 0;
		int b = 0;
		Iterator<Integer> iter = indices.iterator();
		while(iter.hasNext())
		{
			i = iter.next();
			if(i <= n) b++;
			else break;
		}
		HGScore s = calculate_HGT(N, n, B, b);
		s.pValue = s.score;
		return s;
	}
	
	public static double boundMinimalHG(int N, SortedSet<Integer> indices)
	{
		double min = 1;
		SortedSet<Integer> truncated = new TreeSet<Integer>();
		for(int i : indices)
		{
			truncated.add(i);
			int B = truncated.size();
			double hg = HG(N, i, B, B, "upper");
			if(hg < min)
			{
				min = hg;
			}
		}
		return min;
	}
	
	//The minimal HG score can be obtained where b/n < B/N, but then it will not be better than 0.5, so testing other threshold makes sense
	public static boolean checkVector(int N, SortedSet<Integer> indices, double threshold)
	{
		int B = indices.size();
		double r = (double)B/(double)N;
		int n = 0;
		int b = 1;
		Iterator<Integer> iter = indices.iterator();
		while(iter.hasNext())
		{
			n = iter.next();
			if((double)b/(double)n > r)
			{
				double hg = mHG.HG(N, n, B, b, "upper");
				if(hg < threshold)
				{
					double hgt = mHG.HGT(hg, n, N, B, b);
					if(hgt <= threshold)
					{
						return true;
					}
				}
			}
			b++;
		}
		return false;
	}
	
	public static void main(String[] args)throws IOException
	{
//		for(int n=1; n<100; n++)
//		{
//			double exact = logFactorialExact(n);
//			double lower = logStirlingApproximation(n,"lower");
//			double upper = logStirlingApproximation(n, "upper");
//			if(lower > exact | upper < exact)
//			{
//				System.out.println("error for n = " + n + " lower = " + lower + " exact = " + exact + " upper " + upper);
//			}
//		}
//		calculate_HGT(int N, int n, int B, int b)(":-)");
		System.out.println(calculate_HGT(107000, 5966, 46278, 442));
//		int N = 100; int B = 4; int n = 14; int b = 2;
//		System.out.println(logFactorial(15));
//		System.out.println(logStirlingApproximation(15));
//		System.out.println("new " + calculate_HGT(N, n, B, b));
//		System.out.println("old " + HG1(N, n, B, b));
	}
}


