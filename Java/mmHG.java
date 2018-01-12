package statistics;
import java.util.HashMap;
import java.util.Map;

import statistics.HGScore;
import statistics.mHG;

public class mmHG {
	
	static final int MAX_THRESHOLD = 3000;
	
	/**
	 * This method takes as input a permutation (and possibly a cutoff) and returns the mmHG score and the mmHG p-value.
	 * @param perm - a permutation over the numbers 1,...,N
	 * @param nstar - indicates whether the list should be cut at a specific position. Set to -1 if no such cutoff exists.
	 * @returns the mmHG score which is the minimal hyper-geometric tail found. The p-value field contains the Bonferroni correction. 
	 */
	public static HGScore calcScore(int[] perm, int nstar)
	{
		int N = perm.length;
		Map<Integer, Integer> indMapper = new HashMap<Integer, Integer>();
		int M = N < MAX_THRESHOLD ? N : MAX_THRESHOLD;
		for(int k=0; k<N; k++)
		{
			indMapper.put(perm[k], k);
		}
		HGScore min = new HGScore(N,0,0,0,1);
		int maxVal = 0;
		for(int n2=1; n2<M; n2++) //n2 represents the position in the permutation
		{
			if(perm[n2] > maxVal)
			{
				maxVal = perm[n2];
			}
			int b = 0;
			if(nstar == -1)
			{
				for(int n1=1; n1<=maxVal && n1<=M; n1++) //n1 represents the range 1...n1
				{
					if(indMapper.get(n1) < n2)
					{
						b++;
						HGScore s = calcHgtScore(N, n1, n2, b);
						if(s.score < min.score)
						{
							min = s;
						}
					}
				}
			}
			else
			{
				for(int n1=1; n1<=nstar; n1++)
				{
					if(indMapper.get(n1) < n2)
					{
						b++;
					}
				}
				HGScore s = calcHgtScore(N, nstar, n2, b);
				if(s.score < min.score)
				{
					min = s;
				}
			}
		}
		min.pValue = min.score*M*M; //Bonferroni's correction for the p-value
		if(min.pValue > 1)
		{
			min.pValue = 1;
		}
		if(min.score >= 0 && min.score < Double.MIN_VALUE) //handle cases where score=0 due to double precision 
		{
			min.score = Double.MIN_VALUE;
			min.pValue = Double.MIN_VALUE;
		}
		return min;
	}
	
	public static HGScore calcHgtScore(int N, int n1, int n2, int b)
	{
		HGScore s = mHG.calculate_HGT(N, n1, n2, b);
		return s;
	}
	
	public static void main(String[] args) 
	{
		int[] perm = {2,5,1,4,3};
		System.out.println(calcScore(perm, -1));
	}

}
