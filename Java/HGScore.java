package statistics;

public class HGScore {
	
	public int N,B,n,b;
	public double score;
	public double pValue = 1;
	
	public HGScore(int N, int B, int n, int b, double score){
		this.N=N;this.n=n;this.B=B;this.b=b;
		this.score=score;
	}
	
	public String toString(){ 
		String rep = "N "+N+", B "+B+", n "+n+", b "+b+", score "+score;
		if (pValue<1)
			rep += ", p-value "+pValue;
		return rep;
	}
	
	public double calcPvalue(int threshold)
	{
//		if(score < 1e-15/*~10^-7*/) 
//			pValue = score*B;
//		else
		{
			pValue = mHG.calculate_pValue(N, B, score, threshold);
			if (pValue < score)
			{
				pValue = score*B;
			}
			if(pValue > score*B) //bug correction
			{
				pValue = score*B;
			}
		}
		
		return Math.min(pValue,1);
	}
}