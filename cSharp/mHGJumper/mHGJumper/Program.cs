using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mHGJumper
{
    class Program
    {
        static void Main(string[] args)
        {
            var rnd = new Random(0);
            var inputVec = Enumerable.Range(0, 100).Select((v,idx) => rnd.NextDouble() > 0.1+idx/10.0).ToArray();
            var ones = inputVec.Count(v => v);
            var zeros = inputVec.Length - ones;
            mHGJumper.Initialize(ones, zeros);
            var result = mHGJumper.minimumHypergeometric(inputVec);
            Console.WriteLine($"Input vector ones:{ones} zeros:{zeros}");
            Console.WriteLine(string.Join(",", inputVec));
            Console.WriteLine($"p-value={result.Item1:E} threshold={result.Item2} ones in threshold={inputVec.Take(result.Item2).Count(v=>v)}");
            Console.WriteLine("Press any key to end.");
            Console.ReadKey();
        }
    }
}
