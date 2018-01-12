using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mmHG
{
    class Program
    {
        public static void Main(String[] args)
        {
            MHG.SetMaxN(10);
            int[] perm = { 1, 4, 0, 3, 2 };
            Console.WriteLine(MMHG.CalcScore(perm));
            Console.ReadKey();
        }
    }
}
