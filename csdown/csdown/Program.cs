using System;

namespace csdown
{
    class Program
    {
        static int Main(string[] args)
        {
            if (args.Length != 3)
            {
                Console.WriteLine("USAGE: csdown prefix.md assembly.dll outfile.md");
                return 1;
            }
            string inPrefixFileName = args[0];
            string inAssemblyFileName = args[1];
            string outMarkdownFileName = args[2];
            try
            {
                return GenerateMarkdown(inPrefixFileName, inAssemblyFileName, outMarkdownFileName);
            }
            catch (Exception ex)
            {
                Console.WriteLine("EXCEPTION(csdown): {0}", ex);
                return 1;
            }
        }

        static int GenerateMarkdown(string inPrefixFileName, string inAssemblyFileName, string outMarkdownFileName)
        {
            Console.WriteLine("GenerateMarkdown");
            return 0;
        }
    }
}
