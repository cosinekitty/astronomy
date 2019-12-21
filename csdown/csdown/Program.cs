using System;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Xml.Linq;

namespace csdown
{
    class Program
    {
        static int Main(string[] args)
        {
            if (args.Length != 4)
            {
                Console.WriteLine("USAGE: csdown prefix.md assembly.dll assembly.xml outfile.md");
                return 1;
            }
            string inPrefixFileName = args[0];
            string inAssemblyFileName = args[1];
            string inXmlFileName = args[2];
            string outMarkdownFileName = args[3];
            try
            {
                return GenerateMarkdown(inPrefixFileName, inAssemblyFileName, inXmlFileName, outMarkdownFileName);
            }
            catch (Exception ex)
            {
                Console.WriteLine("EXCEPTION(csdown): {0}", ex);
                return 1;
            }
        }

        static int GenerateMarkdown(string inPrefixFileName, string inAssemblyFileName, string inXmlFileName, string outMarkdownFileName)
        {
            Assembly asm = Assembly.LoadFile(Path.GetFullPath(inAssemblyFileName));
            var sb = new StringBuilder(File.ReadAllText(inPrefixFileName));
            XDocument doc = XDocument.Load(inXmlFileName);
            var cinfo = new CodeInfo(doc);
            AppendMarkdown(sb, cinfo, asm);
            File.WriteAllText(outMarkdownFileName, sb.ToString());
            return 0;
        }

        private static void AppendMarkdown(StringBuilder sb, CodeInfo cinfo, Assembly asm)
        {
            Console.WriteLine("Generating C# documentation.");

            sb.AppendLine("<a name=\"functions\"></a>");
            sb.AppendLine("## Functions");
            sb.AppendLine();
            sb.AppendLine("---");
            sb.AppendLine();

            Type astro = asm.GetType("CosineKitty.Astronomy");
            MethodInfo[] funcs = astro.GetMethods()
                .Where(m => m.IsPublic && m.IsStatic)
                .OrderBy(m => m.Name.ToUpperInvariant())
                .ToArray();

            foreach (MethodInfo f in funcs)
            {
                CodeItem item = cinfo.FindMethod(f);

                sb.AppendFormat("<a name=\"Astronomy.{0}\"></a>", f.Name);
                sb.AppendLine();
                sb.AppendFormat("### Astronomy.{0}(", f.Name);
                sb.Append(string.Join(", ", f.GetParameters().Select(p => p.Name)));
                sb.AppendFormat(") &#8658; {0}", TypeMarkdown(f.ReturnType));
                sb.AppendLine();
                sb.AppendLine();
                if (!string.IsNullOrWhiteSpace(item.Summary))
                {
                    sb.AppendLine("**" + item.Summary + "**");
                    sb.AppendLine();
                }

                if (!string.IsNullOrWhiteSpace(item.Remarks))
                {
                    sb.AppendLine(item.Remarks);
                    sb.AppendLine();
                }
            }
        }

        private static string TypeMarkdown(Type t)
        {
            if (t.FullName.StartsWith("CosineKitty."))
                return string.Format("[`{0}`](#{0})", t.Name);

            switch (t.FullName)
            {
                case "System.Double":
                    return "double";

                case "System.String":
                    return "string";

                default:
                    throw new NotImplementedException("Unhandled type " + t.FullName);
            }
        }
    }
}
