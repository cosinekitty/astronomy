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
            return GenerateMarkdown(inPrefixFileName, inAssemblyFileName, inXmlFileName, outMarkdownFileName);
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

            Type astro = asm.GetType("CosineKitty.Astronomy");

            // Document all public constants in the Astronomy class.
            FieldInfo[] constants = astro.GetFields(BindingFlags.Public | BindingFlags.Static)
                .Where(f => f.IsLiteral && !f.IsInitOnly)
                .ToArray();

            AppendSectionHeader(sb, "Constants");
            foreach (FieldInfo c in constants)
                AppendConstantMarkdown(sb, cinfo, c);

            // All the member functions in the Astronomy class go in the "functions" section.

            AppendSectionHeader(sb, "Functions");
            MethodInfo[] funcs = astro.GetMethods()
                .Where(m => m.IsPublic && m.IsStatic)
                .OrderBy(m => m.Name.ToUpperInvariant())
                .ToArray();

            foreach (MethodInfo f in funcs)
                AppendFunctionMarkdown(sb, cinfo, f);

            // The classes other than "Astronomy" are listed one by one in the "classes" section.
            sb.AppendLine("---");
            sb.AppendLine();
            AppendSectionHeader(sb, "Types");
            Type[] typeList = asm.GetExportedTypes()
                .Where(c => c.Name != "Astronomy")
                .OrderBy(c => c.Name.ToUpperInvariant())
                .ToArray();

            foreach (Type c in typeList)
                AppendTypeMarkdown(sb, cinfo, c);
        }

        private static void AppendTypeMarkdown(StringBuilder sb, CodeInfo cinfo, Type type)
        {
            string kind;
            if (type.IsClass)
                kind = "class";
            else if (type.IsEnum)
                kind = "enum";
            else
                kind = "struct";

            CodeItem typeItem = cinfo.FindType(type);

            // Header

            sb.AppendLine("<a name=\"" + type.Name + "\"></a>");
            sb.AppendLine("## `" + kind + " " + type.Name + "`");
            sb.AppendLine();

            if (!string.IsNullOrWhiteSpace(typeItem.Summary))
            {
                sb.AppendLine("**" + typeItem.Summary + "**");
                sb.AppendLine();
            }

            if (!string.IsNullOrWhiteSpace(typeItem.Remarks))
            {
                sb.AppendLine(typeItem.Remarks);
                sb.AppendLine();
            }

            if (type.IsEnum)
            {
                // Dump enum values.
                FieldInfo[] fields = type.GetFields().Where(f => f.IsLiteral).ToArray();
                if (fields.Length > 0)
                {
                    sb.AppendLine("| Value | Description |");
                    sb.AppendLine("| --- | --- |");
                    foreach (FieldInfo f in fields)
                        AppendEnumValueMarkdown(sb, cinfo, f);
                    sb.AppendLine();
                }
            }
            else
            {
                // Dump struct/class fields.
                FieldInfo[] fields = type.GetFields();
                if (fields.Length > 0)
                {
                    sb.AppendLine("| Type | Name | Description |");
                    sb.AppendLine("| --- | --- | --- |");
                    foreach (FieldInfo f in fields)
                        if (f.DeclaringType == type)
                            AppendMemberVariableMarkdown(sb, cinfo, f);
                    sb.AppendLine();
                }
            }

            // Member functions

            MethodInfo[] funcs = type.GetMethods()
                .Where(m => m.IsPublic && m.DeclaringType == type)
                .OrderBy(m => m.Name.ToUpperInvariant())
                .ToArray();

            if (funcs.Length > 0)
            {
                sb.AppendLine("### member functions");
                sb.AppendLine();

                foreach (MethodInfo f in funcs)
                    AppendFunctionMarkdown(sb, cinfo, f);
            }

            sb.AppendLine("---");
            sb.AppendLine();
        }

        private static void AppendMemberVariableMarkdown(StringBuilder sb, CodeInfo cinfo, FieldInfo f)
        {
            CodeItem item = cinfo.FindField(f);
            sb.Append("| ");
            sb.Append(TypeMarkdown(f.FieldType));
            sb.Append(" | `");
            sb.Append(f.Name);
            sb.Append("` | ");
            sb.Append(CodeInfo.Linear(item.Summary));
            sb.AppendLine(" |");
        }

        private static void AppendEnumValueMarkdown(StringBuilder sb, CodeInfo cinfo, FieldInfo f)
        {
            CodeItem item = cinfo.FindField(f);
            sb.Append("| `");
            sb.Append(f.Name);
            sb.Append("` | ");
            sb.Append(CodeInfo.Linear(item.Summary));
            sb.AppendLine(" |");
        }

        private static void AppendSectionHeader(StringBuilder sb, string name)
        {
            sb.AppendLine("<a name=\"" + name.ToLowerInvariant() + "\"></a>");
            sb.AppendLine("## " + name);
            sb.AppendLine();
            sb.AppendLine("---");
            sb.AppendLine();
        }

        private static void AppendConstantMarkdown(StringBuilder sb, CodeInfo cinfo, FieldInfo f)
        {
            CodeItem item = cinfo.FindField(f);
            if (item == null)
                return;

            string fieldType;
            switch (f.FieldType.Name)
            {
                case "Double":
                    fieldType = "double";
                    break;

                default:
                    throw new Exception($"Do not know how to generate markdown for constant type: {f.FieldType.Name}");
            }

            string parentClassName = f.DeclaringType.Name;

            sb.AppendFormat("<a name=\"{0}.{1}\"></a>", parentClassName, f.Name);
            sb.AppendLine();
            sb.AppendFormat("### `const {0} {1}.{2} = {3};`", fieldType, parentClassName, f.Name, f.GetValue(null));
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

        private static void AppendFunctionMarkdown(StringBuilder sb, CodeInfo cinfo, MethodInfo f)
        {
            CodeItem item = cinfo.FindMethod(f);
            if (item == null)
                return;

            ParameterInfo[] parms = f.GetParameters();
            string parentClassName = f.DeclaringType.Name;

            sb.AppendFormat("<a name=\"{0}.{1}\"></a>", parentClassName, f.Name);
            sb.AppendLine();
            sb.AppendFormat("### {0}.{1}(", parentClassName, f.Name);
            sb.Append(string.Join(", ", parms.Select(p => p.Name)));
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

            if (parms.Length > 0)
            {
                sb.AppendLine("| Type | Parameter | Description |");
                sb.AppendLine("| --- | --- | --- |");
                foreach (ParameterInfo p in parms)
                {
                    // | [`astro_rotation_t`](#astro_rotation_t) | `a` |  The first rotation to apply. | 
                    string t = TypeMarkdown(p.ParameterType);
                    string n = "`" + p.Name + "`";
                    string r = item.Params[p.Name];
                    sb.Append("| ");
                    sb.Append(t);
                    sb.Append(" | ");
                    sb.Append(n);
                    sb.Append(" | ");
                    sb.Append(r);
                    sb.AppendLine(" |");
                }
                sb.AppendLine();
            }

            if (!string.IsNullOrWhiteSpace(item.Returns))
            {
                sb.AppendLine("**Returns:** " + item.Returns);
                sb.AppendLine();
            }
        }

        private static string TypeMarkdown(Type t)
        {
            if (t.FullName.StartsWith("CosineKitty."))
                return CodeInfo.InternalLink(t.Name);

            switch (t.FullName)
            {
                case "System.Double":
                    return "`double`";

                case "System.Double[,]":
                    return "`double[3,3]`";     // the two-dimensional arrays are all rotation matrices: 3x3

                case "System.String":
                    return "`string`";

                case "System.Int32":
                    return "`int`";

                case "System.DateTime":
                    return "`DateTime`";

                default:
                    throw new NotImplementedException("Unhandled type " + t.FullName);
            }
        }
    }
}
