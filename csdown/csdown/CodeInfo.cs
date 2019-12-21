using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml.Linq;

namespace csdown
{
    class CodeItem
    {
        public string Summary;
        public string Remarks;
        public string Returns;
        public readonly Dictionary<string, string> Params = new Dictionary<string, string>();
    }

    class CodeInfo
    {
        private readonly Dictionary<string, CodeItem> table = new Dictionary<string, CodeItem>();

        public CodeInfo(XDocument doc)
        {
            XElement members = doc.Root.Element("members");
            foreach (XElement child in members.Elements("member"))
            {
                if (child.Attribute("name") is XAttribute idAttr && idAttr.Value is string id && id.Length > 2 && id[1] == ':')
                {
                    string name = id.Substring(2);

                    var item = new CodeItem();
                    table.Add(name, item);

                    if (child.Element("summary") is XElement summary)
                        item.Summary = Clean(summary.Value);

                    if (child.Element("remarks") is XElement remarks)
                        item.Remarks = Clean(remarks.Value);

                    if (child.Element("returns") is XElement returns)
                        item.Returns = Linear(returns.Value);

                    foreach (XElement p in child.Elements("param"))
                        if (p.Attribute("name") is XAttribute pname)
                            item.Params[pname.Value] = Linear(p.Value);
                }
            }
        }

        private static int LeftSpaceCount(string line)
        {
            if (line == null)
                return 0;

            for (int i = 0; i < line.Length; ++i)
                if (line[i] != ' ')
                    return i;

            return line.Length;
        }

        internal static string Linear(string text)
        {
            if (text == null)
                return "";
            return Regex.Replace(text.Trim(), @"\s+", " ");
        }

        private static string Clean(string text)
        {
            // Find the common amount of space on the left of each line.
            // Remove that amount of space from all lines, but no more.
            string[] lines = text.Split('\n');

            // Remove any blank leading lines and any blank trailing lines.
            // Leave blank interior lines intact.
            int firstNonBlank = -1;
            int lastNonBlank = -1;
            for (int i = 0; i < lines.Length; ++i)
            {
                if (!string.IsNullOrWhiteSpace(lines[i]))
                {
                    lastNonBlank = i;
                    if (firstNonBlank < 0)
                        firstNonBlank = i;
                }
            }

            if (firstNonBlank < 0)
                return "";

            lines = lines.Skip(firstNonBlank).Take(lastNonBlank - firstNonBlank + 1).ToArray();

            if (lines.Length > 0)
            {
                int lspace = -1;
                foreach (string s in lines)
                {
                    if (!string.IsNullOrWhiteSpace(s))
                    {
                        int left = LeftSpaceCount(s);
                        if (lspace < 0 || left < lspace)
                            lspace = left;
                    }
                }
                if (lspace > 0)
                {
                    for (int i = 0; i < lines.Length; ++i)
                        if (string.IsNullOrWhiteSpace(lines[i]))
                            lines[i] = "";
                        else
                            lines[i] = lines[i].Substring(lspace);
                }
                return string.Join("\n", lines);
            }
            return "";
        }

        public CodeItem FindMethod(MethodInfo f)
        {
            string id = f.DeclaringType.FullName + "." + f.Name;
            ParameterInfo[] parms = f.GetParameters();
            if (parms.Length > 0)
                id += "(" + string.Join(",", parms.Select(p => p.ParameterType.FullName)) + ")";

            return table[id];
        }

        internal CodeItem FindEnumValue(FieldInfo f)
        {
            string id = f.DeclaringType.FullName + "." + f.Name;
            return table[id];
        }
    }
}
