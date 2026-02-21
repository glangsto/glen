import re

def parse_idl_file(filename):
    with open(filename, 'r') as f:
        text = f.read()

    # Remove IDL comments (anything after ;)
    text_no_comments = re.sub(r';.*', '', text)

    # Normalize line continuations ($\)
    text_no_comments = text_no_comments.replace('$\\', '')

    # -------------------------------
    # 1. Declared variables (left of =)
    # -------------------------------
    declared_vars = re.findall(r'^\s*([A-Za-z_]\w*)\s*=', text_no_comments, flags=re.MULTILINE)

    # -------------------------------
    # 2. Procedure/function calls
    #    IDL syntax:  print, x
    #    or function calls: sin(x)
    # -------------------------------
    proc_calls = re.findall(r'^\s*([A-Za-z_]\w*)\s*,', text_no_comments, flags=re.MULTILINE)
    func_calls = re.findall(r'([A-Za-z_]\w*)\s*\(', text_no_comments)

    called_funcs = sorted(set(proc_calls + func_calls))

    # -------------------------------
    # 3. Referenced variables
    #    Any identifier not declared and not a function
    # -------------------------------
    # All identifiers
    all_ids = re.findall(r'\b([A-Za-z_!]\w*)\b', text_no_comments)

    # Remove declared vars and function names
    referenced = [
        v for v in all_ids
        if v not in declared_vars and v not in called_funcs
    ]

    # Remove numeric-like tokens and IDL keywords
    keywords = {"if", "then", "begin", "end", "for", "do", "while", "else"}
    referenced = [v for v in referenced if v.lower() not in keywords]

    # Deduplicate while preserving order
    seen = set()
    referenced_vars = []
    for v in referenced:
        if v not in seen:
            seen.add(v)
            referenced_vars.append(v)

    return declared_vars, called_funcs, referenced_vars


# -------------------------------
# Example usage
# -------------------------------
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        print("Supply the line list file")
        exit()

    print("Parsing file: %s" % (filename))
        
    declared, funcs, referenced = parse_idl_file(filename)

    print("Declared variables:", declared)
    print("Function/procedure calls:", funcs)
    print("Referenced variables:", referenced)
