# parse an file with three arrays
#HISTORY
#26Feb16 GIL add function to get individual variables
#26Feb13 GIL initial version
import sys
import re

def extract_idl_arrays(text):
    """
    Extract the VALUES of the first three IDL array assignments.
    Returns a tuple of three Python lists.
    """

    # Remove comments
    text = re.sub(r';.*', '', text)

    # Remove IDL line continuation "$\"
    text = text.replace('$\\', '')

    # Find array assignments:  var = [ ... ]
    array_blocks = re.findall(
        r'=\s*\[(.*?)\]',
        text,
        flags=re.DOTALL
    )

    # Parse each array block into Python lists
    parsed_arrays = []
    for block in array_blocks[:3]:   # Only return the first 3 arrays
        # Split on commas
        items = [x.strip() for x in block.split(',') if x.strip()]

        # Convert numeric strings to floats when possible
        cleaned = []
        for item in items:
            try:
                cleaned.append(float(item))
            except ValueError:
                cleaned.append(item.strip("'"))

        parsed_arrays.append(cleaned)

    # Ensure exactly 3 arrays returned
    while len(parsed_arrays) < 3:
        parsed_arrays.append([])

    return tuple(parsed_arrays)

def extract_idl_values(text):
    """
    Extract the VALUES of the first three IDL variable assignments.
    Returns a list of variables.
    """

    # Remove comments
    text = re.sub(r';.*', '', text)

    # Remove IDL line continuation "$\"
    text = text.replace('$\\', '')

    # Find variable assignments:  var = ...
    variable_blocks = re.findall(
        r'=\s*',
        text,
        flags=re.DOTALL
    )

    nblocks = len( variable_blocks)
    for i, ablock in enumerate( variable_blocks):
        print("%3d: %s" % (i, ablock))
        
    return variable_blocks
    #end of extract idl values
    
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
    with open(filename) as f:
        text = f.read()

    freqs, labels, weights = extract_idl_arrays(text)

    print("freqs  :", freqs)
    print("labels :", labels)
    print("weights:", weights)

    blocks = extract_idl_values( text)
       
