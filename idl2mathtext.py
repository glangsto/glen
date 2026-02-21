#python
#HISTORY
#26Feb15 GIL coplilot convert idl to matplottext

import re

def idl2mathtext(label):
    """
    Convert IDL-style !D...!N subscripts into Matplotlib MathText.
    Example: 'HC!D7!NN 35-34' → r'HC$_7$N 35-34'
    """

    # Replace !Dsub!N with $_sub$
    def repl(match):
        sub = match.group(1)
        return f'$_{sub}$'

    # Convert all !D...!N patterns
    label = re.sub(r'!D(.*?)!N', repl, label)

    return label



if __name__ == "__main__":

    import matplotlib
    import matplotlib.pyplot as plt

    label = "HC!D7!NN 35-34"
    converted = idl2mathtext(label)

    print(converted)
    plt.title(rf'{converted}')

    plt.text(0.5, 0.5, rf'{converted}', fontsize=14)
    plt.show()
    
