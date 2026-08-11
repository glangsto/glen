#python
#HISTORY
#26Aug10 GIL put digits in substring in {}
#26Feb15 GIL coplilot convert idl to matplottext

import re

def idl2mathtext(label):
    """
    Convert IDL-style !D...!N subscripts into Matplotlib MathText.
    Example: 'H!U13!NC!D7!NN 35-34' → r'HC$_7$N 35-34'
    """

    # Replace !Dsub!N with $_sub$
    def replDown(match):
        sub = match.group(1)
        sub = "{%s}" % sub
        return f'$_{sub}$'

    # Replace !Usub!N with $^$
    def replUp(match):
        sub = match.group(1)
        sub = "{%s}" % sub
        return f'$^{sub}$'

#    print("IDL :%s" % (label))
    # Convert all !D...!N patterns
    label = re.sub(r'!D(.*?)!N', replDown, label)
    label = re.sub(r'!U(.*?)!N', replUp, label)
#    print("Math:%s" % (label))

    return label



if __name__ == "__main__":

    import matplotlib
    import matplotlib.pyplot as plt

    label = "H!U13!NC!D7!NN 35-34"
    converted = idl2mathtext(label)

    print(converted)
    plt.title(rf'{converted}')

    plt.text(0.5, 0.5, rf'{converted}', fontsize=14)
    plt.show()
    
