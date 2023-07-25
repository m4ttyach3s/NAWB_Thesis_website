from flask import Flask, render_template, url_for, request, redirect
import pandas as pd
import json
import RNA
import datetime
import os

allowedBases = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']
nucleotides = ['A', 'C', 'U', 'G']
bonds = []

app = Flask(__name__)

@app.route('/', methods=['POST', 'GET'])
def index():
    if request.method == 'POST':
        mll_value = request.form['mll']  # Get the minimal-loop value from the form submission
        content_value = request.form['content']
        try:
            rnaMatrix, db_notation, gen_rna = NussinovAlgorithm(content_value, mll_value)
        except ValueError as e:
            error_message = str(e)
            return render_template('index.html', error_message=error_message)
            # If an error is raised, store the error message
        
        print(gen_rna)
        rnaMatrix_json = json.dumps(rnaMatrix)
        global bonds
        bonds = []
        return render_template('result.html', mll_value=mll_value, content_value=content_value, rnaMatrix_json=rnaMatrix_json, db_notation=db_notation, gen_rna=gen_rna)
    else:
        return render_template('index.html')
    #return render_template('index.html')

def NussinovAlgorithm(s, L):
    """Calls the functions of the Nussinov Algorithm. Prints the number of bonds that are found.

    Args:
        s (str): Input RNA string.
    Returns:
        None: None.
    """
    s = stringCheck(s)
    rnaMatrix = matinit(s)
    nussinovScore(s, L, rnaMatrix)
    traceback(rnaMatrix, 0, len(s)-1, s)
    bonds.sort()
    db_notation = print_dotbracket(s)
    gen_rna = foldRNA(s, db_notation)
    return rnaMatrix, db_notation, gen_rna
    #print(f"\nThere are {len(bonds)} bonds. They are formed in these positions: {bonds}")
    #db_notation = print_dotbracket(s)
    #foldRNA(s, db_notation)
    
def stringCheck(s): # tested
    
    """Checks if the RNA string is a permitted sequence of characters (A, C, G, U) of length greater than 5.

    Args:
        s (str): input RNA string

    Raises:
        ValueError: The input string's length is not correct
        ValueError: The input string has characters that are not allowed

    Returns:
        str: formatted string to upper case
    """
    
    s = str(s) #necessary type-casting
    s = s.upper()

    if len(s) < 6:
        raise ValueError(f"The input {s} given is not correct. It should be of at least 6 characters.")
    
    if not all(base in nucleotides for base in s):
        raise ValueError(f"The input {s} has some characters that are not allowed. Please, only use characters such as A, C, U, G.")
    
    return s

def matinit(s):
    
    """Initializes the matrix with "-1" values. The 1st and 2nd diagonal are initalized with "0" values. Only native data structures are used, such as lists.

    Args:
        s (str): Input RNA string

    Returns:
        list: list that represents a matrix with the initialized values.
    """
    
    # I am going to initialize the list with list comprehensions
    rnaMatrix = [[-1] * len(s) for i in range(len(s))]
    
    for i in range(0, len(s)):
        rnaMatrix[i][i] = 0 # 1st diagonal
        
    for i in range(1, len(s)):
        rnaMatrix[i][i-1] = 0 # 2nd diagonal
        
    return rnaMatrix

def nussinovScore(s, L, rnaMatrix):
    
    """Calculates the score of the RNA matrix, based on the base pairing, the minimal loop length provided and the maximum value found. It uses the Nussinov algorithm DP-matrix score function provided by the book "Biological Sequence Analysis" by R. Durbin.

    Args:
        s (str): Input RNA string
        rnaMatrix (list): Matrix to check of the maximum value.
    Returns:
        None: the matrix is filled with the values.
    """
    L = int(L)
    
    for t in range(1, len(s)):
        for i in range(0, len(s)-t):
            j = i+t # width of the diagonal
                
            case1 = rnaMatrix[i+1][j] # unpaired position i
            case2 = rnaMatrix[i][j-1] # unpaired position j
            
            if i+L<j and complementary(s[i], s[j]):
                case3 = rnaMatrix[i+1][j-1] + 1 #!!!! here i have to check L, because it's the only case with complementary assignment!
                # adds the pair!
            else:
                case3 = 0
                
            case4 = 0
            for k in range(i, j-1):
                case4 = max(case4, rnaMatrix[i][k] + rnaMatrix[k+1][j]) # I don't need an extra variable for the max, since I would assign it to case4 in any case
            # here we combine the substructures
            rnaMatrix[i][j] = max(case1, case2, case3, case4)
    
def complementary(baseOne, baseTwo): # tested
    
    """Checks if the input character can form a Watson-Crick base pair (A-U, C-G) or a Wobble base pair(G-U).
    The alphabet is limited to the characters A, C, G, U.

    Args:
        baseOne (str): first character (or base)
        baseTwo (str): second character (or base)

    Returns:
        bool: True, if they pair; False, otherwise.
    """
    return baseOne+baseTwo in allowedBases

def traceback(rnaMatrix, i, j, s): #tested
    
    """Implements recursively the traceback function to find the best structure for the secondary structure as described by the book "Biological Sequence Analysis by R. Durbin.

    Args:
        rnaMatrix (list): Calculated score matrix for the RNA base pairs.
        i (int): First position of the string.
        j (int): Second position of the string.
        s (str): Input RNA string.

    Returns:
        function: Function with updated values.
    """

    if i >= j: 
        return 0  
    
    if rnaMatrix[i+1][j] == rnaMatrix[i][j]:
        traceback(rnaMatrix, i+1, j, s)
    elif rnaMatrix[i][j-1] == rnaMatrix[i][j]:
        traceback(rnaMatrix, i, j-1, s)
    elif rnaMatrix[i+1][j-1] + complementary(s[i],s[j]) == rnaMatrix[i][j]:
        bonds.append((i, j)) #record base pair i,j
        traceback(rnaMatrix, i+1, j-1, s)
    else:
        for k in range(i+1, j-1): # two tracebacks are needed to reach the end
            if rnaMatrix[i][k] + rnaMatrix[k+1][j] == rnaMatrix[i][j]:
                traceback(rnaMatrix, k+1, j, s)
                traceback(rnaMatrix, i, k, s)
                break

def print_dotbracket(s):
    
    """Calculates and prints the dot-bracket notation, given the input string. The dot-bracket notation only allows the characters "(", ")" and ".", where
    "(" reprents the first base in the pairing and ")" its counterpart to close the pair. "." is a base that does not form any pairing.

    Args:
        s (str): Input string.

    Returns:
        str: formatted dot-bracket notation.
    """
    
    n = len(s)
    dotbracket = ["."] * n
        
    for j in bonds:
        dotbracket[j[0]] = "("
        dotbracket[j[1]] = ")"
    return ''.join(dotbracket)

def foldRNA(s, db_notation):
    """Generates the graphic representation of the given RNA input string.

    Args:
        s (str): Input RNA string
        db_notation (str): Dot-bracket notation string.
        
    Returns:
        string: The function generates an SVG file representing the RNA folding returns the filename.
    """
    current_datetime = datetime.datetime.now()
    formatted_datetime = current_datetime.strftime("%H_%M_%S_%Y_%m_%d")
    file_name = "RNAfold_timestamp"+formatted_datetime+".svg"
    RNA.svg_rna_plot(s, db_notation, "static/images/"+file_name)
    
    return file_name

if __name__ == "__main__":
    app.run(debug=True, port=3000)