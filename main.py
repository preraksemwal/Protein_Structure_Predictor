# AUTHOR            :  Prerak Semwal
# Course Instructor :  Prof. Jaspreet Kaur Dhanjal


protein_sequence = "WHGCITVYWMTV"

# alpha_propensity = {
#     "G" : 0.53,
#     "A" : 1.45,
#     "V" : 1.14,
#     "L" : 1.34,
#     "I" : 1.00,
#     "R" : 0.79,
#     "K" : 1.07,
#     "E" : 1.53,
#     "D" : 0.98,
#     "Q" : 1.17,
#     "N" : 0.73,
#     "T" : 0.82,
#     "S" : 0.79,
#     "C" : 0.77,
#     "M" : 1.20,
#     "F" : 1.12,
#     "Y" : 0.61,
#     "W" : 1.14,
#     "H" : 1.24,
#     "P" : 0.59
# }
    
# beta_propensity = {
#     "G" : 0.81,
#     "A" : 0.97,
#     "V" : 1.65,
#     "L" : 1.22,
#     "I" : 1.60,
#     "R" : 0.90,
#     "K" : 0.74,
#     "E" : 0.26,
#     "D" : 0.80,
#     "Q" : 1.23,
#     "N" : 0.65,
#     "T" : 1.20,
#     "S" : 0.72,
#     "C" : 1.30,
#     "M" : 1.67,
#     "F" : 1.28,
#     "Y" : 1.29,
#     "W" : 1.19,
#     "H" : 0.71,
#     "P" : 0.62
# }

##########################################################   delete whats between
alpha_propensity = {
    "G" : 0.53,
    "A" : 1.45,
    "V" : 1.14,
    "L" : 1.34,
    "I" : 1.00,
    "R" : 0.79,
    "K" : 1.07,
    "E" : 1.53,
    "D" : 0.98,
    "Q" : 1.17,
    "N" : 0.73,
    "T" : 0.82,
    "S" : 0.79,
    "C" : 0.77,
    "M" : 1.20,
    "F" : 1.12,
    "Y" : 0.61,
    "W" : 1.14,
    "H" : 1.24,
    "P" : 0.59
}
    
beta_propensity = {
    "G" : 0.81,
    "A" : 0.97,
    "V" : 1.65,
    "L" : 1.22,
    "I" : 1.60,
    "R" : 0.90,
    "K" : 0.74,
    "E" : 0.26,
    "D" : 0.80,
    "Q" : 1.23,
    "N" : 0.65,
    "T" : 1.20,
    "S" : 0.72,
    "C" : 1.30,
    "M" : 1.17,
    "F" : 1.28,
    "Y" : 1.29,
    "W" : 1.19,
    "H" : 0.71,
    "P" : 0.62
}
##########################################################



prediction   = []    # stores the final structure predicted
found_alphas = []
found_betas  = []


def isAlphaCandidate(seq):  # 'seq' is sequence of length 6 and check if its a contender for alpha-helical structure
    count = 0
    for e in seq:
        if alpha_propensity[e] >= 1:
            count += 1
    if count >= 4 :
        print(seq , " " , "selected")
        return True
    return False

def isBetaCandidate(seq):  # 'seq' is sequence of length 5 and check if its a contender for beta-helical structure
    count = 0
    for e in seq:
        if beta_propensity[e] > 1:
            count += 1
    if count >= 3:
        print(seq , " " , "selected")
        return True
    return False

def expand(i,j):
    if(j-i+1 == 6):     # alpha_helix case
        pass
        # expand leftwards
        
    else:
        pass
    


if __name__ == '__main__':
    
    prediction   = [" "] * len(protein_sequence)
    found_betas  = [" "] * len(protein_sequence)
    found_alphas = [" "] * len(protein_sequence)


    # detecting for alpha-helix
    for i in range(0, len(protein_sequence) - 5):
        seq = protein_sequence[i : i+6]
        if isAlphaCandidate(seq):
            for j in range(i, i+6):
                found_alphas[j] = "H"
            expand(i, i+5)

    print(found_alphas)
    



    # detecting for beta-strand
    for i in range(0, len(protein_sequence) - 4):
        seq = protein_sequence[i : i+5]
        if isBetaCandidate(seq):
            for j in range(i, i+5):
                found_betas[j] = "E"
            # expand(i, i+5)

