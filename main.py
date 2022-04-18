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
found_alphas = []    # stores info of each amino acid about whether ALPHA-helix possible or not
found_betas  = []    # stores info of each amino acid about whether BETA-strand possible or not


def isAlphaCandidate(seq):  # 'seq' is sequence of length 6 and check if its a contender for alpha-helical structure
    if len(seq) == 4:   # during expansion of a selected sequence
        score = 0
        for e in seq:
            score += alpha_propensity[e]
        if score >= 4:
            return True
        return False
    else:
        count = 0
        for e in seq:
            if alpha_propensity[e] >= 1:
                count += 1
        if count >= 4 :
            return True
        return False

def isBetaCandidate(seq):  # 'seq' is sequence of length 5 and check if its a contender for beta-helical structure
    if len(seq) == 4:   # during expansion of a selected sequence
        score = 0
        for e in seq:
            score += beta_propensity[e]
        if score > 4:
            return True
        return False
    else:
        count = 0
        for e in seq:
            if beta_propensity[e] > 1:
                count += 1
        if count >= 3:
            return True
        return False

def expand(i,j):
    if j-i+1 == 6:     # ALPHA HELIX - CASE
        # expand leftwards
        ptr1 = i-1
        while ptr1 >=0 and isAlphaCandidate(protein_sequence[ptr1 : ptr1+4]):
            found_alphas[ptr1] = "H"
            ptr1 -= 1

        # expand rightwards
        ptr2 = j+1
        while ptr2 <= len(protein_sequence)-1 and isAlphaCandidate  (protein_sequence[ptr2-3 : ptr2+1]):
            found_alphas[ptr2] = "H"
            ptr2 += 1
        
    else:               # BETA STRAND - CASE
        # expand leftwards
        ptr1 = i-1
        while ptr1 >=0 and isBetaCandidate(protein_sequence[ptr1 : ptr1+4]):
            found_betas[ptr1] = "S"
            ptr1 -= 1

        # expand rightwards
        ptr2 = j+1
        while ptr2 <= len(protein_sequence)-1 and isBetaCandidate(protein_sequence[ptr2-3 : ptr2+1]):
            found_betas[ptr2] = "S"
            ptr2 += 1
    


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

    # detecting for beta-strand
    for i in range(0, len(protein_sequence) - 4):
        seq = protein_sequence[i : i+5]
        if isBetaCandidate(seq):
            for j in range(i, i+5):
                found_betas[j] = "S"
            expand(i, i+5)

    i = 0
    while i < len(prediction):
        if found_alphas[i] == "H" and found_betas[i] == " ":
            prediction[i] = "H"
            i += 1
        elif found_alphas[i] == " " and found_betas[i] == "S":
            prediction[i] = "S"
            i += 1
        else:
            upto = i
            while upto < len(prediction)-1 and found_alphas[upto] == "H" and found_betas[upto] == "S":
                upto += 1
            if upto > i:    # if upto is != 1 , then upto-th index WON'T have a conflict
                upto -= 1
            alpha_score, beta_score = 0, 0
            for e in protein_sequence[i : upto+1]:
                alpha_score += alpha_propensity[e]
                beta_score += beta_propensity[e]

            if alpha_score > beta_score:
                for j in range(i, upto + 1):
                    prediction[j] = "H"
            else:
                for j in range(i, upto + 1):
                    prediction[j] = "S"
            if upto == i:
                i += 1
            else:
                i = upto


    print("Structure Prediction as Per Chou-Fasman :  ", end = "")
    for e in prediction:
        print(e , end = "")
