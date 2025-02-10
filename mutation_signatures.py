import math
import numpy as np
signatures = {}
mutations = {}
with open('COSMIC_v3.4_SBS_GRCh37.txt') as c:
    data = c.read()
lines = [s.strip().split() for s in data.splitlines()]
for a in range(1, len(lines[0])):
    signatures[lines[0][a]] = []
for b in range(1, len(lines)):
    sum = 0
    for c in range(1, len(lines[b])):
        signatures[lines[0][c]] = c
        sum = sum + float(lines[b][c])
        key = (lines[b][0][0] + lines[b][0][2] + lines[b][0][6], lines[b][0][4])
        if key not in mutations:
            mutations[key] = [lines[b][c]]
        else:
            mutations[key].append(lines[b][c])
    print(sum)



#for i in range(0, len(mutations[('ACA', 'A')])):
for k in signatures.keys():
    signature_vector = {}
    for key in mutations.keys():
        signature_vector[key] = mutations[key][signatures[k] - 1]
    signatures[k] = signature_vector

    

def GetData(filename):
    frequencies = {}
    with open(filename) as f:
        data = f.read()
    lines = [s.strip().split() for s in data.splitlines()]
    lines.pop(0)
    lines.pop(0)
    for line in lines:
            frequencies[line[0], line[1]]= line[5]
    
    return frequencies


def DotProduct(vector1, vector2):
    dot_product = 0
    for key in vector2:
        if key in vector1:
            dot_product +=  (float(vector1[key]) * float(vector2[key]))
    return dot_product



def VectorMagnitude(vector):
    magnitude = 0
    for key in vector:
        magnitude = magnitude + (float(vector[key]) ** 2)
    magnitude = math.sqrt(magnitude)
    return magnitude


def FindAngle(dot_product, magnitude1, magnitude2):
    result = dot_product / (magnitude1 * magnitude2)
    return result


vector1 = GetData('clusters_trinucsignatures.txt')
vector1b = GetData('clusters16_trinucsignatures.txt')
vector1c = GetData('clusters16_UVB_trinucsignatures.txt')
vector1d = GetData('clusters_UVB_trinucsignatures.txt')

similarities = {}
for key in signatures.keys():
    vector2 = signatures[key]
    dot_product = DotProduct(vector1d, vector2)
    magnitude1 = VectorMagnitude(vector1d)
    magnitude2 = VectorMagnitude(vector2)
    cosine = FindAngle(dot_product, magnitude1, magnitude2)
    similarities[key] = cosine
    print(key)
    print(cosine)

vec1 = []
for key in vector2.keys():
    if key not in vector1.keys():
        vec1.append(0)
    else:
        vec1.append(float(vector1[key]))
vector1 = np.array(vec1)
test_similarities = []
for key in signatures.keys():
    vector2 = signatures[key]
    vec2 = []
    for key in vector2.keys():
        vec2.append(float(vector2[key]))
    vector2 = np.array(vec2)
    norm_vec1 = np.linalg.norm(vector1)
    norm_vec2 = np.linalg.norm(vector2)
    dot_product = np.dot(vector1, vector2)
    cosine = dot_product / (norm_vec1 * norm_vec2)
    test_similarities.append(cosine)