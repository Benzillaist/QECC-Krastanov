import numpy as np
from scipy.linalg import circulant


# Class for constructing any CSS codes given matrices
class CSS():
    def __init__(self, H, G=None):
        self.H = H
        self.G = G
        if self.G is None:
            self.G = H

        self.H_height = H.shape[0]
        self.H_width = H.shape[1]

        if H.shape == self.G.shape:
            self.PCM = np.zeros((2 * self.H_height, 2 * self.H_width), dtype=int)
            self.PCM[0:self.H_height, 0:self.H_width] = H
            self.PCM[self.H_height:2 * self.H_height, self.H_width:2 * self.H_width] = self.G

            self.X = self.PCM[:2 * self.H_height, :self.H_width]
            self.Z = self.PCM[:2 * self.H_height, self.H_width:2 * self.H_width]
        else:
            print("Error, H and G matrix widths and heights need to be the same")

    def getPCM(self):
        return self.PCM

    def checkTP(self):
        HG_product_sum = (np.matmul(self.X, self.Z.T) + np.matmul(self.Z, self.X.T)) % 2
        if np.count_nonzero(HG_product_sum) != 0:
            print("Error, invalid H and/or G, twisted product condition not satisfied")
            print(HG_product_sum)
        else:
            print("Valid code")


class BStarter:
    def __init__(self, c_indices, n):
        c_indices = np.array(c_indices)
        if not (c_indices < n).all():
            print("Error, a index is greater or equal to the length of the list")
        else:
            self.c_start = np.zeros(n, dtype=int)
            self.c_start[c_indices] += 1
            self.circulantT = np.array(circulant(self.c_start))
            self.comp = np.zeros((n, 2 * n), dtype=int)
            self.comp[:, n:2 * n] = self.circulantT
            self.comp[:, 0:n] = self.circulantT.T
            self.H_0 = self.comp[0:n, 0:2 * n]

    def getH(self):
        return self.H_0


class UStarter:
    def __init__(self, c_indices, n):
        c_indices = np.array(c_indices)
        if not (c_indices < n).all():
            print("Error, a index is greater or equal to the length of the list")
        else:
            print("2")
