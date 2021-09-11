#Generate a randon variable from Gaussian or constant probability function
import random
import numpy as np

class randVar:
    def __init__(self, mean, bound):#bound in percent
        
        self.min = mean - mean/100* bound
        self.max = mean + mean/100* bound
    
    def rand(self):
        return random.uniform(self.min, self.max)
    
    def randgauss(self):#bound correspond to 2*std_dev
        return np.random.normal(self.mean(), (self.max-self.mean())/2)
    
    def mean(self):
        return (self.min + self.max)/2