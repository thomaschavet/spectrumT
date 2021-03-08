#Generate a randon variable with a constant probability function
import random

class randVar:
    def __init__(self, mean, bound):
        
        self.min = mean - mean/100* bound
        self.max = mean + mean/100* bound
    
    def rand(self):
        return random.uniform(self.min, self.max)
    
    def mean(self):
        return (self.min + self.max)/2