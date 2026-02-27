import numpy as np
import matplotlib.pyplot as plt

def main():
    filename = 'output.txt'
    
    # Load the distribution file 
    data = np.loadtxt(filename)
    
    print(f"Successfully loaded {len(data)} numbers")

    # Create the histogram
    plt.hist(data, bins=100, density=True, alpha=0.7, color='skyblue', edgecolor='black')
    plt.title('Plot of Generated Normal Random Variables')
    plt.xlabel('Value (x)')
    plt.ylabel('Probability Density (Phi(x))')
    plt.grid(axis='y', alpha=0.75)

    plt.show()

if __name__ == "__main__":
    main()