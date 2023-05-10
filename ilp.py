import numpy as np

np.set_printoptions(precision=3)


def readInput(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    # Read the number of variables and constraints
    n, m = (int(x) for x in lines[0].split())

    # Read b vector
    b = [int(x) for x in lines[1].split()]
    b = np.array(b)

    # Read c vector
    c = [int(x) for x in lines[2].split()]
    c = -np.array(c)

    # Read A matrix
    A = []
    for i in range(3, 3 + m):
        A.append([int(x) for x in lines[i].split()])
    A = np.array(A)

    return n, m, b, c, A


def initialTableau(n, m, b, c, A):
    # Create the initial tableau
    T = np.zeros((m + 1, n + m + 1))
    T[0, 1 : n + 1] = c
    T[1:, 1 : n + 1] = A
    T[1:, n + 1 : n + m + 1] = np.eye(m)
    T[1:, 0] = b
    return T


def pivot_row(T, i, j):
    # i is the row index of the pivot element
    # j is the column index of the pivot element
    # use the pivot row to eliminate all other elements in the pivot column

    for k in range(T.shape[0]):
        if k != i:
            factor = T[k, j] / T[i, j]
            T[k, :] = T[k, :] - factor * T[i, :]

    T[i, :] = T[i, :] / T[i, j]


def simplex(T):
    while np.nonzero(T[0, 1:] < 0)[0].size > 0:
        j = 1 + np.nonzero(T[0, 1:] < 0)[0][0]
        u = T[1:, j]
        positive_u = np.nonzero(u > 0)[0]
        if positive_u.size == 0:
            # print("The problem is unbounded")
            return -1

        # Find the minimum ratio
        # avoid division by zero

        ratios = np.zeros(u.size)
        for i in range(u.size):
            if u[i] > 0:
                ratios[i] = T[i + 1, 0] / u[i]
            else:
                ratios[i] = np.inf

        l = 1 + np.argmin(ratios)
        print(
            "pivot row: ",
            l,
            "pivot column: ",
            j,
            "pivot element: ",
            round(T[l, j], 3),
            "\n",
        )
        pivot_row(T, l, j)
        print(T)

    return 0

def dualSimplex(T):
    while(np.nonzero(T[1:,0]<0)[0].size>0):
        i=1+np.nonzero(T[1:,0]<0)[0][0]
        u=T[i,1:]
        negative_u=np.nonzero(u<0)[0]
        if negative_u.size==0:
            print("The problem is unbounded")
            return -1
        ratios=np.zeros(u.size)
        for j in range(u.size):
            if u[j]<0:
                ratios[j]=T[0,j+1]/u[j]
            else:
                ratios[j]=-np.inf
        k=1+np.argmax(ratios)
        print("pivot row: ",i,"pivot column: ",k)
        pivot_row(T,i,k)
        print(T)


def generate_solution(T):
    m = T.shape[0] - 1
    n = T.shape[1] - T.shape[0]

    x = np.zeros(n + m)
    I = np.eye(m)

    for i in range(m):
        for j in range(1, n + m + 1):
            if np.array_equal(I[:, i], T[1:, j]):
                x[j - 1] = T[i + 1, 0]
                break

    return x


def check_integer(T):
    sol = generate_solution(T)
    # check if every value of sol is an integer or not
    floor_sol = np.floor(sol)
    if np.isclose(sol, floor_sol).all():
        return True, floor_sol
    return False, None

def gomoryCut(T):
    m = T.shape[0] - 1
    n = T.shape[1] - T.shape[0]
    for i in range(1,m+1):
        temp=np.zeros(n+m+1)
        for j in range(n+m+1):
            temp[j]=np.floor(T[i,j])-T[i,j]
        if(np.nonzero(temp<0)[0].size>0):
            T=np.vstack((T,temp))
            T=np.hstack(T,np.zeros((T.shape[0],1)))
            T[-1,-1]=1
            return T
    

def gomory(filename):
    n, m, b, c, A = readInput(filename)
    T = initialTableau(n, m, b, c, A)
    print(T)
    status = simplex(T)
    if status == -1:
        print("The problem is unbounded")
        return

    # Check if the solution is integer
    is_integer, sol = check_integer(T)
    if is_integer:
        print("The solution is integer")
        print("The solution is: ", sol)
        return
    
    T=gomoryCut(T)


if __name__ == "__main__":
    gomory("Inputs/c.txt")
