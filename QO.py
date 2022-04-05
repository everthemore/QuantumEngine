class QO:
    def __init__(self, states):
        self.states = states
        self.num_qubits = roof(log_2(states)) # ?

        self.state_to_index = {}
        self.index_to_state = {}
        for i in range(len(states)):
            index = binary(i)
            self.state_to_index[states[i]] = index
            self.index_to_state[index] = i

        # Default initial state
        self.state = np.zeros(2**self.num_qubits)
        self.state[0] = 1

    @property
    def states(self):
        return list(self.states)

    def peek(self, num_samples=10):
        indices = np.random.choice(range(len(states)), p=np.abs(self.state)**2)
        return [self.index_to_state(i) for i in indices]

    def pop(self):
        return self.peek(num_samples=1)

    def __add__(self, state):
        self.state[self.state_to_index[state]] += 1 # ?
        # Normalize state

#-----
myQO = QO(['0','1']) # state = '0' ([1, 0])
myQO += '1' # state = 1/sqrt(2) [1 1]
myQO += '1' # state = 1/norm [1 2]
#-----

#-----
myQO = QO(['00', '01', '10', '11'])
myQO += '11' # state is '00 + 11'
#-----

def entangle(QO1, QO2, state1, state2):
    # Assert both QO1 and QO2 are QOs

    newQO = QO(QO1.states + QO2.states)
    newQO = state1
    newQO += state2

A = QO(['0','1'])
B = QO(['0','1'])
C = entangle(A, B, '00', '11')

#------------------
#  TiqTaqToe?
squares = [QO(['E', 'X', 'O']) for i in range(9)]
pairs = []
def superpose(square1, square2, player):
    # Assert that square1 and square2 are empty

    newPair = entangle(squares[square1], squares[square2], 'E' + player, player + 'E')
    pairs.append(newPair)
